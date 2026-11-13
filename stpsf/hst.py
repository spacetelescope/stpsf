"""
===============
HST Instruments
===============

A re-implementation of the TinyTim optical model for the Hubble Space Telescope, built on poppy and
the STPSF instrument framework. Currently this covers WFC3 (IR and UVIS channels).

WARNING: Work in progress. The goal of this first phase is to reproduce TinyTim's outputs, *not* to improve on
them. Where TinyTim does something questionable, we do the same by default and flag it below.

Optical model
-------------
Each optical element is a separate poppy plane, in this order:

1. ``HSTPupil``: analytic OTA pupil (R = 1.2 m, central obscuration 0.330 R, 4 spiders of width 0.022 R along V2/V3,
   3 mirror support pads). Pixels are antialiased by 4x4 subsampling.
2. ``HSTMirrorMap`` for the primary mirror: TinyTim's measured zonal OPD map.
3. ``HSTMirrorMap`` for the secondary mirror, shifted with field position (the beam footprint moves on the secondary).
4. ``HSTAnnularZernikeWFE``: instrument aberrations, from field-dependent polynomials for Z4-Z8 (plus despace).
5. For WFC3 IR only, ``WFC3IRColdMask`` (a subclass of ``HSTPupil``): the cold pupil mask (outer 0.9755 R, secondary 0.354 R,
   spiders 0.061 R). It is applied in addition to the OTA pupil.

Everything is evaluated in TinyTim's frame: the OTA V2/V3 frame with +V3 along array columns (poppy x) and +V2
along rows (poppy y), rotated clockwise by 135 degrees for both WFC3 channels. The pupil, mirror maps and Zernikes
are all evaluated through the same coordinate transform (``_ota_frame``), so they are always mutually registered;
no poppy rotation plane is used. TinyTim's FFT and phase sign conventions are identical to poppy's, so no flips
are needed.

Data files (not yet in the released STPSF data package)
------------------------------------------------------
Generated from a TinyTim installation by ``dev_utils/make_hst_data.py``, written under ``$STPSF_PATH``:
``HST/HST_{primary,secondary}_mirror_map.fits`` (OPD in meters; the primary header gives the pupil radius in map
pixels), and ``WFC3/`` with ``filters.tsv``, ``filters/*_throughput.fits`` and ``WFC3_tinytim_model.fits``
(field polynomials, static Zernikes, pad geometry, detector kernels, per detector).

Conventions and TinyTim quirks to be aware of
---------------------------------------------
* **Units.** TinyTim Zernike coefficients are RMS waves at 547 nm (``constants.HST_TINYTIM_REFERENCE_WAVELENGTH``).
  The mirror map files here are converted to meters. TinyTim's dumped wave map is in waves at the *current*
  wavelength, not 547 nm.
* **Detector positions are 0-based** (X, Y), as in TinyTim. SIAF calls are made with 1 added.
* **Pupil V2/V3 = minus field V2/V3.** TinyTim's "pupil" coordinates are the negative of the SIAF V2, V3. This
  matters for the secondary mirror shift, which is measured relative to the WFPC2 PC position.
* **Zernike azimuth** is measured from +V2 toward +V3, so "cos" terms align with V2 and "sin" terms with V3.
  Z5 is therefore minus the usual Noll Z6, and Z6 is Noll Z5.
* **Zero-focus quirk.** TinyTim only applies a field-dependent Zernike term if the static table value of that term is
  nonzero. The WFC3 tables use 1e-12 placeholders to switch the terms on, but ``tiny1`` rewrites the Z4 entry with
  ``%f`` (to add the despace), which turns the placeholder into exactly 0. So with zero despace, the
  field-dependent *focus* (about +-9 nm RMS over the detector) is silently dropped; Z5-Z8 are kept. We reproduce this
  by default (``tinytim_skip_zero_terms = True``). With a nonzero despace the field focus is applied. Set
  ``tinytim_skip_zero_terms = False`` to always include it.
* **Half-pixel centering.** TinyTim's 512x512 pupil and wave maps are centered on pixel 256, half a pixel from the
  poppy grid center. We center everything on the poppy grid, which has no effect on PSF intensities. Edge pixels
  differ slightly from TinyTim's because the antialiasing differs.
* **Annular Zernikes use epsilon = 0.33 even for the IR cold mask** (whose obscuration is 0.354), as TinyTim does.
* **Pixel scale.** SIAF HST apertures do not define ``XSciScale``; the mean scale is measured from the SIAF
  transform. TinyTim's nominal scales (0.13 and 0.04 arcsec) are available with ``use_tinytim_pixelscale``.

Validation
----------
The wavefront (pupil, mirror maps, Zernikes, secondary shift) was checked against TinyTim's dumped
``tinytim_wave.fits`` at two IR field positions, agreeing to better than 0.05 nm RMS. See
``stpsf/tests/test_hst.py`` for the tests that run in CI, and the TinyTim notes
(``notes/03_wfc3_ir_porting_guide.md``) for the full validation plan.
"""

import logging
import os
from functools import lru_cache

import astropy.io.fits as fits
import astropy.units as u
import numpy as np
import poppy
from astropy.io import ascii as ioascii
from scipy import ndimage

from . import constants, stpsf_core, utils
from .distortion import _siaf_mean_pixelscale

_log = logging.getLogger('stpsf')


def _in_spiders(v2, v3, half_width):
    """True where a point lies on one of the 4 spiders, which run along the V2 and V3 axes"""
    return (np.abs(v2) < half_width) | (np.abs(v3) < half_width)


class HSTPupil(poppy.AnalyticOpticalElement):
    """Analytic HST OTA entrance pupil in TinyTim's geometry.

    The geometry is defined in the OTA V2/V3 frame, with +V3 along the array columns (poppy x) and
    +V2 along the array rows (poppy y), and is then rotated clockwise by ``rotation`` degrees to
    the camera frame, following TinyTim's convention. Subclasses (e.g. `WFC3IRColdMask`) override
    `_transmission_at` to model other binary masks in the same frame.

    Parameters
    ----------
    radius : float
        Outer radius of the OTA pupil in meters.
    rotation : float
        Clockwise rotation in degrees.
    subsample : int
        Each pixel is evaluated on a subsample x subsample grid of points and averaged.
    """

    default_name = 'HST OTA pupil'

    def __init__(self, radius=constants.HST_PUPIL_DIAMETER / 2, rotation=constants.HST_TINYTIM_WFC3_ROTATION,
                 subsample=4, name=None, **kwargs):
        super().__init__(name=name or self.default_name, planetype=poppy.poppy_core.PlaneType.pupil, **kwargs)
        self.pupil_diam = 2 * radius
        self.radius = radius
        self.rotation_cw = rotation
        self.subsample = subsample

    def _transmission_at(self, v2, v3):
        """Binary transmission at points given in units of the pupil radius, in the OTA V2/V3 frame"""
        r = np.hypot(v2, v3)
        trans = (r <= 1) & (r >= constants.HST_OTA_SECONDARY_RADIUS)
        trans &= ~_in_spiders(v2, v3, constants.HST_OTA_SPIDER_WIDTH / 2)
        for pad_v3, pad_v2, pad_r in constants.HST_OTA_PADS:
            trans &= np.hypot(v3 - pad_v3, v2 - pad_v2) > pad_r
        return trans

    def get_transmission(self, wave):
        y, x = self.get_coordinates(wave)
        n = self.subsample
        pix = wave.pixelscale.to_value(u.m / u.pixel) if hasattr(wave.pixelscale, 'unit') else wave.pixelscale
        offsets = (np.arange(n) + 0.5) / n - 0.5
        total = np.zeros(x.shape)
        for dy in offsets:
            for dx in offsets:
                v2, v3 = _ota_frame(x + dx * pix, y + dy * pix, self.radius, self.rotation_cw)
                total += self._transmission_at(v2, v3)
        return total / n**2


class WFC3IRColdMask(HSTPupil):
    """The WFC3 IR channel's cold pupil mask, a physical element separate from the OTA pupil.

    It sits in the WFC3 IR optical path, undersizing the OTA pupil slightly (outer edge, a larger secondary
    obscuration and wider spiders), and is applied in addition to `HSTPupil`. It is registered with the OTA pupil
    through the shared V2/V3 frame and rotation. Dimensions are in ``constants.HST_WFC3_IR_COLD_MASK``; there are no
    support pads in the mask.
    """

    default_name = 'WFC3 IR cold mask'

    def _transmission_at(self, v2, v3):
        cm = constants.HST_WFC3_IR_COLD_MASK
        r = np.hypot(v2, v3)
        return (r <= cm['outer_radius']) & (r >= cm['secondary_radius']) & ~_in_spiders(v2, v3, cm['spider_width'] / 2)


def _ota_frame(x, y, radius, rotation):
    """Convert pupil-plane coordinates (x, y) in meters to TinyTim's unrotated OTA (V2, V3) frame, in units of
    the pupil radius. The camera frame is rotated clockwise by ``rotation`` degrees relative to that frame."""
    theta = np.deg2rad(rotation)
    xs, ys = x / radius, y / radius
    return xs * np.sin(theta) + ys * np.cos(theta), xs * np.cos(theta) - ys * np.sin(theta)


class HSTMirrorMap(poppy.AnalyticOpticalElement):
    """OPD of one HST OTA mirror, from a TinyTim zonal error map, in the OTA V2/V3 frame rotated like the pupil.

    Parameters
    ----------
    opd_map : ndarray
        Square map in meters of OPD; its center, at index ``n/2``, is the pupil center.
    pupil_radius_pix : float
        Radius of the HST pupil in map pixels.
    shift : tuple of float
        (V2, V3) shift, in units of the pupil radius and in the unrotated OTA frame, of the point of the map that is
        sampled at the pupil center. For the secondary mirror this is the displacement of the beam footprint,
        4906.071 mm * tan(field offset from the WFPC2 PC position) / (267.3 mm / 2), using TinyTim's "pupil"
        (minus field) V2, V3. The map contents therefore move by minus this shift on the pupil, matching TinyTim.

    Notes
    -----
    The map is bilinearly interpolated, and is zero outside its 344 x 344 pixel extent (1.27 pupil radii). Both maps
    share one pixel scale, set by the primary's pupil radius of 135.15 map pixels (from the phase retrieval used to
    derive the primary map). The OPD is added as given, i.e. these maps include the factor of two for reflection.
    """

    def __init__(self, opd_map, pupil_radius_pix, shift=(0.0, 0.0), radius=constants.HST_PUPIL_DIAMETER / 2,
                 rotation=constants.HST_TINYTIM_WFC3_ROTATION, name='HST mirror map', **kwargs):
        super().__init__(name=name, planetype=poppy.poppy_core.PlaneType.pupil, **kwargs)
        self.opd_map = np.asarray(opd_map, dtype=float)
        self.pupil_radius_pix = pupil_radius_pix
        self.shift = shift
        self.radius = radius
        self.rotation_cw = rotation

    def get_transmission(self, wave):
        return np.ones(wave.shape)

    def get_opd(self, wave):
        y, x = self.get_coordinates(wave)
        v2, v3 = _ota_frame(x, y, self.radius, self.rotation_cw)
        center = self.opd_map.shape[0] / 2
        rows = center + (v2 + self.shift[0]) * self.pupil_radius_pix
        cols = center + (v3 + self.shift[1]) * self.pupil_radius_pix
        return ndimage.map_coordinates(self.opd_map, [rows, cols], order=1, mode='constant', cval=0.0)


# Annular Zernikes from TinyTim opd.c, orthonormal on 0.33 < rho < 1. Each entry is
# (normalization, azimuthal order m, radial polynomial coefficients as {power of rho: coefficient}, trig function).
# The azimuth is measured from +V2 toward +V3.
_TT_ZERNIKES = {
    2: (1.8992573, 1, {1: 1.0}, 'cos'),
    3: (1.8992573, 1, {1: 1.0}, 'sin'),
    4: (3.8874443, 0, {2: 1.0, 0: -0.554450}, None),
    5: (2.3137662, 2, {2: 1.0}, 'cos'),
    6: (2.3137662, 2, {2: 1.0}, 'sin'),
    7: (8.3345629, 1, {3: 1.0, 1: -0.673796}, 'cos'),
    8: (8.3345629, 1, {3: 1.0, 1: -0.673796}, 'sin'),
    9: (2.6701691, 3, {3: 1.0}, 'cos'),
    10: (2.6701691, 3, {3: 1.0}, 'sin'),
    11: (16.895979, 0, {4: 1.0, 2: -1.1089, 0: 0.241243}, None),
    12: (12.033645, 2, {4: 1.0, 2: -0.750864}, 'cos'),
    13: (12.033645, 2, {4: 1.0, 2: -0.750864}, 'sin'),
    14: (2.9851527, 4, {4: 1.0}, 'cos'),
    15: (2.9851527, 4, {4: 1.0}, 'sin'),
    16: (36.321412, 1, {5: 1.0, 3: -1.230566, 1: 0.323221}, 'cos'),
    17: (36.321412, 1, {5: 1.0, 3: -1.230566, 1: 0.323221}, 'sin'),
    18: (16.372202, 3, {5: 1.0, 3: -0.8001}, 'cos'),
    19: (16.372202, 3, {5: 1.0, 3: -0.8001}, 'sin'),
    20: (3.2700486, 5, {5: 1.0}, 'cos'),
    21: (3.2700486, 5, {5: 1.0}, 'sin'),
    22: (74.82446, 0, {6: 1.0, 4: -1.663350, 2: 0.803136, 0: -0.104406}, None),
}


def _tinytim_annular_term(j, rho, theta):
    norm, m, radial, trig = _TT_ZERNIKES[j]
    value = norm * sum(c * rho**p for p, c in radial.items())
    if trig is not None:
        value = value * (np.cos(m * theta) if trig == 'cos' else np.sin(m * theta))
    return value


def tinytim_annular_zernike_basis(nterms=22, npix=512, rho=None, theta=None, outside=None):
    """Basis set of TinyTim's annular Zernikes, in the style of `poppy.zernike.zernike_basis`.

    Parameters
    ----------
    nterms : int
        Number of terms, 1 to 22. The returned list is indexed so that element ``j - 1`` is TinyTim term ``j``
        (as for the Noll indexing of `poppy.zernike.zernike_basis`). Element 0 is piston, which is 1 everywhere.
    npix : int
        Size of the square grid, used only if ``rho`` and ``theta`` are not given. The grid then spans
        -1 to 1 pupil radii.
    rho, theta : array_like, optional
        Radius in units of the pupil radius, and azimuth in radians measured from +V2 toward +V3.
    outside : float, optional
        If given, the value used where ``rho > 1`` (as in poppy's basis functions). By default terms
        are evaluated everywhere, which is what TinyTim does.

    Returns
    -------
    list of ndarray
        The terms, orthonormal on the annulus 0.33 < rho < 1.
    """
    if not 1 <= nterms <= len(_TT_ZERNIKES) + 1:
        raise ValueError(f'nterms must be between 1 and {len(_TT_ZERNIKES) + 1}')
    if rho is None or theta is None:
        v3, v2 = np.meshgrid(*[np.linspace(-1, 1, npix)] * 2)
        rho, theta = np.hypot(v2, v3), np.arctan2(v3, v2)
    rho, theta = np.asarray(rho, dtype=float), np.asarray(theta, dtype=float)
    basis = [np.ones(rho.shape)] + [_tinytim_annular_term(j, rho, theta) for j in range(2, nterms + 1)]
    if outside is not None:
        for term in basis:
            term[rho > 1] = outside
    return basis


class HSTAnnularZernikeWFE(poppy.AnalyticOpticalElement):
    """Wavefront error as a sum of TinyTim's annular Zernikes (j = 2-22), in the rotated OTA frame.

    Parameters
    ----------
    coefficients : dict
        Coefficients keyed by TinyTim term number, in waves RMS at 547 nm.

    Notes
    -----
    These are the Mahajan annular Zernikes with obscuration ratio 0.33, orthonormal on 0.33 < rho < 1, as in
    TinyTim's ``opd.c``, provided by `tinytim_annular_zernike_basis`. They are not the Noll/circular Zernikes of
    `poppy.zernike`, and the term numbering and azimuth convention differ (see the module docstring). Term 1
    (piston) is not used.
    """

    def __init__(self, coefficients, radius=constants.HST_PUPIL_DIAMETER / 2,
                 rotation=constants.HST_TINYTIM_WFC3_ROTATION, name='HST SI aberrations', **kwargs):
        super().__init__(name=name, planetype=poppy.poppy_core.PlaneType.pupil, **kwargs)
        self.coefficients = dict(coefficients)
        self.radius = radius
        self.rotation_cw = rotation

    def get_transmission(self, wave):
        return np.ones(wave.shape)

    def get_opd(self, wave):
        y, x = self.get_coordinates(wave)
        v2, v3 = _ota_frame(x, y, self.radius, self.rotation_cw)
        opd = np.zeros(wave.shape)
        used = [j for j, coeff in self.coefficients.items() if coeff != 0]
        if used:
            basis = tinytim_annular_zernike_basis(max(used), rho=np.hypot(v2, v3), theta=np.arctan2(v3, v2))
            for j in used:
                opd += self.coefficients[j] * basis[j - 1]
        return opd * constants.HST_TINYTIM_REFERENCE_WAVELENGTH


@lru_cache(maxsize=2)
def _load_hst_mirror_maps(datapath):
    """Load the HST primary and secondary mirror OPD maps (meters) and the pupil radius in map pixels"""
    with fits.open(os.path.join(datapath, 'HST', 'HST_primary_mirror_map.fits')) as hdul:
        primary = np.array(hdul[0].data, dtype=float)
        radius_pix = hdul[0].header['PUPRADPX']
    with fits.open(os.path.join(datapath, 'HST', 'HST_secondary_mirror_map.fits')) as hdul:
        secondary = np.array(hdul[0].data, dtype=float)
    return primary, secondary, radius_pix


@lru_cache(maxsize=4)
def _load_tinytim_model(path):
    """Load a ``*_tinytim_model.fits`` file (see dev_utils/make_hst_data.py) into a dict keyed by detector name.

    Each entry holds the field polynomial coefficients, static Zernikes, pad geometry and kernels, plus the
    scalar parameters from the FITS headers.
    """
    model = {}
    with fits.open(path) as hdul:
        detectors = sorted({ext.name.split('_')[0] for ext in hdul[1:]})
        for det in detectors:
            hdr = hdul[f'{det}_FIELD'].header
            entry = {
                'field_coeffs': np.array(hdul[f'{det}_FIELD'].data),
                'pads': np.array(hdul[f'{det}_PADS'].data),
                'static_zernikes': np.array(hdul[f'{det}_ZSTATIC'].data),
                'kernels': np.array(hdul[f'{det}_KERNEL'].data),
                'kernel_waves_um': np.array(
                    [hdul[f'{det}_KERNEL'].header[f'WAVE{i + 1}'] for i in range(hdul[f'{det}_KERNEL'].data.shape[0])]
                ),
                'kfield_coeffs': np.array(hdul[f'{det}_KFIELD'].data),
                'kfield_waves_um': np.array(
                    [hdul[f'{det}_KFIELD'].header[f'WAVE{i + 1}'] for i in range(hdul[f'{det}_KFIELD'].data.shape[0])]
                ),
            }
            entry.update({key.lower(): hdr[key] for key in hdr if key not in ('XTENSION', 'BITPIX', 'PCOUNT', 'GCOUNT')
                          and not key.startswith(('NAXIS', 'COMMENT', 'EXTNAME'))})
            model[det] = entry
    return model


class HSTInstrument(stpsf_core.SpaceTelescopeInstrument):
    """
    HSTInstrument contains data and functionality common to HST instruments, such as the OTA pupil geometry
    and the HST SIAF.

    Detector positions for HST instruments are **0-based** pixel indices (X, Y), matching Python indexing and
    TinyTim's input convention. SIAF science coordinates are 1-based, so 1 is added in SIAF calls.

    The telescope model (see the module docstring for details) is shared by all HST instruments: an analytic OTA
    pupil, then the primary and secondary mirror OPD maps, then instrument aberrations, all in TinyTim's rotated
    OTA frame. The JWST OTE exit-pupil inversion and rotation are not applied.

    Attributes
    ----------
    despace_microns : float
        Secondary mirror despace in microns, which adds ``0.011`` waves RMS (547 nm) of focus per micron.
        Nonzero despace also switches on the field-dependent focus term; see ``tinytim_skip_zero_terms``.
    tinytim_skip_zero_terms : bool
        Default True. Reproduces a TinyTim quirk: with zero despace, the field-dependent focus (Z4) of the
        instrument is dropped, because TinyTim only applies a field-dependent term when its static table value is
        nonzero, and the Z4 table value is rewritten as exactly 0. Set False to always apply the field focus.
    include_si_wfe : bool
        Default True. If False, the instrument's field-dependent aberration plane (``HSTAnnularZernikeWFE``) is
        omitted, as for the JWST instruments. The OTA mirror maps are still included.
    """

    telescope = 'HST'
    _default_detector = None  # set in subclasses; used to choose the initial pixel scale

    def __init__(self, name, pixelscale=None):
        self.siaf = stpsf_core.get_siaf_with_caching('HST')
        if pixelscale is None:
            pixelscale = self._default_pixelscale(self._default_detector)

        super().__init__(name, pixelscale=pixelscale)
        self.pupil_radius = constants.HST_PUPIL_DIAMETER / 2 * u.meter
        # Analytic OTA pupil; OPD maps and the IR cold mask are separate planes added in get_optical_system
        self.pupil = HSTPupil(radius=self.pupil_radius.to_value(u.m), rotation=constants.HST_TINYTIM_WFC3_ROTATION)
        self._aperturename = None
        self.despace_microns = 0.0  # secondary mirror despace; adds focus to the aberrations
        self.include_si_wfe = True  # include the instrument's field-dependent aberrations, as for JWST
        self.tinytim_skip_zero_terms = True  # reproduce TinyTim, which drops field terms whose static value is 0

    def _default_pixelscale(self, detector):
        """Default pixel scale in arcsec/pixel for the given detector. Implemented in subclasses."""
        raise NotImplementedError('needs to be subclassed.')

    def _get_pixelscale_from_apername(self, apername):
        """Mean pixel scale of a SIAF aperture. HST apertures do not define XSciScale/YSciScale, so it is
        measured from the sci-to-ideal transform."""
        return _siaf_mean_pixelscale(self.siaf[apername])

    @stpsf_core.SpaceTelescopeInstrument.aperturename.setter
    def aperturename(self, value):
        """Set the SIAF aperture name, with validation."""
        if value not in self.siaf.apertures:
            raise ValueError(f'Aperture name {value} not a valid SIAF aperture name for {self.name}')
        self._aperturename = value

    def _default_detector_position(self):
        """Center of the detector as defined by the SIAF reference pixel, as 0-based (X, Y)."""
        ap = self.siaf[self._aperturename]
        return (int(round(ap.XSciRef)) - 1, int(round(ap.YSciRef)) - 1)

    def _get_fits_header(self, result, options):
        """Populate FITS Header keywords"""
        super()._get_fits_header(result, options)
        result[0].header['DET_X'] = (self.detector_position[0], 'Detector X pixel position (0-based)')
        result[0].header['DET_Y'] = (self.detector_position[1], 'Detector Y pixel position (0-based)')
        result[0].header['DETECTOR'] = (self.detector, 'Detector selected')
        result[0].header['APERNAME'] = (self._aperturename, 'SIAF aperture name')

    def get_optical_system(self, *args, **kwargs):
        optsys = super().get_optical_system(*args, **kwargs)
        # The OTA mirror OPDs follow the pupil. The SI aberrations (added by the base class) come after them.
        primary, secondary = self._get_mirror_map_planes()
        optsys.planes[1:1] = [primary, secondary]
        self._extra_keywords['TEL_WFE'] = (self._rms_opd_nm(self.pupil, primary, secondary),
                                           '[nm] Telescope pupil RMS wavefront error')
        if self.channel == 'IR':
            cold_mask = WFC3IRColdMask(
                radius=self.pupil_radius.to_value(u.m), rotation=constants.HST_TINYTIM_WFC3_ROTATION
            )
            optsys.planes.insert(len(optsys.planes) - 1, cold_mask)  # just before the detector
        return optsys

    def _rms_opd_nm(self, pupil, *opd_planes):
        """RMS OPD in nm of the given planes over the open aperture of the pupil, on a 256-pixel grid"""
        wave = poppy.Wavefront(diam=self.pupil_radius.to_value(u.m) * 2.1, npix=256)
        mask = pupil.get_transmission(wave) > 0.5
        opd = sum(plane.get_opd(wave) for plane in opd_planes)
        return float(np.sqrt(np.mean(opd[mask] ** 2)) * 1e9)

    def _get_pupil_v2v3(self):
        """(V2, V3) in arcsec of the current detector position in TinyTim's "pupil" convention, i.e. minus the
        usual field V2, V3 of the SIAF."""
        ap = self.siaf[self._aperturename]
        v2, v3 = ap.sci_to_tel(self.detector_position[0] + 1, self.detector_position[1] + 1)
        return -v2, -v3

    def _get_mirror_map_planes(self):
        """OPD planes for the OTA primary mirror and the secondary mirror, whose map is shifted with field position"""
        primary_map, secondary_map, radius_pix = _load_hst_mirror_maps(self._STPSF_basepath)
        v2, v3 = self._get_pupil_v2v3()
        v2_ref, v3_ref = constants.HST_SECONDARY_SHIFT_REF_V2V3
        beam_radius_mm = constants.HST_SECONDARY_BEAM_DIAMETER_MM / 2
        shift = tuple(
            constants.HST_SECONDARY_SHIFT_MM_PER_TAN * np.tan(np.deg2rad((a - a_ref) / 3600)) / beam_radius_mm
            for a, a_ref in ((v2, v2_ref), (v3, v3_ref))
        )
        radius = self.pupil_radius.to_value(u.m)
        rotation = constants.HST_TINYTIM_WFC3_ROTATION
        return (
            HSTMirrorMap(primary_map, radius_pix, radius=radius, rotation=rotation, name='HST primary mirror map'),
            HSTMirrorMap(secondary_map, radius_pix, shift=shift, radius=radius, rotation=rotation,
                         name='HST secondary mirror map'),
        )

    def _get_field_zernikes(self):
        """Annular Zernike coefficients (waves at 547 nm, keyed by TinyTim term number) at the detector position.

        Like TinyTim, a field-dependent term is only applied if its static value is nonzero. TinyTim's tiny1 rewrites
        the Z4 table entry as ``%f`` after adding the despace, so the 1e-12 placeholder becomes exactly 0: the
        field-dependent focus is dropped unless the despace is nonzero. Set ``tinytim_skip_zero_terms = False`` to
        always apply it.
        """
        model = self._tinytim_model[self.detector]
        coeffs = {j + 1: float(value) for j, value in enumerate(model['static_zernikes']) if value != 0}
        coeffs[4] = round(coeffs.get(4, 0.0) + constants.HST_DESPACE_FOCUS_WAVES_PER_MICRON * self.despace_microns, 6)
        if coeffs[4] == 0 and not self.tinytim_skip_zero_terms:
            coeffs[4] = 1e-12
        x = self.detector_position[0] / model['scale_d']
        y = (self.detector_position[1] + model['y_offset']) / model['scale_d']
        for k, j in enumerate((4, 5, 6, 7, 8)):
            if coeffs.get(j, 0.0) == 0:
                continue
            field = np.polynomial.polynomial.polyval2d(x, y, model['field_coeffs'][k]) / model['convert']
            coeffs[j] = coeffs.get(j, 0.0) + field
        return coeffs

    def _get_aberrations(self):
        """Field-dependent aberrations of the instrument plus any secondary mirror despace"""
        if not self.include_si_wfe:
            return None
        return HSTAnnularZernikeWFE(
            self._get_field_zernikes(), radius=self.pupil_radius.to_value(u.m),
            rotation=constants.HST_TINYTIM_WFC3_ROTATION, name=f'{self.name} aberrations',
        )

    def _addAdditionalOptics(self, optsys, oversample=2):
        _log.debug(f'   No additional optics added for {self.name}')
        return optsys, False, None


class WFC3(HSTInstrument):
    """
    A class modeling the optics of the Hubble Space Telescope's Wide Field Camera 3.

    WFC3 has two channels, IR and UVIS, with three detectors: ``IR``, ``UVIS1`` and ``UVIS2``.
    The channel is selected automatically, in the same way as for `NIRCam`:

    * Setting the detector always sets the channel.
    * Setting a filter switches to the channel containing that filter, unless ``auto_channel = False``.
      Going from IR to UVIS selects ``UVIS1``; switching between UVIS filters keeps the current UVIS chip.
    * The channel is read-only (set filter or detector instead).

    Changing the detector resets ``detector_position`` to the center of the new detector.
    A custom ``pixelscale`` is preserved when the detector changes.

    The IR channel adds the cold pupil mask (a slightly smaller outer radius, a larger secondary obscuration and
    thicker spiders than the OTA). The UVIS channel has no extra pupil mask. Field-dependent aberrations
    (focus, astigmatism and coma, Z4-Z8) come from 2-D polynomials in detector position, one set per detector;
    UVIS1 uses the same polynomial as UVIS2 with the y position offset by 2048 pixels, as in TinyTim.

    .. warning::
        The model currently reproduces TinyTim, including its quirks. Notably the field-dependent focus is not applied
        unless ``despace_microns`` is nonzero (see `HSTInstrument`). Requires the WFC3 and HST data files described in
        the ``stpsf.hst`` module docstring.

    Parameters
    ----------
    use_tinytim_pixelscale : bool
        If False (default) the pixel scale is the mean scale of the SIAF aperture (about 0.128 arcsec for IR and
        0.0395 arcsec for UVIS). If True use TinyTim's nominal values (0.13 and 0.04 arcsec), which is useful for
        direct comparisons to TinyTim.
    """

    # NPIXELS as (ny, nx), following the convention of SpaceTelescopeInstrument._detector_npixels.
    # The UVIS size follows the position limits of TinyTim (X 0-4095, Y 0-2047).
    _detector_sizes = {'IR': 1014, 'UVIS1': (2048, 4096), 'UVIS2': (2048, 4096)}
    _default_detector = 'IR'
    _DETECTOR_APERTURES = {'IR': 'IIRFIX', 'UVIS1': 'IUVIS1FIX', 'UVIS2': 'IUVIS2FIX'}
    _model_filename = 'WFC3_tinytim_model.fits'

    def __init__(self, use_tinytim_pixelscale=False):
        # Set these up before the superclass __init__, so the overridden filter setter works when called from there
        self.auto_channel = False
        self.use_tinytim_pixelscale = use_tinytim_pixelscale

        self._check_data_present()
        super().__init__('WFC3')

        self._detectors = dict(self._DETECTOR_APERTURES)
        self._detector = 'IR'
        self._aperturename = self._detectors['IR']
        self._detector_npixels = self._detector_sizes['IR']
        self.detector_position = self._default_detector_position()

        self.auto_channel = True

    # ---- data handling
    def _check_data_present(self):
        base = utils.get_stpsf_data_path()
        missing = [d for d in ('HST', 'WFC3') if not os.path.isdir(os.path.join(base, d))]
        if missing:
            raise IOError(
                f"HST/WFC3 data directories {missing} not found in {base}. These are generated from TinyTim files by "
                "dev_utils/make_hst_data.py and are not yet part of the released STPSF data package."
            )

    @property
    def _tinytim_model(self):
        """TinyTim model parameters for all three WFC3 detectors, keyed by detector name"""
        return _load_tinytim_model(os.path.join(self._datapath, self._model_filename))

    def _get_filters(self):
        filter_list, filter_info = super()._get_filters()
        table = ioascii.read(os.path.join(self._STPSF_basepath, self.name, 'filters.tsv'))
        self._filter_channels = {row['filter']: row['channel'] for row in table}
        return filter_list, filter_info

    @property
    def ir_filter_list(self):
        """Filters of the IR channel"""
        return [f for f in self.filter_list if self._filter_channels[f] == 'IR']

    @property
    def uvis_filter_list(self):
        """Filters of the UVIS channel"""
        return [f for f in self.filter_list if self._filter_channels[f] == 'UVIS']

    # ---- channel / detector / pixel scale
    @property
    def channel(self):
        """Current channel, 'IR' or 'UVIS'. Not directly settable; set the filter or detector."""
        return 'IR' if self._detector == 'IR' else 'UVIS'

    @channel.setter
    def channel(self, value):
        raise RuntimeError('WFC3 channel is not directly settable; set filter or detector instead.')

    def _default_pixelscale(self, detector):
        if self.use_tinytim_pixelscale:
            return 0.13 if detector == 'IR' else 0.04
        return self._get_pixelscale_from_apername(self._DETECTOR_APERTURES[detector])

    @stpsf_core.SpaceTelescopeInstrument.detector.setter
    def detector(self, value):
        """Set the detector ('IR', 'UVIS1' or 'UVIS2'). 'UVIS' is accepted as a synonym for 'UVIS1'.

        This switches channel if needed and resets ``detector_position`` to the center of the new detector.
        """
        value = value.upper()
        if value == 'UVIS':
            value = 'UVIS1'
        if value not in self.detector_list:
            raise ValueError('Invalid detector. Valid detector names are: {}'.format(', '.join(self.detector_list)))
        self._set_detector(value)

    def _set_detector(self, value):
        if value == self._detector:
            return
        old_default_scale = self._default_pixelscale(self._detector)
        self._detector = value
        self._detector_npixels = self._detector_sizes[value]
        self._update_aperturename()
        if self.pixelscale == old_default_scale:
            self.pixelscale = self._default_pixelscale(value)
            _log.info(f'WFC3 pixel scale switched to {self.pixelscale:.5f} arcsec/pixel for detector {value}.')
        self.detector_position = self._default_detector_position()

    def _switch_channel(self, channel):
        """Switch to the IR or UVIS channel, keeping the current UVIS chip when already in UVIS."""
        if self.channel == channel:
            return
        _log.debug(f'Automatically changing WFC3 channel to {channel}')
        if channel == 'IR':
            self._set_detector('IR')
        elif channel == 'UVIS':
            self._set_detector('UVIS1')
        else:
            raise ValueError(f'Invalid WFC3 channel name: {channel}')

    def _update_aperturename(self):
        self.aperturename = self._detectors[self._detector]

    @stpsf_core.SpaceTelescopeInstrument.filter.setter
    def filter(self, value):
        """Set the filter, switching channel automatically unless ``auto_channel`` is False"""
        super(WFC3, self.__class__).filter.__set__(self, value)
        if self.auto_channel:
            self._switch_channel(self._filter_channels[self._filter])
