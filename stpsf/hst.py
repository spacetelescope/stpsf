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
   3 mirror support pads). Edge pixels get exact fractional coverage (poppy's ``filled_circle_aa`` for circles, an
   analytic strip formula for the spiders).
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
* **Pupil edges.** TinyTim's 512x512 pupil and wave maps are centered on pixel 256, half a pixel from the
  poppy grid center. We center everything on the poppy grid, which has no effect on PSF intensities. TinyTim draws
  its pupil (256 pixels across) with 8x8 subsampling and then rotates it by 135 degrees with bilinear
  interpolation, which softens the edges. Its open area is about 0.09% larger than the exact geometry (39552 versus
  39517 pixels at R=128); we use the exact geometry.
* **Annular Zernikes use epsilon = 0.33 even for the IR cold mask** (whose obscuration is 0.354), as TinyTim does.
* **Pixel scale.** SIAF HST apertures do not define ``XSciScale``; the mean scale is measured from the SIAF
  transform. TinyTim's nominal scales (0.13 and 0.04 arcsec) are available with ``use_tinytim_pixelscale``.

Validation
----------
The wavefront (pupil mask, mirror maps, Zernikes, secondary shift) was checked against TinyTim's dumped
``tinytim_wave.fits`` at the IR center position (the OPD in nm is independent of wavelength): 34.80 nm RMS versus
TinyTim's 34.84 nm, with a smooth, low-order difference of 0.14 nm RMS (at most 0.3 nm). Within the pupil, a plain
inverse FFT of TinyTim's dumped pupil and wave map reproduces its dumped critical PSF, and swapping in STPSF's OPD
changes the peak by only 0.02%.

Monochromatic WFC3 IR PSFs (1.6 micron, at (507, 507), (50, 50) and (960, 960)) match TinyTim's undistorted ``tiny2``
output to ~0.006% of the peak RMS (and 0.06% maximum at the center position; the test allows 0.5%), with PSF
centroids agreeing to 0.001 pixel. This requires
``detector_oversample >= 8`` (TinyTim's output is the PSF integrated over its pixels), an odd ``fov_pixels`` so the
PSF is centered on a pixel as in TinyTim, and normalizing both images to their own field sum (TinyTim normalizes to
the field it writes, whereas STPSF keeps flux outside the field, about 15% at 3 arcsec). The regression test needs
the reference PSFs from ``make_reference_psfs.sh`` and skips unless ``$TINYTIM_REF`` points at them.

Polychromatic PSFs (set ``tinytim_wavelengths = True``) use TinyTim's wavelengths, which are the points tabulated in the
filter file (generated from ``wfc3ir.fil`` / ``wfc3uvis.fil``), and its weights, throughput * photlam * dlambda. They
reproduce the wavelengths and weights in TinyTim's ``.par`` files to their 6 printed digits for F105W, F125W and F160W
with a G8V spectrum. The resulting PSFs agree with TinyTim's ``tiny2`` output (using the exact pixel scale from the
``.par`` file; the FITS header's rounded scale gives worse agreement) to about 0.003% of the peak RMS, and at most
0.12% (F105W), 0.08% (F125W) and 0.04% (F160W) of the peak at the center position. The tests need ``$TINYTIM`` as well,
for the G8V spectrum. Quirks reproduced: TinyTim's dlambda uses the full
spacing, not half, at the two ends of the table; wavelengths below 1% of the peak weight are skipped; and filters
tabulated at one wavelength are monochromatic. Without ``tinytim_wavelengths``, STPSF's usual ``nlambda`` equal-width
bins and synphot weights are used.

Residual differences from TinyTim (IR, ``tiny2``). The remaining differences are small and depend on the filter,
not on detector position. STPSF puts 0.05-0.13% less of the flux in the core (r < 0.3 arcsec) and 0.03% more
in the wings (r > 1.5 arcsec); the peak pixel is lower by 0.04-0.12%, shorter wavelengths being worse. The PSF
centroids agree to 0.003 pixel, and a fitted scale (<= 0.12%), shift (<= 0.001 pixel) and blur account for little
of it. Tested and ruled out as the cause: flux normalization (also per wavelength), the pupil grid size (TinyTim's
256 pixel pupil gives a *larger* difference), the wavefront (see above), and the
details of TinyTim's resampling (an emulation of its FFT sinc-MTF resampling, which reproduces ``tiny2`` exactly
when fed TinyTim's own pupil, does not move STPSF toward TinyTim). The identified origin is the pupil mask: TinyTim's
open area is 0.09% larger than the exact geometry, from its 8x8 subsampled, bilinearly rotated pupil, and
feeding TinyTim's dumped pupil into the emulation removes the difference completely. Poppy's resampling of the
PSF onto pixels differs from TinyTim's and partly offsets this, so the current agreement should not be taken to
improve with more oversampling: ``detector_oversample=8`` is closest, partly by cancellation, and finer
oversampling moves the peak by about 0.1% in the other direction.

Detector effects: ``DET_DIST`` / ``OVERDIST`` extensions add SIAF geometric distortion (resampled onto the SIAF
detector grid; the SIAF pixel scale of 0.1282 arcsec is the mean of anisotropic 0.1355 x 0.1211 arcsec pixels, so
compare to TinyTim's ``tiny3`` output at that scale, not 0.13) and, for the IR channel, TinyTim's 3x3 IPC kernel
(centre 0.897, edges 0.025, corners 0.0007). Against ``tiny3`` for F160W at three positions, centroids agree to
0.006 pixel, with residuals of about 0.1% of the peak RMS and at most 2.7% at the peak. The largest differences
are near the peak and at the edges. TinyTim and SIAF use different distortion models, and the ``tiny3`` images have
dark bands along some edges that are not present in ``tiny2``, apparently from its distortion mapping. Disable with
``options['add_distortion']`` and ``options['add_ipc']``.

For the UVIS channel the detector step instead convolves with TinyTim's 3x3 charge diffusion kernel (disable with
``options['add_charge_diffusion']``), built as in ``tiny3``: the kernels tabulated at 0.4, 0.6 and 0.8 micron,
interpolated linearly in wavelength and weighted by the spectral weights (retained from the calculation), then
convolved with a field-dependent Gaussian kernel evaluated at the first wavelength (sigma 0.45816 pixel everywhere,
as the tabulated field polynomials are constant). ``tiny3`` quirks reproduced: wavelengths above 0.4 micron all
use the 0.6 and 0.8 micron kernel pair, so shorter ones are extrapolated. As for the IR IPC, the kernel is applied
to the oversampled image after upsampling it. UVIS PSFs match ``tiny2`` to under 0.15% of the peak (rms
< 2e-5 of the peak) for monochromatic and F275W/F438W/F606W/F814W G8V cases on both detectors, and ``tiny3`` to
2-3% of the peak (rms < 5e-4), the difference being in the undersampled core from the distortion resampling.
See ``stpsf/tests/test_hst.py`` for the tests that run in CI, ``dev_utils/run_hst_regression_tests.py`` for the
regression tests against TinyTim output (which need a TinyTim installation), and the TinyTim notes
(``notes/03_wfc3_ir_porting_guide.md``) for the full validation plan.
"""

import logging
import os
from functools import lru_cache

import astropy.io.fits as fits
import astropy.units as u
import numpy as np
import poppy
import synphot
from astropy.io import ascii as ioascii
from scipy import ndimage, signal

from . import constants, detectors, distortion, stpsf_core, utils
from .distortion import _siaf_mean_pixelscale

_log = logging.getLogger('stpsf')


def _sum_of_uniforms_cdf(z, a, b):
    """CDF of ``a*U1 + b*U2``, with U1 and U2 independent and uniform on [-1/2, 1/2] and ``a >= b >= 0``."""
    if b < 1e-9:
        return np.clip(z / a + 0.5, 0, 1)

    def ramp2(t):
        return np.maximum(t, 0) ** 2

    return (ramp2(z + (a + b) / 2) - ramp2(z + (a - b) / 2) - ramp2(z - (a - b) / 2) + ramp2(z - (a + b) / 2)) / (
        2 * a * b
    )


def _strip_coverage(x, y, normal, half_width):
    """Exact fraction of each unit-square pixel covered by an infinite strip through the origin.

    Parameters
    ----------
    x, y : ndarray
        Pixel center coordinates, in pixels.
    normal : (float, float)
        Unit vector (x, y) perpendicular to the strip.
    half_width : float
        Half the strip width, in pixels.
    """
    # The strip covers the pixel where |d0 + s| < half_width, with s the projection onto the normal of a point
    # uniformly distributed within the pixel, relative to its center.
    nx, ny = normal
    a, b = max(abs(nx), abs(ny)), min(abs(nx), abs(ny))
    d0 = nx * x + ny * y
    return _sum_of_uniforms_cdf(half_width - d0, a, b) - _sum_of_uniforms_cdf(-half_width - d0, a, b)


class HSTPupil(poppy.AnalyticOpticalElement):
    """Analytic HST OTA entrance pupil in TinyTim's geometry.

    The geometry is defined in the OTA V2/V3 frame, with +V3 along the array columns (poppy x) and
    +V2 along the array rows (poppy y), and is then rotated clockwise by ``rotation`` degrees to
    the camera frame, following TinyTim's convention. Subclasses (e.g. `WFC3IRColdMask`) override
    `_geometry` to model other masks in the same frame.

    Edge pixels are given their exact fractional coverage: circles use `poppy.geometry.filled_circle_aa` and the
    spiders, which are strips of any orientation, use an analytic formula. The fractions of the different features
    are multiplied, which is only approximate for a pixel touched by more than one feature (e.g. where a pad
    meets a spider).

    Parameters
    ----------
    radius : float
        Outer radius of the OTA pupil in meters.
    rotation : float
        Clockwise rotation in degrees.
    """

    default_name = 'HST OTA pupil'

    def __init__(self, radius=constants.HST_PUPIL_DIAMETER / 2, rotation=constants.HST_TINYTIM_WFC3_ROTATION,
                 name=None, **kwargs):
        super().__init__(name=name or self.default_name, planetype=poppy.poppy_core.PlaneType.pupil, **kwargs)
        self.pupil_diam = 2 * radius
        self.radius = radius
        self.rotation_cw = rotation

    def _geometry(self):
        """Dimensions in units of the pupil radius: outer and secondary radii, spider width, and the pads as
        (V3, V2, radius) tuples, all in the OTA V2/V3 frame"""
        return {
            'outer_radius': 1.0,
            'secondary_radius': constants.HST_OTA_SECONDARY_RADIUS,
            'spider_width': constants.HST_OTA_SPIDER_WIDTH,
            'pads': constants.HST_OTA_PADS,
        }

    def get_transmission(self, wave):
        geom = self._geometry()
        y, x = self.get_coordinates(wave)
        pix = wave.pixelscale.to_value(u.m / u.pixel) if hasattr(wave.pixelscale, 'unit') else wave.pixelscale
        # Work in pixel units so that filled_circle_aa, which assumes unit pixels, is exact
        x, y = x / pix, y / pix
        radius_pix = self.radius / pix
        sin_t, cos_t = np.sin(np.deg2rad(self.rotation_cw)), np.cos(np.deg2rad(self.rotation_cw))

        def disk(xc, yc, radius):
            return poppy.geometry.filled_circle_aa(x.shape, xc, yc, radius * radius_pix, xarray=x, yarray=y)

        trans = disk(0, 0, geom['outer_radius']) * (1 - disk(0, 0, geom['secondary_radius']))
        for pad_v3, pad_v2, pad_r in geom['pads']:
            # the frame conversion is its own inverse, so it also maps pad centers from OTA to camera coordinates
            xc, yc = pad_v2 * sin_t + pad_v3 * cos_t, pad_v2 * cos_t - pad_v3 * sin_t
            trans *= 1 - disk(xc * radius_pix, yc * radius_pix, pad_r)
        half_width = geom['spider_width'] / 2 * radius_pix
        for normal in ((sin_t, cos_t), (cos_t, -sin_t)):  # normals to the spiders along V3 and along V2
            trans *= 1 - _strip_coverage(x, y, normal, half_width)
        return trans


class WFC3IRColdMask(HSTPupil):
    """The WFC3 IR channel's cold pupil mask, a physical element separate from the OTA pupil.

    It sits in the WFC3 IR optical path, undersizing the OTA pupil slightly (outer edge, a larger secondary
    obscuration and wider spiders), and is applied in addition to `HSTPupil`. It is registered with the OTA pupil
    through the shared V2/V3 frame and rotation. Dimensions are in ``constants.HST_WFC3_IR_COLD_MASK``; there are no
    support pads in the mask.
    """

    default_name = 'WFC3 IR cold mask'

    def _geometry(self):
        cm = constants.HST_WFC3_IR_COLD_MASK
        return {
            'outer_radius': cm['outer_radius'],
            'secondary_radius': cm['secondary_radius'],
            'spider_width': cm['spider_width'],
            'pads': (),
        }


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


def _tinytim_weights(wave_angstrom, throughput, photlam, weight_limit=0.01):
    """Wavelengths and weights for a polychromatic PSF, computed as TinyTim does (``spectrum.c``, ``polypsf.c``).

    Parameters
    ----------
    wave_angstrom, throughput : array_like
        The tabulated filter curve; its wavelengths are the ones used for the PSFs.
    photlam : array_like
        The source spectrum in photons per unit wavelength, at ``wave_angstrom``.
    weight_limit : float
        TinyTim's "smart skip": wavelengths with a weight below this fraction of the largest are dropped.

    Returns
    -------
    wavelengths_m, weights : ndarray
        Wavelengths in meters, and weights normalized to sum to 1.

    Notes
    -----
    The weight is throughput * photlam * dlambda, normalized. TinyTim's dlambda is half the span between the two
    neighboring wavelengths for interior points, but the *full* spacing to the single neighbor at either end of the
    table, which we reproduce. The weights are renormalized after skipping low-weight wavelengths.
    """
    wave = np.asarray(wave_angstrom, dtype=float)
    if wave.size == 1:
        return wave * 1e-10, np.ones(1)
    dlambda = np.empty_like(wave)
    dlambda[1:-1] = np.abs(wave[:-2] - wave[2:]) / 2
    dlambda[0], dlambda[-1] = abs(wave[1] - wave[0]), abs(wave[-1] - wave[-2])
    weights = np.asarray(throughput, dtype=float) * np.asarray(photlam, dtype=float) * dlambda
    weights /= weights.sum()
    keep = weights >= weight_limit * weights.max()
    return wave[keep] * 1e-10, weights[keep] / weights[keep].sum()


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
    tinytim_wavelengths : bool
        Default False, which uses STPSF's usual ``nlambda`` equal-width bins and synphot weights. If True, a
        polychromatic calculation uses TinyTim's wavelengths (the points tabulated in the filter file) and TinyTim's
        weights, for reproducing TinyTim; see ``_tinytim_weights``.
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
        self.tinytim_wavelengths = False  # use TinyTim's wavelengths and weights rather than nlambda bins

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

    def _get_weights(self, source=None, nlambda=None, monochromatic=None, verbose=False):
        """Wavelengths and weights for a polychromatic calculation.

        By default this is the generic STPSF behavior: ``nlambda`` bins across the filter, weighted with synphot.
        With ``tinytim_wavelengths`` set, the wavelengths are instead TinyTim's, i.e. the points tabulated in the
        filter throughput file, weighted as TinyTim does (see ``_tinytim_weights``). Filters tabulated at a single
        wavelength are always treated as monochromatic, as in TinyTim.
        """
        if monochromatic is not None:
            return super()._get_weights(source=source, nlambda=nlambda, monochromatic=monochromatic, verbose=verbose)
        with fits.open(self._filters[self.filter].filename) as hdul:
            wave_angstrom, throughput = hdul[1].data['WAVELENGTH'], hdul[1].data['THROUGHPUT']
        if len(wave_angstrom) == 1:
            return _tinytim_weights(wave_angstrom, throughput, [1.0])
        if not self.tinytim_wavelengths or not (source is None or isinstance(source, synphot.SourceSpectrum)):
            return super()._get_weights(source=source, nlambda=nlambda, monochromatic=monochromatic, verbose=verbose)
        if source is None:
            source = synphot.SourceSpectrum(synphot.models.BlackBodyNorm1D, temperature=5700 * u.K)
            _log.info('No source spectrum supplied, therefore defaulting to 5700 K blackbody')
        photlam = source(wave_angstrom * u.AA, flux_unit='photlam').value
        return _tinytim_weights(wave_angstrom, throughput, photlam)

    def _calc_psf_format_output(self, result, options):
        """Format the output HDUList, adding detector effects.

        As for the JWST instruments, if ``options['add_distortion']`` is True (the default) copies of the
        extensions are made with names ending in ``DIST`` (``OVERDIST``, ``DET_DIST``). SIAF geometric distortion
        is applied to these, and then for the IR channel the IPC kernel is applied if ``options['add_ipc']`` is
        True (the default), after binning to detector pixels.
        """
        add_distortion = options.get('add_distortion', True)
        add_ipc = options.get('add_ipc', True)

        if add_distortion:
            _log.info('Adding SIAF distortion')
            n_exts = len(result)
            for ext in range(n_exts):
                hdu_new = fits.ImageHDU(result[ext].data.copy(), result[ext].header.copy())
                hdu_new.header['EXTNAME'] = result[ext].header['EXTNAME'][0:4] + 'DIST'
                hdu_new.data = self._apply_siaf_distortion(hdu_new)
                hdu_new.header.add_history('Applied SIAF geometric distortion')
                result.append(hdu_new)

        super()._calc_psf_format_output(result, options)

        if add_ipc and add_distortion and self.channel == 'IR':
            for extname in ('DET_DIST', 'OVERDIST'):
                if extname in result:
                    self._apply_ir_ipc(result[extname])
        if options.get('add_charge_diffusion', True) and add_distortion and self.channel == 'UVIS':
            kernel = self._uvis_diffusion_kernel(result[0].header)
            for extname in ('DET_DIST', 'OVERDIST'):
                if extname in result:
                    self._apply_kernel(result[extname], kernel, 'WFC3-UVIS', 'charge diffusion')

    def _uvis_diffusion_kernel(self, header):
        """Combined 3x3 UVIS charge diffusion kernel, as TinyTim's ``tiny3`` builds it (``Compute_wfc3_kernel``,
        ``Compute_kernel_xy``): the tabulated kernels interpolated linearly in wavelength and weighted by the
        spectral weights (the ``WAVE#``, ``WGHT#`` keywords of the PSF header), convolved with the field-dependent
        Gaussian kernel evaluated at the *first* wavelength.
        """
        model = self._tinytim_model[self.detector]
        nwaves = header['NWAVES']
        waves_m = np.array([header[f'WAVE{i}'] for i in range(nwaves)])
        weights = np.array([header[f'WGHT{i}'] for i in range(nwaves)])
        waves_um = waves_m * 1e6

        kernel_waves, kernels = model['kernel_waves_um'], model['kernels']
        weighted = np.zeros((3, 3))
        for wave, weight in zip(waves_um, weights):
            # tiny3 quirk, reproduced: any wavelength above the first tabulated one uses the second and third kernels
            iw = 0
            while wave > kernel_waves[iw] and iw < len(kernel_waves) - 2:
                iw += 1
            frac = (wave - kernel_waves[iw]) / (kernel_waves[iw + 1] - kernel_waves[iw])
            kernel = kernels[iw] + frac * (kernels[iw + 1] - kernels[iw])
            weighted += weight * kernel / kernel.sum()

        # field-dependent Gaussian; polynomial coefficients are tabulated at two wavelengths
        x = self.detector_position[0] / model['scale_d'] + 4
        y = (self.detector_position[1] + model['y_offset']) / model['scale_d'] + 4
        sigmas = [np.polynomial.polynomial.polyval2d(x, y, c) for c in model['kfield_coeffs']]
        wave0, kfield_waves = waves_um[0], model['kfield_waves_um']
        if wave0 < kfield_waves[0]:
            sigma = sigmas[0] + (sigmas[1] - sigmas[0]) / (kfield_waves[1] - kfield_waves[0]) * (wave0 - kfield_waves[0])
        else:
            sigma = sigmas[-1]
        if sigma > 0:
            offsets = np.arange(-1, 2)
            r2 = offsets[:, None] ** 2 + offsets[None, :] ** 2
            gauss = np.exp(-0.5 * r2 / sigma**2)
            gauss /= gauss.sum()
        else:
            gauss = np.zeros((3, 3))
            gauss[1, 1] = 1
        return signal.convolve2d(weighted, gauss)

    def _apply_kernel(self, hdu, kernel, instrument, description):
        """Convolve an extension, in place, with a detector kernel defined at detector pixel sampling"""
        oversample = hdu.header['OVERSAMP']
        if oversample != 1:
            kernel = detectors.oversample_ipc_model(kernel, oversample)
        hdu.data = signal.fftconvolve(hdu.data, kernel, mode='same')
        hdu.header['DETKERN'] = (instrument, f'Detector {description} kernel applied')
        hdu.header.add_history(f'Applied {description} kernel, from TinyTim')

    def _apply_siaf_distortion(self, hdu):
        """Resample an image from the ideal frame onto the detector pixel grid of this SIAF aperture.

        The PSF header's DET_X, DET_Y are 0-based, while pysiaf science coordinates are 1-based.
        """
        hdu_tmp = fits.ImageHDU(hdu.data, hdu.header.copy())
        hdu_tmp.header['DET_X'] += 1
        hdu_tmp.header['DET_Y'] += 1
        return distortion.distort_image(fits.HDUList([fits.PrimaryHDU(), hdu_tmp]), ext=1,
                                        aper=self.siaf[self._aperturename])

    def _apply_ir_ipc(self, hdu):
        """Convolve an extension with TinyTim's WFC3/IR interpixel capacitance kernel, in place"""
        kernel = self._tinytim_model['IR']['kernels'][0]  # identical at all tabulated wavelengths
        oversample = hdu.header['OVERSAMP']
        if oversample != 1:
            kernel = detectors.oversample_ipc_model(kernel, oversample)
        hdu.data = signal.fftconvolve(hdu.data, kernel, mode='same')
        hdu.header['IPCINST'] = ('WFC3-IR', 'Instrument whose IPC kernel was applied')
        hdu.header.add_history('Applied detector interpixel capacitance (IPC) model, from TinyTim')

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
