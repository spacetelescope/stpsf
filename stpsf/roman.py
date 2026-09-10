"""
=================
Roman Instruments
=================

WARNING: This model has not yet been validated against other PSF
         simulations and uses several approximations.
"""

import logging
import os.path
import pprint
import re
from pathlib import Path

import astropy.units as u
import numpy as np
import poppy
from astropy.io import fits
from scipy.interpolate import griddata

from . import distortion, utils, stpsf_core, detectors, constants

# Monkey patch: double the cache size over the default maxsize=128
# to prevent slowdowns from increased number of Zernike coefficients
from functools import lru_cache
poppy.zernike.cached_zernike1 = lru_cache(maxsize=256)(poppy.zernike.cached_zernike1.__wrapped__)

_log = logging.getLogger('stpsf')

GRISM_FILTERS = ('GRISM0', 'GRISM1')
PRISM_FILTERS = ('PRISM',)


def _wfi_sci_xy_to_fp(x, y):
    """
    Convert from (x, y) in pixel coordinates to the field point numbering used
    for WFI Zernikes and pupil masks.

    Inverse of _wfi_fp_to_sci_xy().

    Parameters
    ----------
    x, y : float
        Pixel coordinates in Science frame. Values expected to be within 0-4096.

    Returns
    -------
    fp : int
        Field point index from 1 to 25.
    """
    n = WFI.NPIXELS
    vertical = np.round(4 * y / n + 1)
    horizontal = np.round(4 * (n - x) / n)
    fp = int(horizontal * 5 + vertical)
    return fp


def _wfi_fp_to_sci_xy(fp):
    """
    Convert from field point number to Science frame X, Y pixel coordinates.

    Inverse of _wfi_sci_xy_to_fp().

    Parameters
    ----------
    fp : int
        Field point index from 1 to 25.

    Returns
    -------
    (x, y) : tuple
        Pixel coordinates in Science frame. Values expected to be within 0-4096.
    """
    n = WFI.NPIXELS
    y = int(np.mod(fp - 1, 5) * n // 4)
    x = int(4 - (fp - 1) // 5) * n // 4
    return (x, y)


class WavelengthDependenceInterpolator(object):
    """
    WavelengthDependenceInterpolator can be configured with
    `n_zernikes` worth of Zernike coefficients at up to `n_wavelengths`
    wavelengths, and will let you `get_aberration_terms` for any
    wavelength in range interpolated linearly between measured/known
    points.

    For an imaging filter in STPSF's default reference WFI data,
    n_wavelengths=3 and n_zernikes=210. For a grism or prism filter,
    n_wavelengths=16 and n_zernikes=45. See the docstring for
    `build_detector_from_table()` in this file for more information.
    """

    def __init__(self, n_wavelengths, n_zernikes):
        self._n_wavelengths = n_wavelengths
        self._n_zernikes = n_zernikes
        self._aberration_terms = np.zeros((n_wavelengths, n_zernikes), dtype=np.float64)
        self._wavelengths = []

    def set_aberration_terms(self, wavelength, zernike_array):
        """Supply a reference `wavelength` and a `zernike_array`
        (of length `n_zernikes`) where the aberration is known
        """
        n_wavelengths_set = len(self._wavelengths)
        if wavelength not in self._wavelengths and n_wavelengths_set < self._n_wavelengths:
            self._wavelengths.append(wavelength)
            aberration_row_idx = n_wavelengths_set  # which is now index of last row
        elif wavelength in self._wavelengths:
            aberration_row_idx = self._wavelengths.index(wavelength)
        else:
            # can't add more wavelengths without allocating new _aberration_terms array
            raise ValueError(
                f"Already have information at {self._n_wavelengths} "
                "wavelengths (pass larger n_wavelengths to __init__?)")
        if len(zernike_array) != self._n_zernikes:
            raise ValueError(
                f"Expected {self._n_zernikes} aberration terms (pass different "
                "n_zernikes to __init__?)")
        self._aberration_terms[aberration_row_idx] = zernike_array

    def get_aberration_terms(self, wavelength):
        """Return the Zernike coefficients as interpolated for this
        `wavelength`"""
        # return array of length n_zernikes interpolated for this wavelength
        if wavelength in self._wavelengths:
            # aberration known exactly for this wavelength
            aberration_row_idx = self._wavelengths.index(wavelength)
            return self._aberration_terms[aberration_row_idx]
        else:
            # we have to interpolate @ this wavelength
            aberration_terms = griddata(self._wavelengths, self._aberration_terms, wavelength, method='linear')
            if np.any(np.isnan(aberration_terms)):
                if isinstance(wavelength, u.Quantity):
                    wavelength = wavelength.to(u.m).value
                wavelength_closest = np.clip(wavelength, np.min(self._wavelengths), np.max(self._wavelengths))
                _log.warn(
                    'Attempted to get aberrations at wavelength {:.2g} '
                    'outside the range of the reference data; clipping to closest wavelength {:.2g}'.format(
                        wavelength, wavelength_closest
                    )
                )

                aberration_terms = griddata(self._wavelengths, self._aberration_terms, wavelength_closest, method='linear')
            return aberration_terms


class FieldDependentAberration(poppy.ZernikeWFE):
    """
    FieldDependentAberration incorporates aberrations that
    are interpolated in wavelength, x, and y pixel positions by
    computing the Zernike coefficients for a particular wavelength
    and position.

    By default, `get_aberration_terms` will zero out Z1, Z2, and Z3
    (piston, tip, and tilt) as they are not meaningful for telescope
    PSF calculations (the former is irrelevant; the latter two would
    be handled by a distortion solution). Change
    `_omit_piston_tip_tilt` to False to include the Z1-3 terms.
    """
    _omit_piston_tip_tilt = True
    _field_position = None

    def __init__(self, pixel_width, pixel_height, name='Field-dependent Aberration',
                 radius=1.0, oversample=1, interp_order=3):
        self.pixel_width, self.pixel_height = pixel_width, pixel_height
        self.field_position = pixel_width // 2, pixel_height // 2
        self._wavelength_interpolators = {}
        self.pupil_diam = radius * 2.0
        super().__init__(name=name, verbose=True, radius=radius, oversample=oversample, interp_order=interp_order)

    def get_opd(self, wave):
        """Set the Zernike coefficients (for ZernikeWFE.getOPD) based
        on the wavelength of the incoming wavefront and the pixel
        position
        """
        if not isinstance(wave, poppy.Wavefront):
            wavelength = wave
        else:
            wavelength = wave.wavelength
        self.coefficients = self.get_aberration_terms(wavelength) * u.meter
        return super().get_opd(wave)

    @property
    def field_position(self):
        return self._field_position

    @field_position.setter
    def field_position(self, position):
        """Set the x and y pixel position on the detector for which to
        interpolate aberrations"""
        x_pixel, y_pixel = position
        if x_pixel > self.pixel_width or x_pixel < 0:
            raise ValueError('Requested pixel_x position lies outside ' 'the detector width ({})'.format(x_pixel))
        if y_pixel > self.pixel_height or y_pixel < 0:
            raise ValueError('Requested pixel_y position lies outside ' 'the detector height ({})'.format(y_pixel))

        self._field_position = x_pixel, y_pixel

    def add_field_point(self, x_pixel, y_pixel, interpolator):
        """Supply a wavelength-space interpolator for a pixel position
        on the detector"""
        self._wavelength_interpolators[(x_pixel, y_pixel)] = interpolator

    def get_aberration_terms(self, wavelength):
        """Supply the Zernike coefficients for the aberration based on
        the wavelength and pixel position on the detector"""
        if self.field_position in self._wavelength_interpolators:
            # short path: this is a known point
            interpolator = self._wavelength_interpolators[self.field_position]
            coefficients = interpolator.get_aberration_terms(wavelength)
        else:
            # get aberrations at all field points
            field_points, aberration_terms = [], []
            for field_point_coords, point_interpolator in self._wavelength_interpolators.items():
                field_points.append(field_point_coords)
                aberration_terms.append(point_interpolator.get_aberration_terms(wavelength))
            aberration_array = np.asarray(aberration_terms)
            assert len(aberration_array.shape) == 2, (
                'computed aberration array is not 2D ' '(inconsistent number of Zernike terms ' 'at each point?)'
            )
            field_position = tuple(self.field_position)
            coefficients = griddata(np.asarray(field_points), np.asarray(aberration_terms), field_position, method='linear')
            if np.any(np.isnan(coefficients)):
                # Find two closest input grid points
                # [decomposed hypot args are (xdiffs, ydiffs)]
                dist = np.hypot(*(np.asarray(field_points)
                                  - np.asarray(field_position)).T)
                min_dist_indx = np.argsort(dist)[:2]  # keep two closest points
                # Define line between two points, then find orthogonal line at
                # point of interest and intersection of these two lines
                x1, y1 = field_points[min_dist_indx[0]]
                x2, y2 = field_points[min_dist_indx[1]]
                dx = x2 - x1
                dy = y2 - y1
                a = (dy * (field_position[1] - y1) + dx * (field_position[0] - x1)) / (dx * dx + dy * dy)
                closest_interp_point = (x1 + a * dx, y1 + a * dy)
                # Interpolate aberrations to closest interpolated point
                coefficients = griddata(
                    np.asarray(field_points), np.asarray(aberration_terms), closest_interp_point, method='linear'
                )
                # If closest interpolatd point is still outside the input grid,
                # use nearest grid point instead
                if np.any(np.isnan(coefficients)):
                    coefficients = aberration_terms[min_dist_indx[0]]
                    _log.warn(
                        'Attempted to get aberrations at field point {} which is outside the range '
                        'of the reference data; approximating to nearest input grid point'.format(field_position)
                    )
                else:
                    _log.warn(
                        'Attempted to get aberrations at field point {} which is outside the range '
                        'of the reference data; approximating to nearest interpolated point {}'.format(
                            field_position, closest_interp_point
                        )
                    )
                assert not np.any(np.isnan(coefficients)), 'Could not compute aberration ' 'at field point {}'.format(
                    field_position
                )
        if self._omit_piston_tip_tilt:
            _log.debug('Omitting piston/tip/tilt')
            coefficients[:3] = 0.0  # omit piston, tip, and tilt Zernikes
        return coefficients


def _load_wfi_detector_aberrations(filename):
    from astropy.io import ascii

    zernike_table = ascii.read(filename, encoding='utf-8-sig')
    detectors_dict = {}

    det_col = 'sca'
    fp_col = 'fov'
    wave_col = 'wave'
    det_pos_x_col = 'SCI_X'  # already in pixels
    det_pos_y_col = 'SCI_Y'  # already in pixels
    calc_pix = False
    conv_wv_to_m = 1e-9  # nm to m
    # (Zernike units set separately in FieldDependentAberration.get_opd())

    def build_detector_from_table(number, zernike_table):
        """
        Build a FieldDependentAberration optic for a detector using number of
        Zernike terms found per row for specified wavelength and field points.

        In STPSF's default reference WFI data, imaging filters contain
        aberration data at the low, middle, and high wavelengths of their
        bandpasses. There are 16 unique wavelengths with recorded aberration
        data across the WFI's imaging filters. The prism, grism0, and grism1
        filters each contain aberration data for all 16 such wavelengths.

        Each imaging filter contains 25 field points per detector and Zernike
        coefficients up to Z210 for each detector/field point/intra-filter
        wavelength combination. The grism/prism filters contain 5 field points
        per detector and Zernike coefficients up to Z45 for each
        detector/field point/intra-filter wavelength combination.
        """
        n_zernikes = len([c for c in zernike_table.columns
                          if re.match(r'Z\d+', c)])
        single_detector_info = zernike_table[zernike_table[det_col] == number]
        field_points = set(single_detector_info[fp_col])
        detector = FieldDependentAberration(
            WFI.NPIXELS, WFI.NPIXELS, radius=constants.ROMAN_PUPIL_DIAMETER/2,
            name=f"Field Dependent Aberration (WFI{number:02d})"
        )
        for field_id in field_points:
            field_point_rows = single_detector_info[single_detector_info[fp_col] == field_id]
            local_x, local_y = (field_point_rows[0][det_pos_x_col],
                                field_point_rows[0][det_pos_y_col])
            interpolator = build_wavelength_dependence(field_point_rows,
                                                       n_zernikes)

            detector.add_field_point(local_x, local_y, interpolator)
        return detector

    def build_wavelength_dependence(rows, n_zernikes):
        """Build an interpolator object that interpolates `n_zernikes` Zernike
        terms in wavelength space"""
        wavelengths = set(rows[wave_col])
        interpolator = WavelengthDependenceInterpolator(
            n_wavelengths=len(wavelengths),
            n_zernikes=n_zernikes
        )
        for row in rows:
            z = np.array([row['Z{}'.format(idx + 1)]
                          for idx in range(n_zernikes)])
            interpolator.set_aberration_terms(row[wave_col] * conv_wv_to_m, z)

        return interpolator

    detector_ids = set(zernike_table[det_col])
    for detid in detector_ids:
        detid = int(detid) # temp to support g/prism Zernike CSVs
        detectors_dict[f"WFI{detid:02}"] = build_detector_from_table(detid, zernike_table)

    return detectors_dict


class _RomanInstrumentOptionsDict(dict):
    """
    A wrapper for the `RomanInstrument.options` dict that prevents the
    `add_distortion` key from being modified since distortion can't be added to
    PSFs from Roman's WFI (whose input data is already distorted) or Coronagraph
    Instrument (where distortion is not implemented).

    The default value set for `add_distortion` in `RomanInstrument` should be
    something other than True or False that indicates to users that distortion
    can't be toggled in the Roman instrument model.
    """
    def __setitem__(self, key, value):
        if key == 'add_distortion' and value != 'NA':
            _log.warn('Roman: add_distortion disabled for Roman PSFs')
        else:
            super().__setitem__(key, value)


@utils.combine_docstrings
class RomanInstrument(stpsf_core.SpaceTelescopeInstrument):
    """
    RomanInstrument contains data and functionality common to Roman
    instruments, such as setting the pupil shape
    """
    telescope = 'Roman'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.pupil_radius = constants.ROMAN_PUPIL_DIAMETER / 2 * u.meter

        self.siaf = stpsf_core.get_siaf_with_caching('roman')
        self._aperturename = None

        # reassign self.options from dict to _RomanInstrumentOptionsDict,
        # mainly to prevent user modification of the add_distortion key
        self.options['add_distortion'] = 'NA'  # distortion can't be toggled
        self.options = _RomanInstrumentOptionsDict(self.options.copy())

        self.options['jitter'] = 'gaussian'
        self.options['jitter_sigma'] = constants.ROMAN_TYPICAL_LOS_JITTER_PER_AXIS
        # arcsec/axis, see https://github.com/RomanSpaceTelescope/roman-technical-information/tree/main/data/Observatory/MissionandObservatoryTechnicalOverview#telescope-parameters

    def calc_psf(
        self,
        outfile=None,
        source=None,
        nlambda=None,
        monochromatic=None,
        fov_arcsec=None,
        fov_pixels=None,
        oversample=None,
        detector_oversample=None,
        fft_oversample=None,
        overwrite=True,
        display=False,
        save_intermediates=False,
        return_intermediates=False,
        normalize='first',
        add_distortion=None,
        crop_psf=False,
    ):
        """
        Compute a PSF.

        Note that Roman WFI PSFs, unlike those from the JWST instruments,
        always include distortion. Additionally, distortion is not implemented
        in STPSF for the Roman Coronagraph Instrument. As such, the
        `add_distortion` argument of the Roman implementation of `calc_psf()`
        will not affect the simulated PSF.

        Parameters
        ----------
        add_distortion : None
            Included for backward compatibility, but this argument has no
            effect since, as mentioned above, distortion can no longer be
            toggled in Roman PSFs. WFI PSFs natively include distortion effects.
        crop_psf : bool
            Included for API compatibility with the JWST instrument classes,
            but has no effect on the results for Roman WFI PSF calculations.

        """

        if add_distortion is not None:
            _log.warn('Note that the add_distortion argument no longer affects '
                      'Roman instrument simulations. All WFI PSFs natively '
                      'include distortion effects.')

        # Save new keyword to the options dictionary
        self.options['crop_psf'] = crop_psf

        # Run poppy calc_psf
        psf = stpsf_core.SpaceTelescopeInstrument.calc_psf(
            self,
            outfile=outfile,
            source=source,
            nlambda=nlambda,
            monochromatic=monochromatic,
            fov_arcsec=fov_arcsec,
            fov_pixels=fov_pixels,
            oversample=oversample,
            detector_oversample=detector_oversample,
            fft_oversample=fft_oversample,
            overwrite=overwrite,
            display=display,
            save_intermediates=save_intermediates,
            return_intermediates=return_intermediates,
            normalize=normalize,
        )

        return psf

    # slightly different versions of the following two functions from the parent
    # superclass in order to interface with the FieldDependentAberration class
    @property
    def detector_position(self):
        """The pixel position in (X, Y) on the detector"""
        return self._detector_position

    @detector_position.setter
    def detector_position(self, position):
        """Save new detector pixel position into the each detector's
        FieldDependentAberration instance. Update index of matching field
        point's slice in pupil file data cube."""
        try:
            x, y = map(int, position)
        except ValueError:
            raise ValueError('Detector pixel coordinates must be pairs of nonnegative numbers, ' 'not {}'.format(position))
        if x < 0 or y < 0:
            raise ValueError('Detector pixel coordinates must be nonnegative integers')
        if x > self._detector_npixels - 1 or y > self._detector_npixels - 1:
            raise ValueError(
                'The maximum allowed detector pixel ' 'coordinate value is {}'.format(self._detector_npixels - 1)
            )

        for det, aber_obj in self._detectors.items():
            aber_obj.field_position = (int(position[0]), int(position[1]))
        self._detector_position = position

        self._pupil_datacube_index = _wfi_sci_xy_to_fp(*position) - 1
        # subtract 1 because field point indices in the GSFC source data are
        # 1-indexed while the pupil file datacube is 0-indexed

    def _get_aberrations(self):
        """Get the OpticalElement that applies the field-dependent
        optical aberrations. (Called in get_optical_system.)"""
        return self._detectors[self._detector]

    def _get_fits_header(self, result, options):
        """Populate FITS Header keywords"""
        super()._get_fits_header(result, options)
        result[0].header['DETXPIXL'] = (self.detector_position[0], 'X pixel position (for field dependent aberrations)')
        result[0].header['DETYPIXL'] = (self.detector_position[1], 'Y pixel position (for field dependent aberrations)')
        result[0].header['DETECTOR'] = (self.detector, 'Detector selected')

    def _calc_psf_format_output(self, result, options):
        """
        Add detector charge diffusion model to the created 1-extension PSF.

        (As of STPSF v2.1.0, no longer applies distortion given that WFI optical
         models were delivered with inherent distortion effects in Cycle 10. If
         future deliveries arrive sans distortion, see pre-v2.1.0 commit history
         for code that handled distortion.)

        Apply desired formatting to output file:
                 - rebin to detector pixel scale if desired
                 - set up FITS extensions if desired
                 - output either the oversampled, rebinned, or both
        Which image(s) get output depends on the value of the options['output_mode']
        parameter. It may be set to 'Oversampled image' to output just the oversampled image,
        'Detector sampled image' to output just the image binned down onto detector pixels, or
        'Both as FITS extensions' to output the oversampled image as primary HDU and the
        rebinned image as the first image extension. For convenience, the option can be set
        to just 'oversampled', 'detector', or 'both'.

        Modifies the 'result' HDUList object.

        """
        # Set up new extensions for detector charge transfer model
        n_exts = len(result)
        for ext in np.arange(n_exts):
            hdu_new = fits.ImageHDU(result[ext].data, result[ext].header)
            # ^ these will be the PSFs that are edited ^
            result.append(hdu_new)
            ext_new = ext + n_exts
            result[ext_new].header['EXTNAME'] = result[ext].header['EXTNAME'][0:4] + 'DIST'  # change extension name
            _log.debug('Appending new extension {} with EXTNAME = {}'.format(ext_new, result[ext_new].header['EXTNAME']))

        # apply detector charge transfer model
        psf_updated = detectors.apply_detector_charge_diffusion(result,
                                                                options)

        # Edit the variable to match if input didn't request distortion
        # (cannot set result = psf_updated due to return method)
        [result.append(fits.ImageHDU()) for i in np.arange(len(psf_updated) - len(result))]
        for ext in np.arange(len(psf_updated)):
            result[ext] = psf_updated[ext]

        # Rewrite result variable based on output_mode set:
        stpsf_core.SpaceTelescopeInstrument._calc_psf_format_output(self, result, options)


class WFIPupilController:
    """
    This is a helper class for the WFI and is used to swap in
    the correct pupil each time the detector is changed.
    The pupil depends on pupil_mask, detector, and filter;
    pupil_mask is set automatically upon the receipt of the
    detector and filter values selected by the user in the WFI class.
    Users should only interact with this class through the API provided
    in the WFI class.

    Parameters
    ----------
    datapath : string
        Path to STPSF-WFI data files
    """

    def __init__(self, datapath):
        self.set_base_path(datapath)

        self._pupil = None
        self._pupil_mask = None

        # Flag to en-/disable automatic selection of the appropriate pupil_mask
        self._auto_pupil = True

        # Flag to en-/disable automatic selection of the appropriate pupil file
        self._auto_pupil_mask = True

    @property
    def pupil(self):
        """
        The path to the FITS file containing pupil information for the
        detector/filter combination sent from the WFI class. Cannot be
        directly set by the user.
        """
        return self._pupil

    @pupil.setter
    def pupil(self, value):
        raise AttributeError('Pupil cannot be directly specified. ' 'Use lock_pupil() instead.')

    @property
    def pupil_mask(self):
        """
        The corresponding mask for the filter sent from the WFI class.
        (See WFI.pupil_mask_list for a list of valid filters.) Cannot
        be directly set by the user.
        """
        return self._pupil_mask

    @pupil_mask.setter
    def pupil_mask(self, name):
        raise AttributeError('Pupil mask cannot be directly specified. '
                             'Use lock_pupil_mask() instead.')

    def _get_pupil_mask(self, wfi_filter):
        """
        Returns the appropriate pupil mask for a given WFI filter.

        Parameters
        ----------
        wfi_filter : string
            See WFI.filter_list for a list of valid filters.
        """
        wfi_filter = wfi_filter.upper()

        if wfi_filter in GRISM_FILTERS:
            # As of Cycle 10, GRISM0 and GRISM1 are the only filters that share
            # a set of pupil masks with each other
            return 'GRISM'
        else:
            return wfi_filter

    def set_base_path(self, datapath):
        """
        Sets the root directory of the path to STPSF's data files.
        This should be set before this class is used.

        Parameters
        ----------
        datapath : string
            Path to STPSF-WFI data files
        """
        self._datapath = datapath
        self._pupil_basepath = os.path.join(self._datapath, 'pupils')

    def pupil_file_formatter(self, pupil_mask, detector):
        """
        Generate proper pupil filename given a filter and a detector.

        Parameters
        ----------
        pupil_mask : string
            See WFI.filter_list for a list of valid filters. Pupil mask names
            match filter names for all except 'GRISM', which is shared between
            pupil masks 'GRISM0' and 'GRISM1'.

        detector : string
            See WFI.detector_list for a list of valid detectors.
        """
        return (f"RST_WFI_pupil_{pupil_mask.title()}_"
                f"{detector}_allfieldpoints.fits.gz")

    def update_pupil(self, wfi_filter, detector):
        """
        Selects the specific pupil file corresponding with a detector
        and filter combination sent from the WFI class. Also finds and
        indirectly sets the proper pupil_mask in the process.

        Parameters
        ----------
        wfi_filter : string
            See WFI.filter_list for a list of valid filters.

        detector : string
            See WFI.detector_list for a list of valid detectors.
        """
        if not self._auto_pupil:
            _log.info('Automatic pupil selection was locked; '
                      'using user-provided pupil.')
            return

        if self._pupil_basepath is None:
            raise Exception('update_pupil called before setting pupil file path')

        # figure out proper mask based on filter (or use locked mask if enabled)
        pupil_mask = self._get_pupil_mask(wfi_filter) if self._auto_pupil_mask else self.pupil_mask
        pupil = os.path.join(self._pupil_basepath,
                             self.pupil_file_formatter(pupil_mask, detector))

        self._pupil_mask = pupil_mask
        self._pupil = pupil

        _log.info(
            f"Using {'' if self._auto_pupil_mask else 'locked '}"
            f"pupil mask '{pupil_mask}' and detector '{detector}'."
        )

    def lock_pupil(self, pupil_path):
        """
        Prevents the WFIPupilController class from dynamically updating
        the path to the pupil on any changes to the detector or filter
        selected in the WFI class. Instead, the path remains locked on
        whichever `pupil_path` was provided to this method.

        CAUTION: This is non-standard usage of the WFI class and may
        lead to unexpected behavior.

        Parameters
        ----------
        pupil_path : string
            The custom path to your pupil file.
        """
        self._pupil_mask = None
        self._pupil = pupil_path
        self._auto_pupil = False

    def unlock_pupil(self):
        """
        Undoes the effects of lock_pupil() and resets WFIPupilController
        to its default state of updating the pupil whenever a detector
        or filter is changed in the WFI class.
        """
        self._auto_pupil = True

    def lock_pupil_mask(self, pupil_mask):
        """
        Prevents the WFIPupilController class from dynamically updating
        the pupil mask on any changes to the filter selected in the WFI
        class. Instead, the pupil mask remains locked on whichever
        `pupil_mask` was provided to this method.

        CAUTION: This is non-standard usage of the WFI class and may
        lead to unexpected behavior.

        Parameters
        ----------
        filter : string
            See WFI.pupil_mask_list for a list of valid pupil masks.
        """
        if not self._auto_pupil:
            raise Exception('Pupil is locked. Unlock pupil before locking pupil mask.')
        else:
            self._pupil_mask = pupil_mask
            self._auto_pupil_mask = False

    def unlock_pupil_mask(self):
        """
        Undoes the effects of lock_pupil_mask() and resets
        WFIPupilController to its default state of updating the pupil
        mask whenever filter is changed in the WFI class.
        """
        self._auto_pupil_mask = True


class WFI(RomanInstrument):
    """
    WFI represents the Roman mission's Wide Field Imager.

    WARNING: This model has not yet been validated against other PSF
             simulations and uses several approximations.
    """
    NPIXELS = 4096

    def __init__(self):
        # pixel scale is from Roman-AFTA SDT report final version (p. 91)
        # https://roman.ipac.caltech.edu/sims/Param_db.html
        pixelscale = 110e-3  # arcsec/px

        # Initialize the aberrations for parent SpaceTelescopeInstrument class
        self._aberration_files = {}
        self._is_custom_aberration = False
        self._current_aberration_file = ''

        super().__init__('WFI', pixelscale=pixelscale)

        # Initialize the pupil controller
        self._pupil_controller = WFIPupilController(self._datapath)

        self.pupil_mask_list = [fltr for fltr in self.filter_list.copy()
                                if not fltr.startswith('GRISM')]
        self.pupil_mask_list.append(np.str_('GRISM'))  # GRISM0/1 share pupils

        # Define default aberration files for WFI filters
        self._aberration_files = {'custom': None}
        self._aberration_files.update({
            fltr: os.path.join(
                self._datapath,
                'aberrations',
                # f"{fltr}_25fields_Z210_multiwave.csv"
                f"{fltr}_25fields_zernike_Z210_multiwavelength.csv"
                if fltr.startswith('F')
                # else f"{fltr}_5fields_Z45_multiwave.csv")
                else f"{fltr}_5fields_zernike_Z45_multiwavelength.csv")
            for fltr in self.filter_list})

        # Load aberration info from ref files
        self._detector_npixels = self.NPIXELS
        self._load_detector_aberrations(self._aberration_files[self.mode])

        self._opd_file_dict = {
            det: os.path.join(
                self._datapath,
                'aberrations',
                # f"{det}_Z210_high_freq_cube.fits"
                f"RST_WFI{i}_Z210_high_frequency_map.fits")
            # for det in self.detector_list})
            for i, det in enumerate(self.detector_list, 1)}
        self._pupilopd = self._opd_file_dict['WFI01']

        # Set initial detector and position
        self.detector = 'WFI01'
        self.detector_position = (self.NPIXELS // 2, self.NPIXELS // 2)

    def _addAdditionalOptics(self, optsys, **kwargs):
        _log.debug('   No optics added for WFI')
        return optsys, False, None

    def _load_detector_aberrations(self, path):
        """
        Helper function that, given a path to a file containing detector
        aberrations, loads the Zernike values and populates the class'
        dictator list with `FieldDependentAberration` detectors. This
        function achieves this by calling the
        `stpsf.roman._load_wfi_detector_aberrations` function.

        Users should use the `override_aberrations` function to override
        current aberrations.

        Parameters
        ----------
        path : string
            The path to the file containing detector aberrations.
        """
        detectors_dict = _load_wfi_detector_aberrations(path)
        assert len(detectors_dict.keys()) > 0

        self._detectors = detectors_dict
        self._current_aberration_file = path

    def _validate_config(self, **kwargs):
        """
        Validates that the WFI is configured sensibly.

        This mainly consists of selecting the masked or unmasked pupil
        appropriately based on the wavelengths requested.
        """
        assert self.filter is not None, 'filter is None'
        assert self.detector is not None, 'detector is None'
        self._update_pupil()

        assert self.pupil is not None, 'pupil is None'
        super()._validate_config(**kwargs)

    @property
    def pupilopd(self):
        """The file containing high-frequency Zernike information for the
        current detector. Set by the detector setter."""
        return self._pupilopd

    @pupilopd.setter
    def pupilopd(self, value):
        # Only allow direct set of pupilopd when done by parent class in super()
        # (i.e., before WFI class has set a detector)
        if self.detector is None:
            self._pupilopd = value
        else:
            raise ValueError('WFI.pupilopd is set automatically on updates to '
                             'its detector attribute, not manually.')

    def _update_pupil(self, wfi_filter=None, detector=None):
        """
        Chooses proper pupil file. Pupil file geometry depends on field position
        parameterized by detector and field position number within a detector.
        """
        if detector is None:
            detector = self.detector
        if wfi_filter is None:
            wfi_filter = self.filter

        if detector is not None and wfi_filter is not None:
            self._pupil_controller.update_pupil(wfi_filter=wfi_filter,
                                                detector=detector)

    @RomanInstrument.detector.setter
    def detector(self, value):
        """
        The current WFI detector. See WFI.detector_list for valid values.

        Also adjusts the current 1) pupil file (since these depend on filter and
        *detector*), 2) path to throughput files, and 3) OPD file containing
        contributions from high-frequency Zernike terms above Z210.
        """
        if value.upper().startswith('SCA'):  # backward-compatible name assignment
            value = f"WFI{value[-2:]}"
        if value.upper() not in self.detector_list:
            raise ValueError('Invalid detector. Valid detector names are: {}'.format(', '.join(self.detector_list)))

        self._detector = value.upper()

        # Update pupil and high-frequency OPD file locations (conditional needed
        # to exclude call to setter by SpaceTelescopeInstrument.__init__()
        if self._detector is not None:
            self._update_pupil(detector=self._detector)

        # Update throughput file directory
        self._update_aperturename()

    def _update_aperturename(self):
        """Update SIAF aperture name after change in detector or other relevant properties.
        This function handles just inferring the new value of the aperturename.
        See the aperturename.setter function for additional changed based on that.
        """
        self.aperturename = self._detector.replace('SCA', 'WFI') + "_FULL"

    @RomanInstrument.aperturename.setter
    def aperturename(self, value):
        """Set SIAF aperture name to new value, with validation.

        This also updates the pixelscale to the local value for that aperture, for a small precision enhancement.
        """
        try:
            ap = self.siaf[value]
        except KeyError:
            raise ValueError(f'Aperture name {value} not a valid SIAF aperture name for {self.name}')

        # Only update if new value is different
        if self._aperturename != value:
            self._aperturename = value

            # Update pixelscale based on specified aperture name
            self.pixelscale = self._get_pixelscale_from_apername(value)

        # adjust filter throughput file location
        found_parent = False
        filters_info = {}
        for fl, info in self._filters.items():
            filter_path = Path(info.filename)

            # find 'filters', the earliest common parent directory of default
            # stpsf filter file structure (stpsf-data/INSTRUMENT/*filters*)
            # and WFI's filter file structure (stpsf-data/WFI/*filters*/DETECTOR).
            if not found_parent:
                for parent in filter_path.parents:
                    if parent.name == 'filters':
                        found_parent = True
                        break
                    else:
                        continue
            # done this way because SpaceTelescopeInstrument.__init__()
            # (a WFI parent class) sets filter location to the default WFI path
            # with SpaceTelescopeInstrument._get_filters() upon class creation,
            # so the filter files' direct parent directory the first time a
            # detector is assigned in WFI.__init__() (the "filters" directory)
            # is not at the same level as on later changes to the detector attribute
            # (the corresponding "WFINN" child directory of "WFI/filters").

            # save filter info but change throughput file based on new detector
            info_edit = {key : (val if key != 'filename'
                                else str(parent / self._detector / filter_path.name))
                         for key, val in info._asdict().items()}
            filters_info[fl] = stpsf_core.Filter(**info_edit)

        self._filters = filters_info

    @RomanInstrument.filter.setter
    def filter(self, value):
        """
        The current WFI filter. See WFI.filter_list for valid values.
        """
        # Update filter
        value = value.upper()

        if value not in self.filter_list:
            raise ValueError(f"Instrument {self.name} doesn't have a "
                             f"filter called {value}.")

        self._filter = value

        # Update aberrations if self._aberration_files has been initiated (not
        # empty) and if they haven't been locked by user
        if self._aberration_files and not self._is_custom_aberration:
            # identify aberration file for new mode
            aberration_file = self._aberration_files[self._filter]

            # if aberrations are not already loaded for the new mode,
            # load and replace detectors using the new mode's aberration file
            if not os.path.samefile(self._current_aberration_file, aberration_file):
                self._load_detector_aberrations(aberration_file)

        # Update pupil only if detector was previously loaded
        # ( i.e., skip this step when called by super() )
        if self.detector is not None:
            self._update_pupil(wfi_filter=self._filter)

    @property
    def pupil(self):
        """
        The path to the FITS file containing pupil information for the
        detector/filter combination sent from the WFI class. Cannot be
        directly set by the user.
        """
        return self._pupil_controller.pupil

    @pupil.setter
    def pupil(self, value):
        # don't allow pupil to be set until the pupil controller is active. (a
        # parent class tries to set it to None in WFI's preceding super() call)
        if hasattr(self, '_pupil_controller'):
            raise AttributeError('Pupil cannot be directly specified. ' 'Use lock_pupil() instead.')

    @property
    def pupil_mask(self):
        """
        The corresponding mask for the current filter. Cannot be
        directly set by the user.
        """
        return self._pupil_controller.pupil_mask

    @pupil_mask.setter
    def pupil_mask(self, name):
        raise AttributeError('Pupil mask cannot be directly specified. ' 'Use lock_pupil_mask() instead.')

    @property
    def mode(self):
        """
        The current WFI mode. Cannot be directly set by the user.
        """
        return self.filter

    @mode.setter
    def mode(self, value):
        raise AttributeError('WFI mode cannot be directly specified; ' 'it is set by changing filters.')

    def lock_aberrations(self, aberration_path):
        """
        This function loads user provided aberrations from a file and
        locks this instrument to only use the provided aberrations (even
        if the filter or mode change).

        To release the lock and load the default aberrations, use
        unlock_aberrations(). To load new user provided
        aberrations, call this function with the new path.

        To load custom aberrations, please provide a csv file
        containing the detector names, field point positions and Zernike
        values. The file should contain the following column
        names/values (comments in parentheses should not be included):
        - sca (Detector number)
        - wavelength (µm)
        - field_point (field point number/ID for detector and wavelength.
                       starts with 1)
        - local_x (mm, local detector coords)
        - local_y (mm, local detector coords)
        - global_x (mm, global instrument coords)
        - global_y (mm, global instrument coords)
        - axis_local_angle_x (XAN)
        - axis_local_angle_y (YAN)
        - wfe_rms_waves (nm)
        - wfe_pv_waves (waves)
        - Z1 (Zernike phase NOLL coefficients)
        - Z2 (Zernike phase NOLL coefficients)
        - Z3 (Zernike phase NOLL coefficients)
        - Z4 (Zernike phase NOLL coefficients)
        - ...

        Please refer to the default aberration files for examples. If
        you have the STPSF data installed and defined, you can get the
        path to that file by running the following:
        >>> from stpsf import roman
        >>> wfi = roman.WFI()
        >>> print(wfi._aberration_files['imaging'])

        Warning: You should not edit the default files!
        """
        self._load_detector_aberrations(aberration_path)
        self._aberration_files['custom'] = aberration_path
        self._is_custom_aberration = True

    def unlock_aberrations(self):
        """
        Releases the lock on the detector aberration file location
        and loads the default file.
        """
        aberration_path = self._aberration_files[self.mode]
        self._load_detector_aberrations(aberration_path)
        self._aberration_files['custom'] = None
        self._is_custom_aberration = False

    def lock_pupil(self, pupil_path):
        """
        Prevents dynamic updates of the path to the proper pupil file on
        any changes to the selected detector or filter. Instead, the
        path remains locked on whichever `pupil_path` was provided here.

        WARNING: This is non-standard usage of the WFI class and may
        lead to unexpected behavior.

        Parameters
        ----------
        pupil_path : string
            The custom path to your pupil file.
        """
        if os.path.isfile(pupil_path):
            self._pupil_controller.lock_pupil(pupil_path)
        else:
            raise FileNotFoundError(f'{pupil_path} not found.')

        _log.warning('Disabling default pupil selection behavior.')

    def unlock_pupil(self):
        """
        Undoes the effects of lock_pupil() by resetting the class to
        its default state of updating the pupil whenever a detector or
        filter is changed. If necessary, it also sets the proper pupil
        for the current detector/filter combination.
        """
        self._pupil_controller.unlock_pupil()
        self._update_pupil()  # reset pupil
        _log.info('Restoring default pupil selection behavior.')

    def lock_pupil_mask(self, pupil_mask):
        """
        Prevents dynamic updates of the pupil mask on any change to the
        selected filter. Instead, the pupil mask remains locked on
        whichever `pupil_mask` was provided here.

        WARNING: This is non-standard usage of the WFI class and may
        lead to unexpected behavior.

        Parameters
        ----------
        filter : string
            See WFI.pupil_mask_list for a list of valid pupil masks.
        """
        if pupil_mask not in self.pupil_mask_list:
            raise Exception('invalid pupil mask')
        self._pupil_controller.lock_pupil_mask(pupil_mask)
        self._update_pupil()
        _log.warning('Disabling default pupil mask selection behavior.')

    def unlock_pupil_mask(self):
        """
        Undoes the effects of lock_pupil_mask() and resets the class to
        its default state of updating the pupil mask whenever the filter
        is changed.
        """
        self._pupil_controller.unlock_pupil_mask()
        self._update_pupil()  # reset pupil mask
        _log.info('Restoring default pupil mask selection behavior.')
