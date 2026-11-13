"""
===============
HST Instruments
===============

A re-implementation of the TinyTim optical model for the Hubble Space Telescope, built on poppy and
the STPSF instrument framework. Currently this covers WFC3 (IR and UVIS channels).

WARNING: Work in progress. The goal of this first phase is to reproduce TinyTim's outputs.
"""

import logging
import os
from functools import lru_cache

import astropy.io.fits as fits
import astropy.units as u
import numpy as np
from astropy.io import ascii as ioascii

from . import constants, stpsf_core, utils
from .distortion import _siaf_mean_pixelscale

_log = logging.getLogger('stpsf')


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
    """

    telescope = 'HST'
    _default_detector = None  # set in subclasses; used to choose the initial pixel scale

    def __init__(self, name, pixelscale=None):
        self.siaf = stpsf_core.get_siaf_with_caching('HST')
        if pixelscale is None:
            pixelscale = self._default_pixelscale(self._default_detector)

        super().__init__(name, pixelscale=pixelscale)
        self.pupil_radius = constants.HST_PUPIL_DIAMETER / 2 * u.meter
        self._aperturename = None

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
