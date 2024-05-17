import logging

import astropy.units as u
import numpy as np
import pysiaf

from .. import webbpsf_core
from .test_webbpsf import do_test_set_position_from_siaf, do_test_source_offset, generic_output_test

_log = logging.getLogger('test_webbpsf')
_log.addHandler(logging.NullHandler())


# ------------------    NIRSpec Tests    ----------------------------

def test_nirspec():
    return generic_output_test('NIRSpec')
# Use a larger than typical tolerance when testing NIRSpec offsets. The
# pixels are so undersampled (0.1 arcsec!) that it's unreasonable to try for
# better than 1/10th of a pixel precision using default settings.
def test_nirspec_source_offset_00():
    return do_test_source_offset('NIRSpec', theta=0.0, tolerance=0.1, monochromatic=3e-6)

def test_nirspec_source_offset_45():
    return do_test_source_offset('NIRSpec', theta=45.0, tolerance=0.1, monochromatic=3e-6)

def test_nirspec_set_siaf():
    return do_test_set_position_from_siaf('NIRSpec')


def test_nirspec_slit_apertures():
    """Test that we can use slit and aperture names that don't map to a specific detector
    Verify that the V2 and V3 coordinates are reported as expected.
    """
    nrs = webbpsf_core.NIRSpec()

    for apname in ['NRS_FULL_IFU', 'NRS_S200A1_SLIT']:
        nrs.set_position_from_aperture_name(apname)

        ap = pysiaf.Siaf('NIRSpec')[apname]

        assert np.isclose(nrs._tel_coords()[0].to_value(u.arcsec), ap.V2Ref)
        assert np.isclose(nrs._tel_coords()[1].to_value(u.arcsec), ap.V3Ref)


def test_calc_datacube_fast():
    nrs = webbpsf_core.NIRSpec()
    nrs.set_position_from_aperture_name('NRS_FULL_IFU')
    nrs.image_mask = 'IFU'

    waves = np.linspace(3e-6, 5e-6, 3)

    nrs.calc_datacube_fast(waves, fov_pixels=30, oversample=1, compare_methods=True)  # TODO assert success
    


def test_mode_switch():
    """ Test switch between IFU and imaging modes """
    nrs = webbpsf_core.NIRSpec()
    # check mode swith to IFU
    nrs.mode = 'IFU'
    assert 'IFU' in nrs.aperturename
    assert nrs.band == 'PRISM/CLEAR'

    # check switch of which IFU band
    nrs.disperser = 'G395H'
    nrs.filter = 'F290LP'
    assert nrs.band == 'G395H/F290LP'

    # check mode switch back to imaging
    nrs.mode = 'imaging'
    assert 'IFU' not in nrs.aperturename

def test_IFU_wavelengths():
    """ Test computing the wqvelength sampling for a sim IFU cube """
    nrs = webbpsf_core.NIRSpec()
    # check mode swith to IFU
    nrs.mode = 'IFU'
    nrs.disperser = 'G235H'
    nrs.filter = 'F170LP'
    waves = nrs.get_IFU_wavelengths()
    assert isinstance(waves, u.Quantity)

    assert len(waves) > 3000  # there are lots of wavelengths in the high-resolution grating cubes
    # and test we can specify a reduced wavelength sampling:
    for n in (10, 100):
        assert len(nrs.get_IFU_wavelengths(n)) == n
