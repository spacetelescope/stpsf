import numpy as np

from .. import webbpsf_core
from .test_webbpsf import do_test_set_position_from_siaf, do_test_source_offset, generic_output_test

# ------------------    NIRISS Tests    ----------------------------


def test_niriss():
    return generic_output_test('NIRISS')


def test_niriss_source_offset_00():
    return do_test_source_offset('NIRISS', theta=0.0, monochromatic=3.0e-6)


def test_niriss_source_offset_45():
    return do_test_source_offset('NIRISS', theta=45.0, monochromatic=3.0e-6)


def test_niriss_set_siaf():
    return do_test_set_position_from_siaf(
        'NIRISS',
        ['NIS_FP1MIMF', 'NIS_SUB64', 'NIS_SOSSFULL', 'NIS_SOSSTA', 'NIS_AMI1']
    )


def test_niriss_auto_pupil():
    """Test switching between CLEAR and CLEARP
    depending on selected filter or wavelengths
    """

    niriss = webbpsf_core.NIRISS()
    assert niriss.pupil_mask is None

    niriss.filter = 'F277W'
    niriss.calc_psf(nlambda=1)
    assert niriss.pupil_mask == 'CLEARP'

    niriss.filter = 'F090W'
    niriss.calc_psf(nlambda=1)
    assert niriss.pupil_mask is None

    niriss.filter = 'F480M'
    niriss.calc_psf(nlambda=1)
    assert niriss.pupil_mask == 'CLEARP'

    niriss.filter = 'F200W'
    niriss.calc_psf(nlambda=1)
    assert niriss.pupil_mask is None


def test_niriss_gr700xd():
    """
    Smoke-test calculations with the GR700XD custom optic
    present in the system. This is a regression test for
    https://github.com/mperrin/webbpsf/issues/148
    """
    niriss = webbpsf_core.NIRISS()
    niriss.filter = 'CLEAR'
    niriss.pupil_mask = 'GR700XD'
    niriss.calc_psf(monochromatic=1e-6, fov_pixels=2)


def test_niriss_aperturename():
    """Not a lot of options here"""
    niriss = webbpsf_core.NIRISS()
    assert niriss.aperturename == niriss._detectors[niriss.detector], 'Default SIAF aperture is not as expected'

    ref_tel_coords = niriss._tel_coords()

    niriss.aperturename = 'NIS_SUB128'
    assert niriss.detector_position == (64, 64), (
        "Changing to a subarray aperture didn't change the " 'reference pixel coords as expected'
    )
    assert np.any(niriss._tel_coords() != ref_tel_coords), (
        "Changing to a subarray aperture didn't change the V2V3 coords " 'as expected.'
    )
