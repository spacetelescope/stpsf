import os

import numpy as np
import pytest
from astropy.table import Table
from numpy import allclose

from stpsf import measure_fwhm, roman

GRISM_FILTERS = roman.GRISM_FILTERS
PRISM_FILTERS = roman.PRISM_FILTERS


def pupil_path(wfi):
    """
    dynamically generate current pupil path for a given WFI instance
    """
    element = 'GRISM' if wfi.filter in GRISM_FILTERS else wfi.filter

    base = wfi._datapath
    file = f'pupils/RST_WFI_pupil_{element}_WFI{wfi.detector[-2:]}_allfieldpoints.fits.gz'

    return os.path.join(base, file)


def test_WFI_psf():
    """
    Test that instantiating WFI works and can compute a PSF without
    raising any exceptions
    """
    wfi = roman.WFI()
    wfi.calc_psf(fov_pixels=4)


def test_WFI_filters():
    wfi = roman.WFI()
    filter_list = wfi.filter_list

    for fltr in filter_list:
        wfi.filter = fltr
        wfi.calc_psf(fov_pixels=4, oversample=1, nlambda=3)


def test_aberration_detector_position_setter():
    detector = roman.FieldDependentAberration(4096, 4096)

    with pytest.raises(ValueError) as excinfo:
        detector.field_position = (-1, 1)
    assert 'pixel_x' in str(excinfo.value), 'Failed to raise exception for small out-of-bounds ' 'x pixel position'
    with pytest.raises(ValueError) as excinfo:
        detector.field_position = (4096 + 1, 1)
    assert 'pixel_x' in str(excinfo.value), 'Failed to raise exception for large out-of-bounds ' 'x pixel position'
    with pytest.raises(ValueError) as excinfo:
        detector.field_position = (1, -1)
    assert 'pixel_y' in str(excinfo.value), 'Failed to raise exception for small out-of-bounds ' 'y pixel position'
    with pytest.raises(ValueError) as excinfo:
        detector.field_position = (1, 4096 + 1)
    assert 'pixel_y' in str(excinfo.value), 'Failed to raise exception for large out-of-bounds ' 'y pixel position'

    valid_pos = (1.0, 1.0)
    detector.field_position = valid_pos
    assert detector._field_position == valid_pos, (
        'Setting field position through setter did not ' 'update private `_field_position` value'
    )


def test_WFI_fwhm():
    """
    Test that computed PSFs are physically realistic, at least relatively.
    Loose test...
    """
    wfi = roman.WFI()

    wfi.pupilopd = None
    wfi.options['jitter'] = None

    wfi.filter = 'F062'
    fwhm_f062 = measure_fwhm(wfi.calc_psf(oversample=6))

    wfi.filter = 'F184'
    fwhm_f184 = measure_fwhm(wfi.calc_psf(oversample=6))

    assert 4.0 > fwhm_f184 / fwhm_f062 > 2.0


def test_wfi_pupil_field_dependence():
    """Test that it automatically loads different pupil based on detector and detector position
     - The filename is selected based on detector and filter.
     - The index into the datacube in that file is selected based on detector position
    The expectation is there are a total of 25 field points (5x5) across each detector, with a particular
    numbering and arrangement as used in the 2026 prelaunch data delivery.
    Furthermore there are distinct sets of such files for each filter.

    This test only iterates over a subset of filters and detectors and positions, because that
    suffices to test the relevant functionalities in general.

    """

    elements = ['F087', "F184", 'PRISM']
    dets = ['WFI01', 'WFI09', 'WFI18']
    corners = np.asarray( [(0,0), (0, 1), (1,0), (1,1)]) * 4090
    corner_field_point_expected = [21, 25, 1, 5]

    for element in elements:
        for det in dets:
            wfi = roman.WFI()
            wfi.filter = element
            wfi.detector = det
            assert wfi.detector_position == (2048, 2048), "Default detector position should be middle of the detector"
            assert wfi.pupil_datacube_index == 13, "Default pupil datacube index should be middle of the 25 field points"

            for (x,y), fp  in zip(corners, corner_field_point_expected):
                wfi.detector_position = (x, y)

                # Check that the wfi.pupil element updates as expected
                assert det in wfi.pupil, "The pupil filename ought to update based on detector "
                assert element in wfi.pupil, "The pupil filename ought to update based on filter"
                assert wfi.pupil_datacube_index == fp, ("The pupil file datacube index did not match the"+
                                                        f" expected field point: got {wfi.pupil_datacube_index} instead of {fp}")

                # Check that if we toggle off the auto_pupil feature then the pupil and datacube index stop changing.
                selected_pupil_file = wfi.pupil
                selected_pupil_cube_index = wfi.pupil_datacube_index
                wfi.auto_pupil = False
                wfi.detector_position = (512, 1536)
                assert wfi.pupil_datacube_index == selected_pupil_cube_index, "Pupil file should not change with new det pos if auto_pupil is False"
                assert wfi.pupil == selected_pupil_file, "Pupil file should not change with new det pos"
                wfi.detector = 'WFI03'
                assert wfi.pupil == selected_pupil_file, "Pupil file should not change with new detector if auto_pupil is False"
                wfi.detector = det   # Reset before continuing the next part of the test
                wfi.detector_position = (x, y)

                # now re-enable auto_pupil and verify that it works as expected,
                #   first for a change in detector position, then for a change in detector
                wfi.auto_pupil = True
                wfi.detector_position = (512, 1536)
                assert wfi.pupil_datacube_index != selected_pupil_cube_index, "Pupil file should change with new det pos if auto_pupil is True"
                assert wfi.pupil == selected_pupil_file, "Pupil file should not change with new det pos "
                wfi.detector = 'WFI03'
                assert wfi.pupil != selected_pupil_file, "Pupil file should change with new detector if auto_pupil is True"
                wfi.detector = det   # Reset to current detector before continuing the for loop over field positions




def test_WFI_detector_position_setter():
    wfi = roman.WFI()
    wfi.detector = 'WFI01'
    valid_pos = (4000, 1000)
    wfi.detector_position = valid_pos
    assert wfi._detectors[wfi._detector].field_position == valid_pos, (
        'Setting field position through Instrument.detector_position did not update field_position '
        "for the detector's aberration optic"
    )
    assert wfi.detector_position == valid_pos, "`detector_position` getter doesn't reflect " 'assignment to setter'


def test_WFI_includes_aberrations():
    wfi = roman.WFI()
    wfi.detector = 'WFI01'
    osys = wfi.get_optical_system()
    assert isinstance(osys[2], roman.FieldDependentAberration), (
        'Third plane of Roman WFI optical system should be the ' 'field dependent aberration virtual optic'
    )


def test_swapping_modes(wfi=None):
    if wfi is None:
        wfi = roman.WFI()

    tests = [
        # [filter, mode, pupil_file]
        ['F146', 'F146', pupil_path],
        ['F213', 'F213', pupil_path],
        [PRISM_FILTERS[0], PRISM_FILTERS[0], pupil_path],
        [GRISM_FILTERS[0], GRISM_FILTERS[0], pupil_path],
    ]

    for test_filter, test_mode, test_pupil in tests:
        wfi.filter = test_filter

        fail_str = f"failed on {test_filter}, {test_mode}, " f"{test_pupil(wfi).split('/')[-1]}"

        assert wfi.filter == test_filter, fail_str
        assert wfi.mode == test_mode, fail_str
        assert wfi._current_aberration_file == wfi._aberration_files[test_mode], fail_str
        assert wfi.pupil == test_pupil(wfi), fail_str


def test_custom_aberrations():
    wfi = roman.WFI()

    # Use GRISM0 aberration_file for testing
    test_aberration_file = wfi._aberration_files['GRISM0']

    # Test override
    # -------------
    wfi.lock_aberrations(test_aberration_file)

    for fltr in wfi.filter_list:
        wfi.filter = fltr
        assert wfi._current_aberration_file == test_aberration_file, 'Filter change caused override to fail'

    # Test Release Override
    # ---------------------
    wfi.unlock_aberrations()
    assert wfi._aberration_files['custom'] is None, 'Custom aberration file not deleted on override release.'
    test_swapping_modes(wfi)


def test_WFI_limits_interpolation_range():
    wfi = roman.WFI()
    det = wfi._detectors['WFI01']
    det.get_aberration_terms(1.29e-6)
    det.field_position = (0, 0)
    det.get_aberration_terms(1.29e-6)

    with pytest.raises(ValueError) as excinfo:
        det.field_position = (500000, 0)
    assert 'Requested pixel_x position' in str(
        excinfo.value
    ), 'FieldDependentAberration did not error on out-of-bounds field point'

    with pytest.raises(ValueError) as excinfo:
        det.field_position = (-1, 0)
    assert 'Requested pixel_x position' in str(
        excinfo.value
    ), 'FieldDependentAberration did not error on out-of-bounds field point'

    with pytest.raises(ValueError) as excinfo:
        det.field_position = (0, 500000)
    assert 'Requested pixel_y position' in str(
        excinfo.value
    ), 'FieldDependentAberration did not error on out-of-bounds field point'

    with pytest.raises(ValueError) as excinfo:
        det.field_position = (0, -1)
    assert 'Requested pixel_y position' in str(
        excinfo.value
    ), 'FieldDependentAberration did not error on out-of-bounds field point'

    det.field_position = (2048, 2048)

    # Get min and max valid wavelengths from aberration file
    zern = Table.read(wfi._aberration_files[wfi.mode], format='ascii.csv')
    min_wv = zern['wavelength'][0] * 1e-6  # convert from micron to meter
    max_wv = zern['wavelength'][-1] * 1e-6

    # Test that get_aberration_terms() uses an approximated wavelength when
    # called with an out-of-bounds wavelength.
    too_lo_wv = min_wv * 0.9
    too_hi_wv = max_wv / 0.9
    valid_wv = np.mean([min_wv, max_wv])

    assert allclose(
        det.get_aberration_terms(min_wv), det.get_aberration_terms(too_lo_wv)
    ), 'Aberration below wavelength range did not return closest value.'

    assert allclose(
        det.get_aberration_terms(max_wv), det.get_aberration_terms(too_hi_wv)
    ), 'Aberration above wavelength range did not return closest value.'

    # Test border pixels outside the ref data. In Cycle 10, (32, 0) is the first
    # pixel, so we check if (0, 0) is approximated to it as the nearest point.
    det.field_position = (0, 0)
    coefficients_outlier = det.get_aberration_terms(valid_wv)

    det.field_position = (32, 0)
    coefficients_data = det.get_aberration_terms(valid_wv)

    assert np.allclose(coefficients_outlier, coefficients_data), (
        'nearest point extrapolation ' 'failed for outlier field point'
    )


def test_WFI_auto_aperturename_and_pixelscale():
    wfi = roman.WFI()
    assert wfi.aperturename == "WFI01_FULL", "Aperture name should match detector"
    wfi.detector = 'WFI12'
    assert wfi.aperturename == "WFI12_FULL", "Aperture name should match detector"


    aperture = wfi.siaf[wfi.aperturename]
    assert wfi.pixelscale == (aperture.XSciScale + aperture.YSciScale)/2, "Pixel scale should match the SIAF for that aperturename"


# -------- Test functions for the (very limited) Roman Coronagraph implementation below here


def test_coronagraph_detector_position():
    """Test existence of the Coronagraph detector position etc, and that you can't set it."""
    cor = roman.RomanCoronagraph()

    valid_pos = (512, 512)
    assert cor.detector_position == valid_pos, "Coronagraph detector position isn't as expected"

    with pytest.raises(RuntimeError) as excinfo:
        cor.detector_position = valid_pos
    assert 'not adjustable' in str(excinfo.value), (
        'Failed to raise exception for' 'trying to change Coronagraph detector position.'
    )


def test_coronagraph_psf(display=False):
    """
    Just test that instantiating RomanCoronagraph works and can compute a PSF
    without raising any exceptions
    """
    char_spc = roman.RomanCoronagraph()
    char_spc.mode = 'CHARSPC_F660'

    # print('Reading instrument data from {:s}'.format(charspc._STPSF_basepath)
    # print('Filter list: {:}'.format(charspc.filter_list))

    monopsf = char_spc.calc_psf(nlambda=1, display=False)
    if display:
        roman.poppy.display_psf(monopsf)

# -------- Test utility functions for WFI field of view position coordinate conversions


def test_sci_xy_to_fp():
    """Test corners match expectations"""
    corner_data = ( ((0,0), 21), # lower left
                    ((0, 4096), 25), # upper left
                    ((4096, 0), 1), # lower right
                    ((4096, 4096), 5) # upper right
                  )
    for (x, y), fp in corner_data:
        assert roman._wfi_sci_xy_to_fp(x, y) == fp  # test xy->fp
        assert roman._wfi_sci_xy_to_fp( *roman._wfi_fp_to_sci_xy(fp)) == fp , f'round trip error for fp {fp}'  # test round trip
