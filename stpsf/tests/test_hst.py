"""Tests for the HST / WFC3 models. Those that need the HST reference data skip if it is missing from STPSF_PATH.

Comparisons with TinyTim's own output need a TinyTim installation, and are in
``dev_utils/run_hst_regression_tests.py`` instead.
"""

import numpy as np
import poppy
import pytest

from .. import hst
import stpsf


def _pupil_map(cls=hst.HSTPupil, **kwargs):
    # 512 pixel grid with the pupil radius spanning 128 pixels, as in TinyTim
    wave = poppy.Wavefront(diam=4.8, npix=512)
    wave.pixelscale = 1.2 / 128
    return cls(**kwargs).get_transmission(wave)


def test_hst_pupil_geometry():
    """The rotated OTA pupil has pads at 346 and 104 degrees, and spiders on the diagonals (from TinyTim)"""
    pupil = _pupil_map()
    center = 255.5
    for angle_deg in (346.4, 104.0):
        x = center + 0.886 * 128 * np.cos(np.deg2rad(angle_deg))
        y = center + 0.886 * 128 * np.sin(np.deg2rad(angle_deg))
        assert pupil[int(round(y)), int(round(x))] < 0.05, f'expected a pad at {angle_deg} degrees'
    # a point on the +x axis at mid radius is clear, but a point on a diagonal is in a spider
    assert pupil[256, int(center + 0.6 * 128)] > 0.99
    assert pupil[int(center + 0.6 * 128 * 0.7071), int(center + 0.6 * 128 * 0.7071)] < 0.1


def test_hst_pupil_exact_coverage():
    """Edge pixels have exact fractional coverage, so the area is independent of how the aperture sits on the grid"""
    expected = (1 - 0.33**2) - 4 * (1 - 0.33) * 0.022 / np.pi - 3 * 0.065**2  # annulus, spiders, pads
    areas = []
    for npix in (255, 256, 400, 1024):
        diam = 2.4 * npix / 256 * 1.1
        wave = poppy.Wavefront(diam=diam, npix=npix)
        areas.append(hst.HSTPupil().get_transmission(wave).sum() * (diam / npix) ** 2 / (np.pi * 1.2**2))
    assert np.allclose(areas, expected, atol=1e-3)


def test_hst_pupil_area():
    """The OTA pupil times the IR cold mask matches the TinyTim pupil dump area (39552 pixels at R=128)"""
    ota = _pupil_map()
    cold = _pupil_map(hst.WFC3IRColdMask)
    assert np.isclose((ota * cold).sum(), 39552, rtol=1.5e-3)


def test_tinytim_annular_zernikes_orthonormal():
    """TinyTim's annular Zernikes (terms 1-22) are orthonormal over the annulus 0.33 < rho < 1"""
    x = np.linspace(-1, 1, 801)
    v3, v2 = np.meshgrid(x, x)
    rho = np.hypot(v2, v3)
    annulus = (rho > 0.33) & (rho < 1)
    basis = np.array([term[annulus] for term in hst.tinytim_annular_zernike_basis(22, rho=rho, theta=np.arctan2(v3, v2))])
    gram = basis @ basis.T / annulus.sum()
    assert np.allclose(gram, np.eye(22), atol=0.02)


def test_hst_mirror_maps_and_field_zernikes():
    """With the data present, reproduce TinyTim's WFC3 IR wavefront (checked against its dumped wave map)"""
    try:
        wfc3 = stpsf.WFC3()
    except IOError:
        pytest.skip('HST data not available')
    wfc3.filter = 'F160W'
    wfc3.detector_position = (507, 507)
    primary, secondary = wfc3._get_mirror_map_planes()
    pupil = hst.WFC3IRColdMask()
    assert np.isclose(wfc3._rms_opd_nm(pupil, primary, secondary), 21.3, atol=1.5)  # maps alone; ~19 nm after mean removal
    assert np.isclose(wfc3._rms_opd_nm(pupil, primary, secondary, wfc3._get_aberrations()), 34.8, atol=1.5)
    # TinyTim drops the field-dependent focus unless a despace is set
    assert wfc3._get_field_zernikes()[4] == 0
    wfc3.despace_microns = 10
    assert np.isclose(wfc3._get_field_zernikes()[4], 0.11 + 0.0075, atol=2e-3)

    wfc3.include_si_wfe = False
    assert wfc3._get_aberrations() is None
    assert 'WFC3 aberrations' not in [plane.name for plane in wfc3.get_optical_system(fov_pixels=16).planes]


def test_wfc3_default_weights():
    """Default weights use nlambda bins across the filter; a filter tabulated at one wavelength is monochromatic"""
    try:
        wfc3 = stpsf.WFC3()
    except IOError:
        pytest.skip('HST data not available')
    wfc3.filter = 'F160W'
    wavelengths, weights = wfc3._get_weights()
    assert len(wavelengths) == 21
    assert np.isclose(weights.sum(), 1)
    wfc3.filter = 'F164N'
    wavelengths, weights = wfc3._get_weights()
    assert np.allclose(wavelengths * 1e6, 1.64336, atol=1e-4)
    assert np.allclose(weights, 1)


def test_tinytim_weights():
    """TinyTim weights are throughput * photlam * dlambda, with full-spacing dlambda at the ends, then smart-skip"""
    wave = np.array([1000.0, 1100.0, 1300.0, 1400.0])
    wavelengths, weights = hst._tinytim_weights(wave, np.ones(4), np.ones(4))
    assert np.allclose(wavelengths, wave * 1e-10)
    assert np.allclose(weights, np.array([100, 150, 150, 100]) / 500)
    # a wavelength with less than 1% of the peak weight is dropped and the rest renormalized
    wavelengths, weights = hst._tinytim_weights(wave, np.array([1, 1, 1, 0.005]), np.ones(4))
    assert len(wavelengths) == 3
    assert np.isclose(weights.sum(), 1)


def test_wfc3_ir_distortion_and_ipc_options():
    """Distortion adds *DIST extensions, and IPC (IR only) spreads flux to neighbors while conserving it"""
    try:
        wfc3 = stpsf.WFC3()
    except IOError:
        pytest.skip('HST data not available')
    wfc3.filter = 'F160W'
    mono = dict(monochromatic=1.6e-6, fov_pixels=21, oversample=4)
    wfc3.options['add_distortion'] = False
    plain = wfc3.calc_psf(**mono)
    assert [ext.name for ext in plain] == ['OVERSAMP', 'DET_SAMP']
    wfc3.options['add_distortion'] = True
    both = wfc3.calc_psf(**mono)
    wfc3.options['add_ipc'] = False
    no_ipc = wfc3.calc_psf(**mono)
    assert [ext.name for ext in both] == ['OVERSAMP', 'DET_SAMP', 'OVERDIST', 'DET_DIST']
    assert 'IPCINST' in both['DET_DIST'].header
    assert 'IPCINST' not in no_ipc['DET_DIST'].header
    assert both['DET_DIST'].data.max() < no_ipc['DET_DIST'].data.max()
    assert np.isclose(both['DET_DIST'].data.sum(), no_ipc['DET_DIST'].data.sum(), rtol=2e-3)


def test_wfc3_uvis_charge_diffusion():
    """The UVIS kernel is normalized, follows the wavelength, and conserves flux; it can be turned off"""
    try:
        wfc3 = stpsf.WFC3()
    except IOError:
        pytest.skip('HST data not available')
    wfc3.detector = 'UVIS2'
    wfc3.filter = 'F606W'
    mono = dict(fov_pixels=21, oversample=4)
    psf_400 = wfc3.calc_psf(monochromatic=0.4e-6, **mono)
    kernel_400 = wfc3._uvis_diffusion_kernel(psf_400[0].header)
    psf_800 = wfc3.calc_psf(monochromatic=0.8e-6, **mono)
    kernel_800 = wfc3._uvis_diffusion_kernel(psf_800[0].header)
    assert np.isclose(kernel_400.sum(), 1) and np.isclose(kernel_800.sum(), 1)
    assert kernel_400[2, 2] < kernel_800[2, 2]  # less diffusion at long wavelengths
    wfc3.options['add_charge_diffusion'] = False
    no_diffusion = wfc3.calc_psf(monochromatic=0.4e-6, **mono)
    assert 'DETKERN' in psf_400['DET_DIST'].header
    assert 'DETKERN' not in no_diffusion['DET_DIST'].header
    assert psf_400['DET_DIST'].data.max() < no_diffusion['DET_DIST'].data.max()
    # some flux leaves the small field
    assert np.isclose(psf_400['DET_DIST'].data.sum(), no_diffusion['DET_DIST'].data.sum(), rtol=2e-2)
