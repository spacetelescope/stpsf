"""Tests for the HST / WFC3 models that do not require the STPSF data package."""

import os

import numpy as np
import poppy
import pytest
from astropy.io import fits

from .. import hst
import stpsf


def _pupil_map(cls=hst.HSTPupil, **kwargs):
    # 512 pixel grid with the pupil radius spanning 128 pixels, as in TinyTim
    wave = poppy.Wavefront(diam=4.8, npix=512)
    wave.pixelscale = 1.2 / 128
    return cls(**kwargs).get_transmission(wave)


def test_hst_pupil_geometry():
    """The rotated OTA pupil has pads at 346 and 104 degrees, and spiders on the diagonals (from TinyTim)"""
    pupil = _pupil_map(subsample=2)
    center = 255.5
    for angle_deg in (346.4, 104.0):
        x = center + 0.886 * 128 * np.cos(np.deg2rad(angle_deg))
        y = center + 0.886 * 128 * np.sin(np.deg2rad(angle_deg))
        assert pupil[int(round(y)), int(round(x))] == 0, f'expected a pad at {angle_deg} degrees'
    # a point on the +x axis at mid radius is clear, but a point on a diagonal is in a spider
    assert pupil[256, int(center + 0.6 * 128)] == 1
    assert pupil[int(center + 0.6 * 128 * 0.7071), int(center + 0.6 * 128 * 0.7071)] == 0


def test_hst_pupil_area():
    """The OTA pupil times the IR cold mask matches the TinyTim pupil dump area (39552 pixels at R=128)"""
    ota = _pupil_map()
    cold = _pupil_map(hst.WFC3IRColdMask)
    assert np.isclose((ota * cold).sum(), 39552, rtol=3e-3)
    assert np.isclose(ota.sum() / (np.pi * 128**2), 1 - 0.33**2, atol=0.03)


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


@pytest.mark.parametrize(('x', 'y', 'label'), [(507, 507, 'ctr'), (50, 50, 'll'), (960, 960, 'ur')])
def test_wfc3_ir_monochromatic_vs_tinytim(x, y, label):
    """A monochromatic WFC3 IR PSF matches TinyTim's tiny2 output (needs $TINYTIM_REF, from make_reference_psfs.sh)"""
    ref_dir = os.environ.get('TINYTIM_REF')
    if not ref_dir:
        pytest.skip('TinyTim reference PSFs not available ($TINYTIM_REF)')
    with fits.open(os.path.join(ref_dir, f'ir_mono1600_{label}00_psf.fits')) as hdul:
        tt = hdul[0].data.astype(float)
        pixscale = hdul[0].header['PIXSCALE']
    n = tt.shape[0]
    try:
        wfc3 = stpsf.WFC3()
    except IOError:
        pytest.skip('HST data not available')
    wfc3.filter = 'F160W'
    wfc3.detector_position = (x, y)
    wfc3.pixelscale = pixscale
    # An odd field of view centers the PSF on a pixel, as in TinyTim. TinyTim normalizes to the sum of its own field.
    psf = wfc3.calc_psf(monochromatic=1.6e-6, fov_pixels=n + 1, oversample=8)['DET_SAMP'].data[:n, :n]
    diff = psf / psf.sum() - tt / tt.sum()
    peak = tt.max() / tt.sum()
    assert np.abs(diff).max() / peak < 0.005
    assert np.sqrt(np.mean(diff**2)) / peak < 5e-4
