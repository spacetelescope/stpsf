"""Tests for the HST / WFC3 models that do not require the STPSF data package."""

import os

import numpy as np
import poppy
import pytest
import synphot
from astropy import units as u
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


def _tinytim_inputs():
    """Reference PSF directory, and the G8V spectrum TinyTim used, as a synphot source; skips if unavailable"""
    ref_dir, tinytim_dir = os.environ.get('TINYTIM_REF'), os.environ.get('TINYTIM')
    if not (ref_dir and tinytim_dir):
        pytest.skip('TinyTim reference PSFs and install not available ($TINYTIM_REF, $TINYTIM)')
    spec = np.loadtxt(os.path.join(tinytim_dir, 'bpgs_53.dat'), skiprows=7)
    g8v = synphot.SourceSpectrum(
        synphot.models.Empirical1D, points=spec[:, 0] * u.AA, lookup_table=spec[:, 1] * synphot.units.PHOTLAM
    )
    return ref_dir, g8v


def _par_value(par_file, comment):
    with open(par_file) as f:
        return next(line.split()[0] for line in f if comment in line)


@pytest.mark.parametrize('filt', ['F105W', 'F125W', 'F160W'])
def test_wfc3_ir_weights_vs_tinytim(filt):
    """The wavelengths and weights match those TinyTim wrote to its .par file (which has 6 digits)"""
    ref_dir, g8v = _tinytim_inputs()
    with open(os.path.join(ref_dir, f'ir_{filt.lower()}_ctr.par')) as f:
        expected = np.array([line.split()[:2] for line in f if '# Wavelength' in line], dtype=float)
    try:
        wfc3 = stpsf.WFC3()
    except IOError:
        pytest.skip('HST data not available')
    wfc3.filter = filt
    wfc3.tinytim_wavelengths = True
    wavelengths, weights = wfc3._get_weights(source=g8v)
    assert np.allclose(wavelengths * 1e6, expected[:, 0], atol=1e-5)
    assert np.allclose(weights, expected[:, 1], atol=2e-6)


@pytest.mark.parametrize(
    ('filt', 'x', 'y', 'label'), [('F160W', 507, 507, 'ctr'), ('F160W', 50, 50, 'll'), ('F105W', 960, 960, 'ur')]
)
def test_wfc3_ir_polychromatic_vs_tinytim(filt, x, y, label):
    """A polychromatic WFC3 IR PSF with TinyTim's wavelengths matches its tiny2 output for a G8V star"""
    ref_dir, g8v = _tinytim_inputs()
    root = os.path.join(ref_dir, f'ir_{filt.lower()}_{label}')
    tt = fits.getdata(root + '00_psf.fits').astype(float)
    n = tt.shape[0]
    try:
        wfc3 = stpsf.WFC3()
    except IOError:
        pytest.skip('HST data not available')
    wfc3.filter = filt
    wfc3.detector_position = (x, y)
    wfc3.pixelscale = float(_par_value(root + '.par', 'Integrated PSF scaling'))  # the FITS header is rounded
    wfc3.tinytim_wavelengths = True
    psf = wfc3.calc_psf(source=g8v, fov_pixels=n + 1, oversample=8)['DET_SAMP'].data[:n, :n]
    diff = psf / psf.sum() - tt / tt.sum()
    peak = tt.max() / tt.sum()
    assert np.abs(diff).max() / peak < 0.005
    assert np.sqrt(np.mean(diff**2)) / peak < 5e-4
