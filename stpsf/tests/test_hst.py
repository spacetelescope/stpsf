"""Tests for the HST / WFC3 models that do not require the STPSF data package."""

import os

import matplotlib.colors
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


def _centroid(img):
    y, x = np.indices(img.shape)
    return (img * x).sum() / img.sum(), (img * y).sum() / img.sum()


def _plot_vs_tinytim(tt, psf, title):
    """Show the TinyTim PSF, the STPSF PSF, their difference and the radial profiles, all normalized to unit sum.

    The first two panels share a log color scale; the difference is in percent of the TinyTim peak.
    """
    import matplotlib.pyplot as plt

    tt_n, psf_n = tt / tt.sum(), psf / psf.sum()
    peak = tt_n.max()
    diff_pct = 100 * (psf_n - tt_n) / peak
    cen_tt, cen_psf = _centroid(tt_n), _centroid(psf_n)

    fig, axes = plt.subplots(1, 4, figsize=(18, 4.5))
    norm = matplotlib.colors.LogNorm(vmin=peak * 1e-5, vmax=peak)
    for ax, img, name in ((axes[0], tt_n, 'TinyTim'), (axes[1], psf_n, 'STPSF')):
        im = ax.imshow(np.clip(img, peak * 1e-5, None), norm=norm, origin='lower')
        ax.plot(*_centroid(img), 'r+', ms=10)
        ax.set_title(f'{name} (log)')
        fig.colorbar(im, ax=ax, fraction=0.046)
    lim = np.abs(diff_pct).max()
    im = axes[2].imshow(diff_pct, origin='lower', cmap='RdBu_r', vmin=-lim, vmax=lim)
    axes[2].set_title(f'STPSF - TinyTim, % of peak\nrms {np.sqrt(np.mean(diff_pct**2)):.3f}%, max {lim:.2f}%')
    fig.colorbar(im, ax=axes[2], fraction=0.046)

    # Fourth panel: azimuthally averaged profiles about the TinyTim centroid, and the centroid offset
    y, x = np.indices(tt.shape)
    r = np.hypot(x - cen_tt[0], y - cen_tt[1])
    bins = np.arange(0, r.max(), 1.0)
    idx = np.digitize(r.ravel(), bins)
    for img, name in ((tt_n, 'TinyTim'), (psf_n, 'STPSF')):
        prof = [img.ravel()[idx == i].mean() for i in range(1, len(bins))]
        axes[3].semilogy(bins[:-1] + 0.5, prof, label=name)
    axes[3].set_xlabel('radius (pixels)')
    axes[3].set_title(f'Radial profile\ncentroid offset ({cen_psf[0] - cen_tt[0]:+.4f}, {cen_psf[1] - cen_tt[1]:+.4f}) px')
    axes[3].legend()
    fig.suptitle(title)
    fig.tight_layout()
    plt.show()


@pytest.mark.parametrize(('x', 'y', 'label'), [(507, 507, 'ctr'), (50, 50, 'll'), (960, 960, 'ur')])
def test_wfc3_ir_monochromatic_vs_tinytim(x, y, label, plot=False):
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
    if plot:
        _plot_vs_tinytim(tt, psf, f'WFC3/IR 1.6 micron monochromatic, {label} {x, y}: tiny2')
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
def test_wfc3_ir_polychromatic_vs_tinytim(filt, x, y, label, plot=False):
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
    if plot:
        _plot_vs_tinytim(tt, psf, f'WFC3/IR {filt} G8V, {label} {x, y}: tiny2')
    diff = psf / psf.sum() - tt / tt.sum()
    peak = tt.max() / tt.sum()
    assert np.abs(diff).max() / peak < 0.005
    assert np.sqrt(np.mean(diff**2)) / peak < 5e-4


@pytest.mark.parametrize(('x', 'y', 'label'), [(507, 507, 'ctr'), (50, 50, 'll'), (960, 960, 'ur')])
def test_wfc3_ir_detector_vs_tinytim(x, y, label, plot=False):
    """DET_DIST (SIAF distortion and IPC) matches the tiny3 detector image, to the accuracy that the SIAF
    and TinyTim distortion models agree: centroids within 0.02 pixel, residuals under 0.2% rms and 3% max of peak"""
    ref_dir, g8v = _tinytim_inputs()
    tt = fits.getdata(os.path.join(ref_dir, f'ir_f160w_{label}00.fits')).astype(float)
    n = tt.shape[0]
    try:
        wfc3 = stpsf.WFC3()
    except IOError:
        pytest.skip('HST data not available')
    wfc3.filter = 'F160W'
    wfc3.detector_position = (x, y)
    wfc3.tinytim_wavelengths = True
    psf = wfc3.calc_psf(source=g8v, fov_pixels=n, oversample=8)['DET_DIST'].data
    assert psf.shape == tt.shape
    if plot:
        _plot_vs_tinytim(tt, psf, f'WFC3/IR F160W G8V, {label} {x, y}: tiny3 (SIAF distortion + IPC)')
    assert np.allclose(_centroid(psf), _centroid(tt), atol=0.02)
    diff = psf / psf.sum() - tt / tt.sum()
    peak = tt.max() / tt.sum()
    assert np.abs(diff).max() / peak < 0.03
    assert np.sqrt(np.mean(diff**2)) / peak < 2e-3


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
