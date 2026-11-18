"""Regression tests of the HST WFC3 models against TinyTim, for development; not part of the CI test suite.

These compare STPSF's WFC3 PSFs, pupils, wavefronts and spectral weights with the output of TinyTim itself (see
``stpsf/hst.py`` for the results and the differences that remain), so they need a TinyTim installation and a set of
reference PSFs that are not distributed with STPSF. The tests that need no TinyTim files are in
``stpsf/tests/test_hst.py``.

Setup:

1. ``$TINYTIM`` points at the TinyTim installation (with ``tiny2`` compiled and ``bpgs_53.dat``).
2. ``$TINYTIM_REF`` points at a directory of reference PSFs, made with::

       dev_utils/hst_make_tinytim_reference_psfs.sh -S -o /tmp/tinytim_ref

3. ``$STPSF_PATH`` points at the STPSF data, including the HST files made by ``make_hst_data.py``.

Run with pytest, naming the file explicitly (it is outside the test paths)::

    pytest dev_utils/run_hst_regression_tests.py -p no:warnings

The tests skip if the TinyTim inputs or the HST data are missing. Each test comparing PSFs takes ``plot=True`` to
show TinyTim, STPSF and their difference, for use from a notebook::

    sys.path.insert(0, '<stpsf>/dev_utils')
    import run_hst_regression_tests as regression
    regression.test_wfc3_ir_monochromatic_vs_tinytim(507, 507, 'ctr', plot=True)
"""

import csv
import os
import re
import subprocess

import matplotlib.colors
import numpy as np
import poppy
import pytest
import synphot
from astropy import units as u
from astropy.io import fits

import stpsf


def _centroid(img):
    y, x = np.indices(img.shape)
    return (img * x).sum() / img.sum(), (img * y).sum() / img.sum()


def _assert_psfs_match(psf, tt, max_tol, rms_tol):
    """Compare unit-sum PSFs; the max and rms differences, as fractions of the TinyTim peak, must be under the tolerances"""
    diff = psf / psf.sum() - tt / tt.sum()
    peak = tt.max() / tt.sum()
    max_diff, rms_diff = np.abs(diff).max() / peak, np.sqrt(np.mean(diff**2)) / peak
    print(f'PSF difference vs TinyTim: max {max_diff:.2e} (tol {max_tol:.0e}), rms {rms_diff:.2e} (tol {rms_tol:.0e})')
    assert max_diff < max_tol
    assert rms_diff < rms_tol


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
    for img, name, style in ((tt_n, 'TinyTim', '-'), (psf_n, 'STPSF', '--')):
        prof = [img.ravel()[idx == i].mean() for i in range(1, len(bins))]
        axes[3].semilogy(bins[:-1] + 0.5, prof, style, label=name)
    axes[3].set_xlabel('radius (pixels)')
    axes[3].set_title(f'Radial profile\ncentroid offset ({cen_psf[0] - cen_tt[0]:+.4f}, {cen_psf[1] - cen_tt[1]:+.4f}) px')
    axes[3].legend()
    fig.suptitle(title)
    fig.tight_layout()
    plt.show()


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


def _wfc3(det, x, y, filt):
    """A WFC3 set up for the given detector ('IR', 'UVIS1' or 'UVIS2'), position and filter"""
    try:
        wfc3 = stpsf.WFC3()
    except IOError:
        pytest.skip('HST data not available')
    if det != 'IR':
        wfc3.detector = det
    wfc3.filter = filt
    wfc3.detector_position = (x, y)
    return wfc3


_UVIS_CASES = [
    ('UVIS2', 2048, 1024, 'ctr'), ('UVIS2', 100, 100, 'll'), ('UVIS1', 2048, 1024, 'ctr'), ('UVIS1', 4000, 1950, 'ur')
]


@pytest.mark.parametrize(('x', 'y', 'label'), [(507, 507, 'ctr'), (50, 50, 'll'), (960, 960, 'ur')])
def test_wfc3_ir_monochromatic_vs_tinytim(x, y, label, plot=False):
    """A monochromatic WFC3 IR PSF matches TinyTim's tiny2 output (needs $TINYTIM_REF, from hst_make_tinytim_reference_psfs.sh)"""
    ref_dir = os.environ.get('TINYTIM_REF')
    if not ref_dir:
        pytest.skip('TinyTim reference PSFs not available ($TINYTIM_REF)')
    with fits.open(os.path.join(ref_dir, f'ir_mono1600_{label}00_psf.fits')) as hdul:
        tt = hdul[0].data.astype(float)
        pixscale = hdul[0].header['PIXSCALE']
    n = tt.shape[0]
    wfc3 = _wfc3('IR', x, y, 'F160W')
    wfc3.pixelscale = pixscale
    # An odd field of view centers the PSF on a pixel, as in TinyTim. TinyTim normalizes to the sum of its own field.
    psf = wfc3.calc_psf(monochromatic=1.6e-6, fov_pixels=n + 1, oversample=8)['DET_SAMP'].data[:n, :n]
    if plot:
        _plot_vs_tinytim(tt, psf, f'WFC3/IR 1.6 micron monochromatic, {label} {x, y}: tiny2')
    _assert_psfs_match(psf, tt, 1e-3, 1e-4)


def test_wfc3_field_zernikes_vs_tinytim_manifest():
    """The field-dependent Z4-Z8 (waves) at every position in the TinyTim reference set match those tiny1 reported
    (4 decimals); Z4 is applied only with ``tinytim_skip_zero_terms = False``, see ``_get_field_zernikes``"""
    ref_dir, _ = _tinytim_inputs()
    with open(os.path.join(ref_dir, 'MANIFEST.csv')) as f:
        rows = list(csv.DictReader(f))
    try:
        wfc3 = stpsf.WFC3()
    except IOError:
        pytest.skip('HST data not available')
    for row in rows:
        if row['camera'] != 'ir':
            wfc3.detector = row['camera'].upper()
        else:
            wfc3.filter = 'F160W'
        wfc3.detector_position = (int(row['x']), int(row['y']))
        expected = [float(row[key]) for key in ('Z4_focus', 'Z5_xastig', 'Z6_yastig', 'Z7_xcoma', 'Z8_ycoma')]
        wfc3.tinytim_skip_zero_terms = True
        assert wfc3._get_field_zernikes().get(4, 0) == 0, row['root']
        wfc3.tinytim_skip_zero_terms = False
        zernikes = wfc3._get_field_zernikes()
        assert np.allclose([zernikes[j] for j in range(4, 9)], expected, atol=1e-4), row['root']


@pytest.mark.parametrize(
    ('root', 'det', 'filt', 'x', 'y'),
    [('ir_mono1600_ctr', 'IR', 'F160W', 507, 507), ('ir_mono1600_ur', 'IR', 'F160W', 960, 960),
     ('uvis2_mono600_ctr', 'UVIS2', 'F606W', 2048, 1024), ('uvis1_mono600_ur', 'UVIS1', 'F606W', 4000, 1950)],
)
def test_wfc3_pupil_and_opd_vs_tinytim(root, det, filt, x, y, tmp_path):
    """The pupil transmission and total wavefront (OPD) match the maps TinyTim's ``tiny2`` writes out, on its grid.

    The pupil agrees in open area to 0.1% and in mean absolute transmission to 0.2%, the differences being at the
    edges (TinyTim's pupil is drawn with 8x8 subsampling and then bilinearly rotated); the wavefronts, which also
    check the orientation of the field terms and mirror maps, agree to 0.05 nm rms after removing piston.
    """
    ref_dir, _ = _tinytim_inputs()
    tiny2 = os.path.join(os.environ['TINYTIM'], 'tiny2')
    if not os.path.exists(tiny2):
        pytest.skip('TinyTim executables not available ($TINYTIM/tiny2)')
    with open(os.path.join(ref_dir, f'{root}.par')) as f:
        par = f.read()
    par = re.sub(r'(?m)^0(\s+# Write out)', r'1\1', par)  # also write the pupil, wave and critical PSF maps
    par = re.sub(r'(?m)^\S+(\s+# Output PSF file rootname)', r'wf\1', par)
    (tmp_path / 'wf.par').write_text(par)
    subprocess.run([tiny2, 'wf.par'], cwd=tmp_path, check=True, capture_output=True)
    wavelength_um, _, dim = next(line.split()[:3] for line in par.splitlines() if '# Wavelength 1' in line)
    dim = int(dim)
    tt_pupil = fits.getdata(tmp_path / 'tinytim_pupil.fits').astype(float)
    tt_opd = fits.getdata(tmp_path / 'tinytim_wave.fits').astype(float) * float(wavelength_um) * 1e3  # nm

    wfc3 = _wfc3(det, x, y, filt)
    # TinyTim's pupil has radius dim/4 pixels with the optical axis at pixel dim/2
    dx = wfc3.pupil_radius.to_value(u.m) / (dim / 4)
    wave = poppy.Wavefront(diam=(dim + 1) * dx, npix=dim + 1)
    pupil, opd = np.ones(wave.shape), np.zeros(wave.shape)
    for plane in wfc3.get_optical_system(fov_pixels=16).planes:
        if not isinstance(plane, poppy.Detector):
            pupil *= plane.get_transmission(wave)
            opd += plane.get_opd(wave)
    pupil, opd = pupil[:dim, :dim], opd[:dim, :dim] * 1e9

    assert np.isclose(pupil.sum(), tt_pupil.sum(), rtol=1.5e-3)
    assert np.abs(pupil - tt_pupil).mean() < 1.5e-3
    assert np.abs(pupil - tt_pupil).max() < 0.3
    both = (pupil > 0.5) & (tt_pupil > 0.5)
    diff = (opd - opd[both].mean() - (tt_opd - tt_opd[both].mean()))[both]
    assert np.sqrt(np.mean(diff**2)) < 0.05
    assert np.abs(diff).max() < 0.3
    assert np.isclose(np.std(opd[both]), np.std(tt_opd[both]), atol=0.05)


@pytest.mark.parametrize('filt', ['F105W', 'F125W', 'F160W'])
def test_wfc3_ir_weights_vs_tinytim(filt):
    """The wavelengths and weights match those TinyTim wrote to its .par file (which has 6 digits)"""
    ref_dir, g8v = _tinytim_inputs()
    with open(os.path.join(ref_dir, f'ir_{filt.lower()}_ctr.par')) as f:
        expected = np.array([line.split()[:2] for line in f if '# Wavelength' in line], dtype=float)
    wfc3 = _wfc3('IR', 507, 507, filt)
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
    wfc3 = _wfc3('IR', x, y, filt)
    wfc3.pixelscale = float(_par_value(root + '.par', 'Integrated PSF scaling'))  # the FITS header is rounded
    wfc3.tinytim_wavelengths = True
    psf = wfc3.calc_psf(source=g8v, fov_pixels=n + 1, oversample=8)['DET_SAMP'].data[:n, :n]
    if plot:
        _plot_vs_tinytim(tt, psf, f'WFC3/IR {filt} G8V, {label} {x, y}: tiny2')
    _assert_psfs_match(psf, tt, 2.5e-3, 1e-4)


@pytest.mark.parametrize(('x', 'y', 'label'), [(507, 507, 'ctr'), (50, 50, 'll'), (960, 960, 'ur')])
def test_wfc3_ir_detector_vs_tinytim(x, y, label, plot=False):
    """DET_DIST (SIAF distortion and IPC) matches the tiny3 detector image, to the accuracy that the SIAF
    and TinyTim distortion models agree: centroids within 0.02 pixel, residuals under 0.15% rms and 3% max of peak"""
    ref_dir, g8v = _tinytim_inputs()
    tt = fits.getdata(os.path.join(ref_dir, f'ir_f160w_{label}00.fits')).astype(float)
    n = tt.shape[0]
    wfc3 = _wfc3('IR', x, y, 'F160W')
    wfc3.tinytim_wavelengths = True
    psf = wfc3.calc_psf(source=g8v, fov_pixels=n, oversample=8)['DET_DIST'].data
    assert psf.shape == tt.shape
    if plot:
        _plot_vs_tinytim(tt, psf, f'WFC3/IR F160W G8V, {label} {x, y}: tiny3 (SIAF distortion + IPC)')
    assert np.allclose(_centroid(psf), _centroid(tt), atol=0.02)
    _assert_psfs_match(psf, tt, 3e-2, 1.5e-3)


@pytest.mark.parametrize(('det', 'x', 'y', 'label'), _UVIS_CASES)
def test_wfc3_uvis_monochromatic_vs_tinytim(det, x, y, label, plot=False):
    """A monochromatic WFC3 UVIS PSF matches TinyTim's tiny2 output (needs $TINYTIM_REF)"""
    ref_dir = os.environ.get('TINYTIM_REF')
    if not ref_dir:
        pytest.skip('TinyTim reference PSFs not available ($TINYTIM_REF)')
    root = os.path.join(ref_dir, f'{det.lower()}_mono600_{label}')
    tt = fits.getdata(root + '00_psf.fits').astype(float)
    n = tt.shape[0]
    wfc3 = _wfc3(det, x, y, 'F606W')
    wfc3.pixelscale = float(_par_value(root + '.par', 'Integrated PSF scaling'))
    psf = wfc3.calc_psf(monochromatic=0.6e-6, fov_pixels=n + 1, oversample=8)['DET_SAMP'].data[:n, :n]
    if plot:
        _plot_vs_tinytim(tt, psf, f'WFC3/{det} 0.6 micron monochromatic, {label} {x, y}: tiny2')
    _assert_psfs_match(psf, tt, 1e-3, 1e-4)


@pytest.mark.parametrize(
    ('filt', 'det', 'x', 'y', 'label'),
    [('F275W', 'UVIS2', 2048, 1024, 'ctr'), ('F438W', 'UVIS1', 2048, 1024, 'ctr'), ('F606W', 'UVIS1', 4000, 1950, 'ur'),
     ('F814W', 'UVIS2', 2048, 1024, 'ctr')],
)
def test_wfc3_uvis_polychromatic_vs_tinytim(filt, det, x, y, label, plot=False):
    """Polychromatic UVIS PSFs with TinyTim's wavelengths match tiny2 (tiny3, with distortion and charge diffusion,
    to the level that the SIAF and TinyTim distortion models agree, as the PSF is undersampled)"""
    ref_dir, g8v = _tinytim_inputs()
    root = os.path.join(ref_dir, f'{det.lower()}_{filt.lower()}_{label}')
    tt = fits.getdata(root + '00_psf.fits').astype(float)
    t3 = fits.getdata(root + '00.fits').astype(float)
    n = tt.shape[0]
    wfc3 = _wfc3(det, x, y, filt)
    wfc3.tinytim_wavelengths = True
    wfc3.pixelscale = float(_par_value(root + '.par', 'Integrated PSF scaling'))
    psf = wfc3.calc_psf(source=g8v, fov_pixels=n + 1, oversample=8)['DET_SAMP'].data[:n, :n]
    if plot:
        _plot_vs_tinytim(tt, psf, f'WFC3/{det} {filt} G8V, {label} {x, y}: tiny2')
    _assert_psfs_match(psf, tt, 2.5e-3, 5e-5)

    wfc3 = _wfc3(det, x, y, filt)  # the detector image uses the SIAF pixel scale
    wfc3.tinytim_wavelengths = True
    psf = wfc3.calc_psf(source=g8v, fov_pixels=t3.shape[0], oversample=8)['DET_DIST'].data
    assert psf.shape == t3.shape
    if plot:
        _plot_vs_tinytim(t3, psf, f'WFC3/{det} {filt} G8V, {label} {x, y}: tiny3 (SIAF distortion + charge diffusion)')
    _assert_psfs_match(psf, t3, 5e-2, 1e-3)
