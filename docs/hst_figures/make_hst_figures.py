"""Make the figures used in docs/hst.rst and docs/hst_validation_vs_tinytim.rst.

Usage, from anywhere::

    python make_hst_figures.py examples      # example PSFs; needs only STPSF and its data
    python make_hst_figures.py validation    # comparisons with Tiny Tim; needs $TINYTIM and $TINYTIM_REF
    python make_hst_figures.py               # both

The validation figures are made with the ``plot=True`` option of the tests in ``dev_utils/run_hst_regression_tests.py``,
which also describes how to set up Tiny Tim and the reference PSFs. Figures are written to this directory.
"""

import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

import stpsf

OUTDIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(OUTDIR, '..', '..', 'dev_utils'))


def make_example_figures():
    """Example PSFs for the overview page"""
    # Three PSFs: IR and UVIS polychromatic with detector effects, and a monochromatic oversampled UVIS PSF
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    wfc3 = stpsf.WFC3()
    wfc3.filter = 'F160W'
    psf = wfc3.calc_psf(fov_pixels=61, oversample=4)
    stpsf.display_psf(psf, ext='DET_DIST', ax=axes[0], title='WFC3/IR F160W', imagecrop=3.5)
    wfc3.filter = 'F606W'
    psf = wfc3.calc_psf(fov_pixels=101, oversample=4)
    stpsf.display_psf(psf, ext='DET_DIST', ax=axes[1], title='WFC3/UVIS1 F606W', imagecrop=2.0)
    psf = wfc3.calc_psf(fov_pixels=101, oversample=4, monochromatic=0.6e-6)
    stpsf.display_psf(psf, ext='OVERSAMP', ax=axes[2], title='UVIS1 0.6 micron, oversampled', imagecrop=2.0)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, 'hst_wfc3_example_psfs.png'), dpi=110, bbox_inches='tight')
    plt.close(fig)

    # All four extensions of one PSF, as in usage.ipynb
    wfc3 = stpsf.WFC3()
    wfc3.filter = 'F160W'
    psf = wfc3.calc_psf(fov_pixels=41, oversample=4)
    fig, axes = plt.subplots(figsize=(12, 3), ncols=4)
    for ext in range(len(psf)):
        stpsf.display_psf(
            psf, ext=ext, ax=axes[ext], title=f'Ext {ext}: {psf[ext].header["EXTNAME"]}', imagecrop=2, colorbar=False
        )
    fig.savefig(os.path.join(OUTDIR, 'hst_wfc3_extensions.png'), dpi=110, bbox_inches='tight')
    plt.close(fig)


def make_validation_figures():
    """Comparisons with Tiny Tim, saved from the plots made by the regression tests"""
    import run_hst_regression_tests as regression

    names = []

    def save_and_close():
        plt.gcf().savefig(os.path.join(OUTDIR, f'{names.pop(0)}.png'), dpi=110, bbox_inches='tight')
        plt.close('all')

    plt.show = save_and_close  # the test plotting function calls plt.show() once per comparison

    names[:] = ['validation_ir_f160w']
    regression.test_wfc3_ir_polychromatic_vs_tinytim('F160W', 507, 507, 'ctr', plot=True)
    names[:] = ['validation_ir_detector']
    regression.test_wfc3_ir_detector_vs_tinytim(507, 507, 'ctr', plot=True)
    # this test plots the undistorted (tiny2) and then the distorted (tiny3) comparison
    names[:] = ['validation_uvis_f606w_tiny2', 'validation_uvis_f606w_tiny3']
    regression.test_wfc3_uvis_polychromatic_vs_tinytim('F606W', 'UVIS1', 4000, 1950, 'ur', plot=True)


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if which in ('examples', 'all'):
        make_example_figures()
    if which in ('validation', 'all'):
        make_validation_figures()
