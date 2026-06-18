import pytest
import numpy as np
from astropy.io import fits

from .. import distortion
from .. import stpsf_core


def test_apply_distortion_skew():
    """
    Test the PSF axis is skewed appropriately by the apply_distortion function.

    Check that a rectangle of 1s run through the FGS1 distortion will be skewed (since FGS1 is heavily skewed in a
    known way). We'll check that the top left corner of the rectangle is higher up than the top right corner of the
    rectangle by checking the indices where the 1s begin/end.

    """

    # Create a baseline PSF to have shape/header keywords correct
    fgs = stpsf_core.FGS()
    fgs.detector = 'FGS1'
    fgs.options['output_mode'] = 'Oversampled image'
    psf = fgs.calc_psf(add_distortion=False)

    # Set up new extensions (from stpsf_core.JWInstrument._calc_psf_format_output)
    n_exts = len(psf)
    for ext in np.arange(n_exts):
        hdu_new = fits.ImageHDU(psf[ext].data, psf[ext].header)  # these will be the PSFs that are edited
        psf.append(hdu_new)
        ext_new = ext + n_exts
        psf[ext_new].header['EXTNAME'] = psf[ext].header['EXTNAME'][0:4] + 'DIST'  # change extension name

    # Run data through the distortion function
    psf_siaf = distortion.apply_distortion(psf)

    # Rebin data to get 3rd extension
    fgs.options['output_mode'] = 'Both extensions'
    fgs.options['detector_oversample'] = psf[0].header['DET_SAMP']
    stpsf_core.SpaceTelescopeInstrument._calc_psf_format_output(fgs, result=psf_siaf, options=fgs.options)

    # Test the slope of the rectangle
    for ext in [2, 3]:
        left = psf_siaf[ext].data[:, 0]  # isolate the far left column
        right = psf_siaf[ext].data[:, -1]  # isolate the far right column

        indexes_left = [i for i, x in enumerate(left) if x != 0.0]  # find the indices of the rectangle in the left col
        indexes_right = [i for i, x in enumerate(right) if x != 0.0]  # find the indices of the rectangle in the right

        top_of_left = np.min(indexes_left)  # find the index of the top left corner of the rectangle
        top_of_right = np.min(indexes_right)  # find the index of the top right corner of the rectangle

        # Assert that the top of left > top of right due to FGS1 skew
        assert top_of_left > top_of_right, 'FGS PSF does not have expected skew after distortion application'


def test_apply_distortion_pixel_scale():
    """
    Test the pixel scale is changed by the apply_distortion function.

    Create a fake data set that has rows of constant value, so row 0 is all 0s, row 1 is all 1s, etc. Then distort it
    via apply_distortion(), which will change both the pixel scale and skew the data. If there was no skew in the
    data, there'd be a constant x and y pixel scale change, which wouldn't really affect the x direction since the
    function would be blending the same values together, and it would affect the y direction, but by the same amount
    (since again, the function is blending the same values at each index (ie blending 0 and 1 across the entire row).

    So subtract the linear shape caused by the skew out of the row and then check that across a specific row,
    the newly pixel-scale distorted ("blended") values are approximately equal.

    Use FGS1 for its large pixel scale change in this test

    """

    # Create a baseline PSF to have shape/header keywords correct
    fgs = stpsf_core.FGS()
    fgs.detector = 'FGS1'
    fgs.options['output_mode'] = 'Oversampled image'
    psf = fgs.calc_psf(add_distortion=False)

    # Replace data with a fake image of row values equal to the row number
    data = np.zeros_like(psf[0].data)
    ny, nx = data.shape
    for i in np.arange(ny):
        data[i, :] = i

    # Replace the data in the PSF with the fake image
    psf[0].data = data

    # Set up new extensions (from stpsf_core.JWInstrument._calc_psf_format_output)
    n_exts = len(psf)
    for ext in np.arange(n_exts):
        hdu_new = fits.ImageHDU(psf[ext].data, psf[ext].header)  # these will be the PSFs that are edited
        psf.append(hdu_new)
        ext_new = ext + n_exts
        psf[ext_new].header['EXTNAME'] = psf[ext].header['EXTNAME'][0:4] + 'DIST'  # change extension name

    # Run data through the distortion function
    psf_siaf = distortion.apply_distortion(psf)

    # Rebin data to get 3rd extension (DET_DIST)
    fgs.options['output_mode'] = 'Both extensions'
    fgs.options['detector_oversample'] = psf[0].header['DET_SAMP']
    stpsf_core.SpaceTelescopeInstrument._calc_psf_format_output(fgs, result=psf_siaf, options=fgs.options)

    # Test that the change caused by the pixel distortion is approximately constant along the row
    # Choosing to check the 20th row.
    i = 20
    ext = 3

    # Crop off the edges due to skew / rotation that brings in 0s from beyond edge of detector
    psf_arr = psf_siaf[ext].data[5:-5, 5:-5]
    ncol = psf_arr.shape[1]
    inds = np.arange(ncol)

    # Model the skew with a basic linear function
    slope, intercept = np.polyfit(inds, psf_arr[i, :], 1)
    linear = (slope * inds) + intercept

    # Create a new 1D array that's your 20th row with the linear skew subtracted out
    final = psf_arr[i, :] - linear

    # Check the difference between adjacent values is the same to 1 decimal place
    diff = final[:-1] - final[1:]
    assert pytest.approx(diff, abs=0.1) == 0, 'FGS PSF does not have expected pixel scale distortion for adjacent pixels'

    # Check that the difference between the first and last value is also the same to 1 decimal
    assert pytest.approx(final[-1], abs=0.1) == final[0], (
        'FGS PSF does not have expected pixel scale distortion in the ' 'entire row'
    )


def test_apply_rotation_error():
    """Test that the apply_rotation function raises an error for NIRSpec and MIRI PSFs"""

    # Create a PSF
    for inst in [stpsf_core.NIRSpec(), stpsf_core.MIRI()]:
        psf = inst.calc_psf(nlambda=1)  # done for speed

        # Test that running this function will raise a ValueError
        with pytest.raises(ValueError) as excinfo:
            distortion.apply_rotation(psf)
        assert 'ValueError' in str(excinfo), 'NIRSpec & MIRI PSFs should not be able to run through apply_rotation'


def test_distortion_with_custom_pixscale():
    """Verifies the distortion model works properly even if the pixel scale is changed to
    a nonstandard value for the calculation. This tests/verifies the fix in PR 669:
        https://github.com/spacetelescope/stpsf/pull/669
    """

    miri = stpsf_core.MIRI()
    miri.pixelscale = 0.061
    psf = miri.calc_psf(fov_arcsec=2)

    # A symptom of the prior bug was the total sum of a distorted PSF would be very
    # discrepant from the sum of the undistorted PSF. So verif that symptom is not the case:

    assert np.isclose(psf[0].data.sum(), psf[3].data.sum(), rtol=0.001)
    assert np.isclose(psf[1].data.sum(), psf[3].data.sum(), rtol=0.001)


def test_distortion_linear_anisotropy():
    """
    Test that distort_image is correctly reproducing the expected anisotropy in pixel scale between the x and y axes

    To do this, we generate a 2D gaussian of known "ideal" dimensions and check that the correct FWHMs are recovered after distortion is applied.

    See https://github.com/spacetelescope/stpsf/issues/148
    """
    from copy import deepcopy
    from scipy.interpolate import interp1d

    def gaussian_2d(X, Y, mu, A=1.0, sigma=1.0):
        x0, y0 = mu
        r2 = (X - x0) ** 2 + (Y - y0) ** 2
        return A * np.exp(-0.5 * r2 / sigma ** 2)

    def measure_2d_fwhms(image, mu):
        xsli = image[mu[1], :]
        ysli = image[:, mu[0]]

        xyv = np.arange(len(xsli))

        x2 = interp1d(xsli[mu[0]:], xyv[mu[0]:])(0.5)
        x1 = interp1d(xsli[:mu[0]], xyv[:mu[0]])(0.5)
        xfwhm = x2 - x1

        y2 = interp1d(ysli[mu[1]:], xyv[mu[1]:])(0.5)
        y1 = interp1d(ysli[:mu[1]], xyv[:mu[1]])(0.5)
        yfwhm = y2 - y1

        return xfwhm, yfwhm

    nrc = stpsf_core.NIRCam()
    nrc.filter = 'F444W'
    nrc.pupil_mask = 'MASKRND'
    nrc.image_mask = 'MASK335R'
    nrc.aperturename = 'NRCA5_MASK335R'
    nrc.set_position_from_aperture_name(nrc.aperturename)

    hdul = nrc.calc_psf(add_distortion=False, monochromatic=4.4e-6,
                        fov_pixels=120)  # We're primarily running this to get the right header info

    ext = 0

    # Generate a broad 2D gaussian where X/Y scale can be measured more robustly
    fwhm = 100 * hdul[ext].header['OVERSAMP']
    sigma = fwhm / np.sqrt(8. * np.log(2.))
    mu = np.round(((np.array(hdul[ext].data.shape)[::-1]) - 1) / 2.).astype(
        int)  # Center on nearest pixel to geometric array center
    Y, X = np.indices(hdul[ext].data.shape, dtype=np.float32)
    G = gaussian_2d(X, Y, mu, sigma=sigma)

    # Replace PSF model image in the PSF HDUList with the gaussian
    hdul[ext].data = G

    # Create simplified SIAF aperture to be sure only anisotropy between X and Y scales are relevant.
    aper = deepcopy(nrc.siaf[nrc.aperturename])
    for key in aper.__dict__:
        if key.startswith('Idl2Sci') and aper.__dict__[key] is not None:
            keyinv = key.replace('Idl2Sci', 'Sci2Idl')
            if key not in ['Idl2SciX10', 'Idl2SciY11']:  # Zero all coeffs besides X10 and Y11
                aper.__dict__[key] = aper.__dict__[keyinv] = 0.
            else:  # Otherwise, ensure reversible transform
                aper.__dict__[keyinv] = 1. / aper.__dict__[key]

    G_distorted = distortion.distort_image(hdul, ext=ext, to_frame='sci', fill_value=0, aper=aper)

    xfwhm, yfwhm = measure_2d_fwhms(G_distorted, mu)

    xfwhm_expected = fwhm * nrc.pixelscale / aper.XSciScale
    yfwhm_expected = fwhm * nrc.pixelscale / aper.YSciScale

    assert pytest.approx(xfwhm, rel=1e-5) == xfwhm_expected, (
        'Distorted image does not have the expected X-axis pixel scale'
    )

    assert pytest.approx(yfwhm, rel=1e-5) == yfwhm_expected, (
        'Distorted image does not have the expected Y-axis pixel scale'
    )


def test_distortion_pixel_coords_precisely_match_siaf(plot=False):
    """Test that the distortion code can interpolate + resample pixels onto a grid
    that precisely matches the actual pixel Ideal frame coordinates as specified in SIAF.

    This tests proper handling of various aspects of the SIAF usage and pixel coordinate setup

    See https://github.com/spacetelescope/stpsf/issues/148
    """

    # Setup a sim instance on a subarray
    nrc = stpsf_core.NIRCam()
    nrc.filter = 'F444W'
    nrc.aperturename = 'NRCA5_MASK335R'
    aper = nrc._detector_geom_info.aperture
    aper_npix = aper.XSciSize
    nrc.detector_position = (aper_npix/2, aper_npix/2)  # enforce exact centering for this test, even if the actual reference location is elsewhere in the aperture

    # Compute a PSF, without distortion for now
    psf = nrc.calc_psf(fov_pixels = aper_npix, oversample=1, add_distortion=False, monochromatic=4e-6)

    ext = 1  # DET_SAMP
    fill_value = 0

    ### Compute distortion, including inference of pixel coordinates in Science frame
    psf_dist, xnew_idl, ynew_idl = distortion.distort_image(psf, ext, to_frame='sci', fill_value=fill_value, aper=aper, return_coords=True)

    # Get pixel indices in SIAF Science frame
    # Recall that SIAF pixel indices are 1-based. This is off by 1 from 0-based Python array indices.
    # as per Colin Cox in JWST-STScI-001550
    #  "All pixel counting for JWST starts with (1,1) being the central point within the first pixel."
    y_aper_sci, x_aper_sci = np.indices((aper_npix, aper_npix))
    x_aper_sci += 1
    y_aper_sci += 1

    # Now convert those to Ideal frame
    x_aper_idl, y_aper_idl = aper.convert(x_aper_sci, y_aper_sci, 'sci', 'idl')

    # Check the coords used in distort_image do match the ones defined in the SIAF
    assert np.allclose(xnew_idl, x_aper_idl), "X coordinates used in distort_image ought to precisely match SIAF for this setup."
    assert np.allclose(ynew_idl, y_aper_idl), "Y coordinates used in distort_image ought to precisely match SIAF for this setup."

    if plot:
        import matplotlib, matplotlib.pyplot as plt

        #--- Compute regular grid of pixels relative to the SciRef location
        cen = (aper_npix-1)/2+1
        regular_ideal_grid_x = (x_aper_sci - aper.XSciRef) * nrc.pixelscale
        regular_ideal_grid_y = (y_aper_sci - aper.YSciRef) * nrc.pixelscale

        fig, ax = plt.subplots(figsize=(12,12))
        nevery = 10
        ax.scatter(regular_ideal_grid_x[::nevery, ::nevery], regular_ideal_grid_y[::nevery, ::nevery],
                   marker='+', label='Regular grid at uniform pixelscale')
        ax.scatter(xnew_idl[::nevery, ::nevery], ynew_idl[::nevery, ::nevery],
                   marker='+', label='x/y new grid used in distort_image')
        ax.scatter(x_aper_idl[::nevery, ::nevery], y_aper_idl[::nevery, ::nevery],
                  marker='+', label='aperture pixel Ideal coords from SIAF', s=100, zorder=-10)

        cornerx, cornery = aper.corners('idl')   # Note, corners are the **outer** corners of the aperture,
                                                 # offset by half a pixel outwards relative to the pixel centers

        ax.plot(cornerx[[0,1,2,3,0]], cornery[[0,1,2,3,0]],  # use extra indices to close the square
                color='black',
                label='Aperture Border')

        ax.set_xlim(5, 11)
        ax.set_ylim(-11.5, -5.5)

        ax.legend(framealpha=1.0, loc='upper right')
        ax.set_title(f"Pixel sampling locations in the Ideal frame for {aper.AperName}\nShowing every {nevery}th pixel per axis, and zoomed in on a corner",
                    )
