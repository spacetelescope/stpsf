import copy

import astropy.io.fits as fits
import numpy as np
import pysiaf
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import rotate

import stpsf.stpsf_core


def _get_default_siaf(instrument, aper_name):
    """
    Create instance of pysiaf for the input instrument and aperture
    to be used later to pull SIAF values like distortion polynomial
    coefficients and rotation.

    Parameters
    ----------
    instrument : str
        The name of the instrument
    aper_name : str
        The name of the specific instrument aperture

    Returns
    -------
    aper : instance of pysiaf
    """

    # Create new naming because SIAF requires special capitalization
    if instrument == 'NIRCAM':
        siaf_name = 'NIRCam'
    elif instrument == 'NIRSPEC':
        siaf_name = 'NIRSpec'
    else:
        siaf_name = instrument

    # Select a single SIAF aperture
    if instrument == 'WFI':
        siaf = pysiaf.Siaf('Roman')
        aper = siaf[aper_name]
    else:
        siaf = stpsf.stpsf_core.get_siaf_with_caching(siaf_name)
        aper = siaf.apertures[aper_name]

    return aper


# Functions for applying distortion from SIAF polynomials
def distort_image(
    hdulist_or_filename,
    ext=0,
    to_frame='sci',
    fill_value=0,
    return_coords=False,
    aper=None,
):
    """Distort an image

    Apply SIAF instrument distortion to an image that is assumed to be in
    its ideal coordinates. The header information should contain the relevant
    SIAF point information, such as SI instrument, aperture name, pixel scale,
    detector oversampling, and detector position ('sci' coords).

    This function then transforms the image to the new coordinate system using
    scipy's RegularGridInterpolator (linear interpolation).

    Parameters
    ----------
    hdulist_or_filename : str or HDUList
        A PSF from STPSF, either as an HDUlist object or as a filename
    ext : int
        Extension of HDUList to perform distortion on.
    fill_value : float or None
        Value used to fill in any blank space by the skewed PSF. Default = 0.
        If set to None, values outside the domain are extrapolated.
    to_frame : str
        Requested type of output coordinate frame. DEPRECATED, now must be 'sci'
            * 'sci': pixels, in conventional DMS axes orientation
    return_coords : bool
        In addition to returning the final image, setting this to True
        will return the full set of new coordinates. Output will then
        be (psf_new, xnew, ynew), where all three array have the same
        shape.
    aper : None or :mod:`pysiaf.Aperture`
        Option to pass the SIAF aperture if it is already known or
        specified to save time on generating a new one. If set to None,
        then automatically determines a new `pysiaf` aperture based on
        information stored in the header.
    """

    # Read in input PSF
    if isinstance(hdulist_or_filename, str):
        hdu_list = fits.open(hdulist_or_filename)
    elif isinstance(hdulist_or_filename, fits.HDUList):
        hdu_list = hdulist_or_filename
    else:
        raise ValueError('input must be a filename or HDUlist')

    if aper is None:
        # Log instrument and detector names
        instrument = hdu_list[0].header['INSTRUME'].upper().strip()

        if instrument == 'WFI':
            aper_name = hdu_list[0].header['DETECTOR'] + '_FULL'
        else:
            aper_name = hdu_list[0].header['APERNAME'].upper()

        # Pull default values
        aper = _get_default_siaf(instrument, aper_name)

    # Pixel scale information
    ny, nx = hdu_list[ext].shape
    pixelscale = hdu_list[ext].header['PIXELSCL']  # the pixel scale carries the over-sample value

    # Subtle issue: How to handle oversampling and pixel scale depends on context, and
    # in particular whether the user may have overridden the normal pixelscale with some custom
    # value. If the simulated pixelscale does not match the actual SIAF pixelscale, then the
    # oversampling relative to the SCI coordinate frame will not simply be the OVERSAMP keyword
    # value. To accomodate that edge case, compute here the oversampling as the ratio of the
    # computed pixelscale relative to the nominal average scale used in SIAF.
    osamp = 0.5*(aper.XSciScale + aper.YSciScale) / pixelscale

    # Get 'sci' reference location where PSF is observed
    xcen_sci = hdu_list[ext].header['DET_X']  # center x location in pixels ('sci')
    ycen_sci = hdu_list[ext].header['DET_Y']  # center y location in pixels ('sci')

    # Convert the PSF center point from pixels to arcseconds using pysiaf
    xcen_idl, ycen_idl = aper.sci_to_idl(xcen_sci, ycen_sci)

    # ###############################################
    # Create an array of indices (in pixels) for where the PSF is located on the detector
    nx_half, ny_half = ((nx - 1) / 2.0, (ny - 1) / 2.0)
    xlin_pix = np.linspace(-1 * nx_half, nx_half, nx)
    ylin_pix = np.linspace(-1 * ny_half, ny_half, ny)
    xarr, yarr = np.meshgrid(xlin_pix, ylin_pix)

    # ###############################################
    # Create an array of indices (in pixels) that the final data will be interpolated onto
    if to_frame == 'sci':
        xnew = xarr / osamp + xcen_sci
        ynew = yarr / osamp + ycen_sci
    else:
        raise NotImplementedError('Only transforms to SCI frame are now supported for distortions.')

    # Convert requested coordinates to 'idl' coordinates
    xnew_idl, ynew_idl = aper.convert(xnew, ynew, to_frame, 'idl')

    # ###############################################
    # Interpolate using Regular Grid Interpolator
    # First we compute the coords in Ideal frame for the input PSF sim on a regular grid
    xvals = xlin_pix * pixelscale + xcen_idl
    yvals = ylin_pix * pixelscale + ycen_idl
    func = RegularGridInterpolator(
        (yvals, xvals), hdu_list[ext].data, method='linear', bounds_error=False, fill_value=fill_value
    )

    # Create an array of (yidl, xidl) values to interpolate onto
    pts = np.array([ynew_idl.flatten(), xnew_idl.flatten()]).transpose()
    # And then we interpolate from the input sampling onto the desired output sampling
    psf_new = func(pts).reshape(xnew.shape)

    if return_coords:
        return (psf_new, xnew_idl, ynew_idl)
    else:
        return psf_new


def apply_distortion(hdulist_or_filename=None, fill_value=0):
    """
    Apply a distortion to the input PSF. The distortion comes from the SIAF 4-5 degree polynomial
    (depending on the instrument). This function pulls and applies the SIAF polynomial values
    using pysiaf package, which ensures the most up-to-date values will be called.

    Parameters
    ----------
    hdulist_or_filename :
        A PSF from STPSF, either as an HDUlist object or as a filename
    fill_value : float
        Value used to fill in any blank space by the skewed PSF. Default = 0

    Returns
    -------
    psf : HDUlist object
        PSF with distortion applied from SIAF polynomial
    """

    # Read in input PSF
    if isinstance(hdulist_or_filename, str):
        hdu_list = fits.open(hdulist_or_filename)
    elif isinstance(hdulist_or_filename, fits.HDUList):
        hdu_list = hdulist_or_filename
    else:
        raise ValueError('input must be a filename or HDUlist')

    # Create a copy of the PSF
    psf = copy.deepcopy(hdu_list)
    ext = 1  # edit the oversampled PSF (OVERDIST extension)

    # Log instrument and detector names
    instrument = hdu_list[0].header['INSTRUME'].upper().strip()
    if instrument == 'WFI':
        aper_name = hdu_list[0].header['DETECTOR'] + '_FULL'
    else:
        aper_name = hdu_list[0].header['APERNAME'].upper()

    # Pull default values
    aper = _get_default_siaf(instrument, aper_name)

    # Distort grid through interpolation
    psf_new = distort_image(psf, ext, to_frame='sci', fill_value=fill_value, aper=aper)

    # Apply data to correct extensions
    psf[ext].data = psf_new

    # Set new header keywords
    psf[ext].header['DISTORT'] = ('True', 'SIAF distortion coefficients applied')
    psf[ext].header['SIAF_VER'] = (pysiaf.JWST_PRD_VERSION, 'SIAF PRD version used')

    degree = int(getattr(aper, 'Sci2IdlDeg'))
    number_of_coefficients = int((degree + 1) * (degree + 2) / 2)
    all_keys = aper.__dict__.keys()
    for axis in ['X', 'Y']:
        coeff_keys = np.sort(np.array([c for c in all_keys if 'Idl2Sci' + axis in c]))
        coeff = np.array([getattr(aper, c) for c in coeff_keys[0:number_of_coefficients]])
        for i in range(len(coeff)):
            key = 'COEF_{}'.format(coeff_keys[i][-3:])
            psf[ext].header[key] = (coeff[i], 'SIAF distortion coefficient for {}'.format(coeff_keys[i]))

    return psf


# Function for applying Rotation to NIRCam, NIRISS, and FGS


def apply_rotation(hdulist_or_filename=None, rotate_value=None, crop=True):
    """
    Apply the detector's rotation to the PSF. This is for NIRCam, NIRISS, and FGS.
    MIRI and NIRSpec's large rotation is already added inside STPSF's calculations.

    Parameters
    ----------
    hdulist_or_filename :
        A PSF from STPSF, either as an HDUlist object or as a filename
    rotate_value : float
        Rotation in degrees that PSF needs to be. If set to None, function
        will pull the most up to date SIAF value. Default = None.
    crop : bool
        True or False to crop the PSF so it matches the size of the input
        PSF (e.g. so they could be more easily compared).

    Returns
    -------
    psf : HDUlist object
        PSF with rotation applied from SIAF values
    """

    # Read in input PSF
    if isinstance(hdulist_or_filename, str):
        hdu_list = fits.open(hdulist_or_filename)
    elif isinstance(hdulist_or_filename, fits.HDUList):
        hdu_list = hdulist_or_filename
    else:
        raise ValueError('input must be a filename or HDUlist')

    # Create a copy of the PSF
    psf = copy.deepcopy(hdu_list)

    # Log instrument and detector names
    instrument = hdu_list[0].header['INSTRUME'].upper().strip()
    if instrument == 'WFI':
        aper_name = hdu_list[0].header['DETECTOR'] + '_FULL'
    else:
        aper_name = hdu_list[0].header['APERNAME'].upper()

    if instrument in ['MIRI', 'NIRSPEC']:
        raise ValueError("{}'s rotation is already included in STPSF and " "shouldn't be added again.".format(instrument))
    if instrument == 'WFI':
        raise ValueError('Rotation not necessary for {:} as pupil are aligned with detectors (to confirm).'.format(instrument))

    # Set rotation value if not already set by a keyword argument
    if rotate_value is None:
        aper = _get_default_siaf(instrument, aper_name)
        rotate_value = getattr(aper, 'V3IdlYAngle')  # the angle to rotate the PSF in degrees

    # If crop = True, then reshape must be False - so invert this keyword
    reshape = np.invert(crop)

    ext = 1  # edit the oversampled PSF (OVERDIST extension)

    psf_new = rotate(psf[ext].data, rotate_value, reshape=reshape)

    # Apply data to correct extensions
    psf[ext].data = psf_new

    # Set new header keyword
    psf[ext].header['ROTATION'] = (rotate_value, 'PSF rotated to match detector rotation')

    return psf
