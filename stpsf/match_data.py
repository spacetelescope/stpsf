# Functions to match or fit PSFs to observed JWST data
import astropy
import astropy.io.fits as fits
import pysiaf

import stpsf


def setup_sim_to_match_file(filename_or_HDUList_or_datamodel, verbose=True, plot=False, choice='closest'):
    """Setup a stpsf Instrument instance matched to a given dataset

    The input can flexibly be either:
     - a string filename, either of a JWST FITS file or a Roman ASDF file
     - a FITS HDUList instance, for JWST data
     - a roman_datamodels.DataModel instance, for Roman data

    Parameters
    ----------
    filename_or_HDUlist_or_datamodel : str or astropy.io.fits HDUList or Roman Datamodel
        file to load
    verbose : bool
        be more verbose?
    plot : bool
        plot?
    choice : string
        Method to choose which OPD file to use, e.g. 'before', 'after', or 'closest', for
        JWST data. Not currently relevant for Roman.

    Returns
    -------
    an STSPF instrument instance for one of JWST NIRCam, NIRISS, NIRSpec, MIRI, FGS, or Roman WFI,
    configured with filter, detector, and other relevant properties configured to match the input data.

    """

    # Handle checking for JWST and Roman data in a flexible way that allows optional dependencies to be absent
    # we may not be running in an environment that has both the JWST and Roman package collections installed

    try:
        import roman_datamodels as rdm
        _HAVE_ROMAN = True
    except ImportError:
        _HAVE_ROMAN = False
        rdm = None
    try:
        import jwst.datamodels
        _HAVE_JWST = True
    except ImportError:
        _HAVE_JWST = False

    # If the input is a Roman data model, or a filename for an ASDF file, then try to set up a Roman sim
    if (_HAVE_ROMAN and (isinstance(filename_or_HDUList_or_datamodel, rdm.DataModel) or
        (isinstance(filename_or_HDUList_or_datamodel, str) and filename_or_HDUList_or_datamodel.endswith('asdf')))):
        return _setup_sim_to_match_file_roman(filename_or_HDUList_or_datamodel, verbose=verbose, plot=plot, choice=choice)
    elif _HAVE_JWST:
        return _setup_sim_to_match_file_jwst(filename_or_HDUList_or_datamodel, verbose=verbose, plot=plot, choice=choice)
    else:
        raise ImportError("Neither Roman nor JWST data models are available.")


def _setup_sim_to_match_file_jwst(model, verbose=True, plot=False, choice='closest'):
    """JWST implementation for setup_sim_to_match_file
    """
    import jwst.datamodels as jdm
    model = jdm.JwstDataModel(model)

    exptype = model.meta.exposure.type
    filt = model.meta.instrument.filter
    pupil = model.meta.instrument.pupil
    apername = model.meta.aperture.name
    channel = model.meta.instrument.channel
    band = model.meta.instrument.band
    inst = stpsf.instrument(model.meta.instrument.name)

    if inst.name == 'MIRI' and exptype == 'MIR_MRS':
        print("MIRI MRS exposure detected; configuring for IFU mode")
        inst.mode = 'IFU'
        # There is no FILTER keyword for MRS, so don't set filter to anything.
    elif inst.name == 'MIRI' and filt == 'P750L':
        # stpsf doesn't model the MIRI LRS prism spectral response
        print('Please note, stpsf does not currently model the LRS spectral response. Setting filter to F770W instead.')
        inst.filter = 'F770W'
    elif (inst.name == 'NIRCam') and (pupil[0] == 'F') and (pupil[-1] in ['N', 'M']):
        # These NIRCam filters are physically in the pupil wheel, but still act as filters.
        # Grab the filter name from the PUPIL keyword in this case.
        inst.filter = pupil
    elif (inst.name == 'NIRISS') and (filt == 'CLEAR'):
        # For NIRISS, 6 out of 12 filters are in the pupil wheel, which mean if FILTER=CLEAR,
        # PUPIL keyword will point to the actual filter. [S. T. Sohn Feb 13, 2025]
        inst.filter = pupil
    else:
        inst.filter = filt
    inst.set_position_from_aperture_name(apername)

    dateobs = astropy.time.Time(model.meta.observation.date + 'T' + model.meta.observation.time)
    inst.load_wss_opd_by_date(dateobs, verbose=verbose, plot=plot, choice=choice)

    # per-instrument specializations
    if inst.name == 'NIRCam':
        if pupil.startswith('MASK'):
            if pupil == 'MASKBAR':
                # the FITS header is just 'BAR' but the value needed in stpsf is either
                # 'MASKLWB' or "MASKSWB' depending on channel.
                inst.pupil_mask = 'MASKLWB' if channel == 'LONG' else 'MASKSWB'
            else:
                inst.pupil_mask = pupil
            if 'CORONMSK' in model.meta.instrument:
                inst.image_mask = model.meta.instrument.coronmsk.replace('MASKA', 'MASK')  # note, have to modify the value slightly for
                # consistency with the labels used in stpsf
            # The apername keyword is not always correct for cases with dual-channel coronagraphy
            # in some such cases, APERNAME != PPS_APER. Let's ensure we have the proper apername for this channel:
            apername = get_nrc_coron_apname(model)
            inst.set_position_from_aperture_name(apername)

        elif pupil != 'CLEAR' and not pupil.startswith('F'):  # no action needed for these
            # note that filters in the pupil wheel were handled already above
            inst.pupil_mask = pupil

    elif inst.name == 'MIRI':
        if exptype == 'MIR_MRS':
            ch = channel
            band_lookup = {'SHORT': 'A', 'MEDIUM': 'B', 'LONG': 'C'}
            inst.band = str(ch) + band_lookup[band]

        elif inst.filter in ['F1065C', 'F1140C', 'F1550C']:
            inst.image_mask = 'FQPM' + inst.filter[1:5]
        elif inst.filter == 'F2300C':
            inst.image_mask = 'LYOT2300'
        elif filt == 'P750L':
            inst.pupil_mask = 'P750L'

        if apername == 'MIRIM_SLIT':
            inst.image_mask = 'LRS slit'

    elif inst.name == 'NIRISS':
        if pupil == 'NRM': # else could be CLEARP for KPI observations
            inst.pupil_mask = 'MASK_NRM'

    # TODO add other per-instrument keyword checks

    if verbose:
        print(
            f"""
Configured simulation instrument for:
    Instrument: {inst.name}
    Filter: {inst.filter}
    Detector: {inst.detector}
    Apername: {inst.aperturename}
    Det. Pos.: {inst.detector_position} {'in subarray' if "FULL" not in inst.aperturename else ""}
    Image plane mask: {inst.image_mask}
    Pupil plane mask: {inst.pupil_mask}
    """
        )

    return inst


def _setup_sim_to_match_file_roman(filename_or_datamodel, verbose=True, plot=False, choice='closest'):
    """Roman implementation for setup_sim_to_match_file
    """
    import roman_datamodels as rdm
    if isinstance(filename_or_datamodel, str):
        if verbose:
            print(f'Setting up sim to match {filename_or_datamodel}')
        datamodel = rdm.open(filename_or_datamodel)
    elif isinstance(filename_or_datamodel, rdm.DataModel):
        datamodel = filename_or_datamodel
        if verbose:
            print('Setting up sim to match provided Roman Datamodel object')
    else:
        raise ValueError("Don't know how to load a Roman data model from that type of input")

    inst = stpsf.WFI()
    inst.detector = datamodel.meta.instrument.detector

    filter_or_disperser = datamodel.meta.instrument.optical_element
    # special case: STPSF has both  'GRISM0' and 'GRISM1' to handle the 0th undispersed light and 1st order spectra.
    # if we are given a grism file, by default let's assume the user wants to model the 1st order spectra
    if filter_or_disperser.upper() == 'GRISM':
        filter_or_disperser = 'GRISM1'
    inst.filter = filter_or_disperser

    if verbose:
        print(
            f"""
Configured simulation instrument for:
    Instrument: {inst.name}
    Filter: {inst.filter}
    Detector: {inst.detector}
    Apername: {inst.aperturename}
    Det. Pos.: {inst.detector_position} {'in subarray' if "FULL" not in inst.aperturename else ""}
    """
        )

    return inst


def get_nrc_coron_apname(input):
    """Get NIRCam coronagraph aperture name from header or data model

    Handles edge cases for dual-channel coronagraphy.

    By Jarron Leisenring originally in stpsf_ext, copied here by permission

    Parameters
    ==========
    input : fits.header.Header or datamodels.DataModel
        Input header or data model
    """

    if isinstance(input, (fits.header.Header)):
        # Aperture names
        apname = input['APERNAME']
        apname_pps = input['PPS_APER']
        subarray = input['SUBARRAY']
    else:
        # Data model meta info
        meta = input.meta

        # Aperture names
        apname = meta.aperture.name
        apname_pps = meta.aperture.pps_name
        subarray = meta.subarray.name

    # print(apname, apname_pps, subarray)

    # No need to do anything if the aperture names are the same
    # Also skip if MASK not in apname_pps
    if ((apname == apname_pps) or ('MASK' not in apname_pps)) and ('400X256' not in subarray):
        apname_new = apname
    else:
        # Should only get here if coron mask and apname doesn't match PPS
        apname_str_split = apname.split('_')
        sca = apname_str_split[0]
        image_mask = get_nrc_coron_mask_from_pps_apername(apname_pps)

        # Get subarray info
        # Sometimes apname erroneously has 'FULL' in it
        # So, first for subarray info in apname_pps
        if ('400X256' in apname_pps) or ('400X256' in subarray):
            apn0 = f'{sca}_400X256'
        elif 'FULL' in apname_pps:
            apn0 = f'{sca}_FULL'
        else:
            apn0 = sca

        apname_new = f'{apn0}_{image_mask}'

        # Append filter or NARROW if needed
        pps_str_arr = apname_pps.split('_')
        last_str = pps_str_arr[-1]
        # Look for filter specified in PPS aperture name
        if ('_F1' in apname_pps) or ('_F2' in apname_pps) or ('_F3' in apname_pps) or ('_F4' in apname_pps):
            # Find all instances of "_"
            inds = [pos for pos, char in enumerate(apname_pps) if char == '_']
            # Filter is always appended to end, but can have different string sizes (F322W2)
            filter = apname_pps[inds[-1] + 1:]
            apname_new += f'_{filter}'
        elif last_str == 'NARROW':
            apname_new += '_NARROW'
        elif ('TAMASK' in apname_pps) and ('WB' in apname_pps[-1]):
            apname_new += '_WEDGE_BAR'
        elif ('TAMASK' in apname_pps) and (apname_pps[-1] == 'R'):
            apname_new += '_WEDGE_RND'

    # print(apname_new)

    # If apname_new doesn't exit, we need to fall back to apname
    # even if it may not fully make sense.
    if apname_new in pysiaf.Siaf('NIRCam').apernames:
        return apname_new
    else:
        return apname


def get_nrc_coron_mask_from_pps_apername(apname_pps):
    """Get NIRCam coronagraph mask name from PPS aperture name

    The PPS aperture name is of the form:
        NRC[A/B][1-5]_[FULL]_[TA][MASK]
    where MASK is the name of the coronagraphic mask used.

    For target acquisition apertures the mask name can be
    prependend with "TA" (eg., TAMASK335R).

    Return '' if MASK not in input aperture name.
    """

    if 'MASK' not in apname_pps:
        return ''

    pps_str_arr = apname_pps.split('_')
    for s in pps_str_arr:
        if 'MASK' in s:
            image_mask = s
            break

    # Special case for TA apertures
    if 'TA' in image_mask:
        # Remove TA from mask name
        image_mask = image_mask.replace('TA', '')

        # Remove FS from mask name
        if 'FS' in image_mask:
            image_mask = image_mask.replace('FS', '')

        # Remove trailing S or L from LWB and SWB TA apertures
        if ('WB' in image_mask) and (image_mask[-1] == 'S' or image_mask[-1] == 'L'):
            image_mask = image_mask[:-1]

    return image_mask
