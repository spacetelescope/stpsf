import pytest
import webbpsf
from astropy.table import Table


@pytest.mark.remote_data
def test_load_mast_opd():
    """Test the machinery for loading an OPD from MAST.

    NOTE THIS TEST WILL DOWNLOAD A FILE FROM MAST.
    """

    nrc = webbpsf.NIRCam()
    nrc.load_wss_opd_by_date('2022-07-30')

    assert nrc.pupilopd[0].data.shape == (1024, 1024), "OPD data does not have expected dimensions"
    assert nrc.pupilopd[0].header['CORR_ID'] == 'O2022073001', "Missing expected correction ID in header"


@pytest.mark.remote_data
def test_load_mast_opd_larger_npix(npix=2048):
    """Test the machinery for loading an OPD from MAST, and rescaling to a larger pupil
    sampling for larger PSF simulations.

    NOTE THIS TEST WILL DOWNLOAD A FILE FROM MAST.
    """

    nrc = webbpsf.NIRCam()
    nrc.pupil = f'jwst_pupil_RevW_npix{npix}.fits.gz'
    nrc.load_wss_opd_by_date('2022-07-30')

    assert nrc.pupilopd[0].data.shape == (npix, npix), "OPD data does not have expected dimensions"
    assert nrc.pupilopd[0].header['CORR_ID'] == 'O2022073001', "Missing expected correction ID in header"
    osys = nrc.get_optical_system()
    assert osys.planes[0].amplitude.shape == (npix, npix), \
            "Optical system model pupil amplitude does not have expected dimensions"
    assert osys.planes[0].opd.shape == (npix, npix), "Optical system model pupil OPD does not have expected dimensions"


@pytest.mark.remote_data
def test_query_wfsc_images(test_download=False):
    """ Test that we can query MAST for images

    Usually we do not test the download because it's slow, and not
    necessary to test that we can download a large file in every CI run.
    """

    filetable = webbpsf.mast_wss.query_wfsc_images_latest()
    assert len(filetable) > 0, "Query should have nonzero results"
    assert isinstance(filetable, Table), "Query should return a table"

    filetable2 = webbpsf.mast_wss.query_wfsc_images_by_program(1160, 1)
    assert len(filetable2) > 0, "Query should have nonzero results"
    assert isinstance(filetable2, Table), "Query should return a table"

    if test_download:
        file_list = download_wfsc_images(1160, 1)
