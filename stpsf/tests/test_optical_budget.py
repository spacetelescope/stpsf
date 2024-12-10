import stpsf


def test_visualize_wfe_budget():
    """Basic test that the visual optical budget functionality at least runs without raising an error or exception
    No actual checking of the output plots is performed.
    """
    nrc = stpsf.NIRCam()
    nrc.visualize_wfe_budget()
