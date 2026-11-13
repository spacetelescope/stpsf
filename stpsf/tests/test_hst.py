"""Tests for the HST / WFC3 models that do not require the STPSF data package."""

import numpy as np
import poppy

from .. import hst


def _pupil_map(**kwargs):
    # 512 pixel grid with the pupil radius spanning 128 pixels, as in TinyTim
    wave = poppy.Wavefront(diam=4.8, npix=512)
    wave.pixelscale = 1.2 / 128
    return hst.HSTPupil(**kwargs).get_transmission(wave)


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
    cold = _pupil_map(cold_mask=True)
    assert np.isclose((ota * cold).sum(), 39552, rtol=3e-3)
    assert np.isclose(ota.sum() / (np.pi * 128**2), 1 - 0.33**2, atol=0.03)
