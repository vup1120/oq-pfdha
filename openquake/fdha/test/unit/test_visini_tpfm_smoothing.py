# -*- coding: utf-8 -*-
"""
Regression tests for the TPFm along-strike smoothing window of
Visini2025SecondaryFD.compute_tpfm_from_scaling.

The half-width of the smoothing window is r = 0.5 * distance / L(M); it
requires a fault-length estimate from the scaling relation. WC1994 (the
DEFAULT scaling model) exposes SRL(M) as ``get_surface_rupture_length``,
which the probe chain historically missed - so the smoothing was silently
never active with the default configuration.
"""
import numpy as np
import pytest

from openquake.fdha.secondary_surf_displ.visini2025 import Visini2025SecondaryFD

pytestmark = [pytest.mark.unit, pytest.mark.visini2025]

MAG = 7.0
DIP = 60.0
POS = 0.25  # on the sloping part of the mean-slip profile


def _manual_tpfm(distance_km):
    """Replicate the documented TPFm construction for WC1994/normal."""
    ad = 10 ** (-4.45 + 0.63 * MAG)   # WC94 normal AD (m)
    md = 10 ** (-5.90 + 0.89 * MAG)   # WC94 normal MD (m)
    L_km = 10 ** (-2.01 + 0.50 * MAG)  # WC94 normal SRL (km)

    x = np.linspace(0.0, 1.0, 1001)
    tri = np.maximum(1.0 - 2.0 * np.abs(x - 0.5), 0.0)
    tap = np.sqrt(np.sin(np.pi * x))
    T = 0.5 * (md * tri + 1.311 * ad * tap) * np.sin(np.radians(DIP))

    r = min(0.5 * distance_km / L_km, 0.5) if distance_km > 0 else 0.0
    if r == 0.0:
        return float(np.interp(POS, x, T))
    lo, hi = max(0.0, POS - r), min(1.0, POS + r)
    mask = (x >= lo) & (x <= hi)
    return float(np.mean(T[mask]))


def test_wc1994_smoothing_window_is_active():
    """A nonzero site distance must widen the window and change TPFm."""
    fd = Visini2025SecondaryFD()
    t0, _ = fd.compute_tpfm_from_scaling(
        MAG, style="normal", model="WC1994",
        norm_pos=POS, distance=0.0, dip=DIP)
    t5, _ = fd.compute_tpfm_from_scaling(
        MAG, style="normal", model="WC1994",
        norm_pos=POS, distance=5.0, dip=DIP)
    assert t5 != pytest.approx(t0), (
        "distance-dependent smoothing must alter TPFm away from the "
        "point value with the default WC1994 scaling model")


@pytest.mark.parametrize("distance_km", [0.0, 2.0, 5.0, 15.0])
def test_wc1994_tpfm_matches_manual_window_mean(distance_km):
    fd = Visini2025SecondaryFD()
    tpfm, sigma = fd.compute_tpfm_from_scaling(
        MAG, style="normal", model="WC1994",
        norm_pos=POS, distance=distance_km, dip=DIP)
    assert tpfm == pytest.approx(_manual_tpfm(distance_km), rel=1e-12)
    assert sigma == pytest.approx(0.33)  # WC94 normal AD log10 sigma


@pytest.mark.parametrize("model", ["THINGBAIJAM2017", "LEONARD2010"])
def test_other_scalers_still_smooth(model):
    """The pre-existing scalers keep their smoothing behaviour."""
    fd = Visini2025SecondaryFD()
    t0, _ = fd.compute_tpfm_from_scaling(
        MAG, style="normal", model=model,
        norm_pos=POS, distance=0.0, dip=DIP)
    t5, _ = fd.compute_tpfm_from_scaling(
        MAG, style="normal", model=model,
        norm_pos=POS, distance=5.0, dip=DIP)
    assert t5 != pytest.approx(t0)
