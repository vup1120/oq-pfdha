# -*- coding: utf-8 -*-
"""Pytest wrapper for the Takao et al. (2013) paper-figure reproductions.

Validates the Takao (2013) model chain against the paper's own worked
examples and against the Fig. 10 / Fig. 11 curves extracted from the
paper PDF's vector graphics (``extract_reference.py`` →
``reference/fig*_points.csv``). See README.md for the agreement summary.

Published text anchors:

- P1p(6.8) = 0.784, P2p(6.8) = 0.751, all eight placement counts of the
  Mw 6.8 worked example, rupture length 21.38 → 21 km;
- case (a): nu(0.01 m | Mw 6.6, R 30,000) = 1.2e-5 /yr;
- case (b): four-fault sum at 0.01 m = 7.0e-7 /yr.

Digitized-curve tolerances reflect the observed agreement (medians ~1-3%,
maxima ~5-10% of annual rate over up to 8 decades) plus headroom.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

from openquake.pfd.secondary_surf_displ import Takao2013SecondaryFD

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import reproduce_fig11b  # noqa: E402
import reproduce_fig11a  # noqa: E402
import reproduce_fig10  # noqa: E402


@pytest.fixture(scope="module")
def fig11b():
    return reproduce_fig11b.agreement()


@pytest.fixture(scope="module")
def fig11a():
    result, _, _ = reproduce_fig11a.agreement()
    return result


@pytest.fixture(scope="module")
def fig10():
    result, _, _ = reproduce_fig10.agreement()
    return result


# ------------------------------------------------------------- Fig 11 (b)

def test_11b_p1p_matches_paper(fig11b):
    assert fig11b["P1p_M6.8"]["computed"] == pytest.approx(
        reproduce_fig11b.PAPER_P1P, rel=1e-3)

def test_11b_sum_at_001m_matches_paper(fig11b):
    assert fig11b["sum_at_0.01m"]["computed"] == pytest.approx(
        reproduce_fig11b.PAPER_SUM_001M, rel=0.05)

def test_11b_fault2_fault3_cross(fig11b):
    assert fig11b["fault2_fault3_curves_cross"]

def test_11b_recurrence_scaling(fig11b):
    assert fig11b["fault1_over_fault2"] == pytest.approx(10.0, rel=1e-9)

def test_11b_digitized_curves(fig11b):
    for name, st in fig11b["fig11b_match"].items():
        assert st["n"] >= 20, name
        assert st["median_abs_relerr"] < 0.05, (name, st)
        assert st["max_abs_relerr"] < 0.10, (name, st)


# ------------------------------------------------------------- Fig 11 (a)

def test_11a_p1p_p2p_match_paper(fig11a):
    assert fig11a["P1p_M6.8"]["computed"] == pytest.approx(
        reproduce_fig11a.PAPER_P1P, rel=1e-3)
    assert fig11a["P2p_M6.8"]["computed"] == pytest.approx(
        reproduce_fig11a.PAPER_P2P, rel=1e-3)

def test_11a_placement_counts_match_paper(fig11a):
    assert fig11a["placement_counts_match_paper"]
    assert fig11a["rupture_length_M6.8_km"] == 21

def test_11a_text_anchor(fig11a):
    assert fig11a["nu_001m_M6.6_R30000"]["computed"] == pytest.approx(
        reproduce_fig11a.PAPER_NU_001_M66_R30000, rel=0.03)

def test_11a_digitized_curves(fig11a):
    for color in ("gray", "black"):
        for name, st in fig11a[f"fig11a_{color}_match"].items():
            assert st["n"] >= 20, (color, name)
            assert st["median_abs_relerr"] < 0.05, (color, name, st)
            assert st["max_abs_relerr"] < 0.12, (color, name, st)


# ------------------------------------------------------------- Fig 10

def test_10_ad_lognormal_exact(fig10):
    st = fig10["AD"]
    assert st["peak_x_ratio"] == pytest.approx(1.0, abs=0.01)
    assert st["peak_h_ratio"] == pytest.approx(1.0, abs=0.005)
    assert st["max_relerr_above_5pct_peak"] < 0.03

def test_10_d_convolution_exact(fig10):
    # no free scale: validates the AD-mass x gamma composition of
    # Takao2013PrimaryFD against the paper's plotted D density
    st = fig10["D"]
    assert st["peak_x_ratio"] == pytest.approx(1.0, abs=0.01)
    assert st["peak_h_ratio"] == pytest.approx(1.0, abs=0.005)
    assert st["max_relerr_above_5pct_peak"] < 0.03

def test_10_dad_gamma_shape(fig10):
    st = fig10["D_AD"]
    assert st["peak_x_ratio"] == pytest.approx(1.0, abs=0.03)
    assert st["max_shape_relerr_above_10pct_peak"] < 0.05


# ------------------------------------------- Takao2013SecondaryFD internals

@pytest.mark.parametrize("norm_disp_type,c90", [("MD", 0.55), ("AD", 1.9)])
@pytest.mark.parametrize("r", [0.0, 2.0, 5.0, 10.0, 20.0])
def test_eq15_16_are_the_90pct_level(norm_disp_type, c90, r):
    model = Takao2013SecondaryFD()
    level = c90 * np.exp(-0.17 * r)
    prob = model.get_prob_norm_displ(np.array([level]), r, norm_disp_type)
    assert prob[0, 0] == pytest.approx(0.10, abs=1e-12)


def test_invalid_norm_disp_type_raises():
    model = Takao2013SecondaryFD()
    with pytest.raises(ValueError):
        model.get_prob(np.array([0.1]), 6.8, 5.0, "XX")
