# -*- coding: utf-8 -*-
"""Moss et al. (2013) P(sr | Mw) against the redraw in GIRS-2022-05 Fig. 3.2.

GIRS-2022-05 (Moss et al., 2022) reproduces the four Moss et al. (2013)
curves in its Figure 3.2 "after Moss et al., 2013". They were digitised by
``digitize_girs_fig3_2.py`` into ``reference/girs_fig3_2.csv`` (axis
calibration residual <= 0.23 % in P and 0.006 in Mw). The primary check is
against the original paper's Figure 3 (``test_moss2013_fig3.py``); this one
records how the widely used report redraw relates to it.

Tolerances: 2 percentage points maximum, 1 point median. The plotted
stroke is ~6 px wide (~1 point in P); observed agreement is <= 1.1 points.

The redrawn reverse soft-soil curve is wrong: it is the logistic
z = -2.78 + 0.284 Mw, whereas the original Figure 3, its legend, the
report's own Eq. 3.5, its Appendix C MATLAB code and the library all give
z = -6.2548 + 0.8308 Mw; they differ by up to 31 points. The comparison is
kept as a strict xfail so that the report erratum stays on record.
"""
import csv
from pathlib import Path

import numpy as np
import pytest

from openquake.fdha.primary_surf_rup import Moss2013PrimarySR

pytestmark = pytest.mark.benchmark

REFERENCE = Path(__file__).parent / "reference" / "girs_fig3_2.csv"
VS30 = {"stiff": 760.0, "soft": 200.0}   # either side of the 600 m/s split
MAX_ABS = 0.02
MEDIAN_ABS = 0.01


def _curve(panel, curve):
    with open(REFERENCE) as f:
        rows = [r for r in csv.DictReader(f)
                if r["panel"] == panel and r["curve"] == curve]
    mw = np.array([float(r["mw"]) for r in rows])
    p = np.array([float(r["p"]) for r in rows])
    return mw, p


def _abs_diff(panel, curve):
    mw, p = _curve(panel, curve)
    model = np.asarray(Moss2013PrimarySR().get_prob(
        mw, style=panel, vs30=VS30[curve]))
    return np.abs(model - p)


@pytest.mark.parametrize("panel,curve", [
    ("reverse", "stiff"),
    ("strike-slip", "stiff"),
    ("strike-slip", "soft"),
])
def test_curve_matches_digitised_figure(panel, curve):
    d = _abs_diff(panel, curve)
    assert d.size >= 100
    assert np.median(d) < MEDIAN_ABS, (panel, curve, np.median(d))
    assert d.max() < MAX_ABS, (panel, curve, d.max())


@pytest.mark.parametrize("curve,a,b", [
    ("stiff", -11.4071, 1.8465),
    ("soft", -12.2908, 1.9520),
])
def test_strike_slip_coefficients_recovered_from_figure(curve, a, b):
    """A logistic fit to the digitised curve recovers the coded
    coefficients: the figure pins the strike-slip pair, not just the
    curve shape."""
    mw, p = _curve("strike-slip", curve)
    b_fit, a_fit = np.polyfit(mw, np.log(p / (1.0 - p)), 1)
    assert a_fit == pytest.approx(a, rel=0.02)
    assert b_fit == pytest.approx(b, rel=0.02)


@pytest.mark.xfail(strict=True, reason=(
    "erratum in GIRS-2022-05 Fig. 3.2 (reverse, soft soil): the redraw is "
    "z = -2.78 + 0.284 Mw, whereas Moss et al. (2013) Fig. 3, the report's "
    "Eq. 3.5 and the library are z = -6.2548 + 0.8308 Mw"))
def test_reverse_soft_curve_matches_digitised_figure():
    d = _abs_diff("reverse", "soft")
    assert d.max() < MAX_ABS, d.max()
