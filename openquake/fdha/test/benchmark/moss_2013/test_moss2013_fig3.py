# -*- coding: utf-8 -*-
"""Moss et al. (2013) P(sr | Mw) against the paper's own Figure 3.

The four published curves (reverse / strike-slip x stiff / soft) were
digitised from the open-access copy by ``digitize_moss2013_fig3.py`` into
``reference/moss2013_fig3.csv`` (axis calibration residual <= 0.25 % in P
and 0.007 in Mw). The coefficients printed in the figure legend are checked
exactly in ``unit/test_moss2013_psr.py``; this benchmark confirms that the
plotted curves are those equations and that the library reproduces them.

Tolerances: 2 percentage points maximum, 1 point median. The source image
is a 150 ppi JPEG in which one pixel is ~0.4 points of P; observed
agreement is <= 0.9 points.
"""
import csv
from pathlib import Path

import numpy as np
import pytest

from openquake.fdha.primary_surf_rup import Moss2013PrimarySR

pytestmark = pytest.mark.benchmark

REFERENCE = Path(__file__).parent / "reference" / "moss2013_fig3.csv"
VS30 = {"stiff": 760.0, "soft": 200.0}   # either side of the 600 m/s split
MAX_ABS = 0.02
MEDIAN_ABS = 0.01

# Figure 3 legend: z = a + b Mw, P = 1 / (1 + exp(-z))
LEGEND = {
    ("reverse", "stiff"): (-13.9745, 2.1395),
    ("reverse", "soft"): (-6.2548, 0.8308),
    ("strike-slip", "stiff"): (-11.4071, 1.8465),
    ("strike-slip", "soft"): (-12.2908, 1.9520),
}


def _curve(panel, curve):
    with open(REFERENCE) as f:
        rows = [r for r in csv.DictReader(f)
                if r["panel"] == panel and r["curve"] == curve]
    return (np.array([float(r["mw"]) for r in rows]),
            np.array([float(r["p"]) for r in rows]))


@pytest.mark.parametrize("panel,curve", sorted(LEGEND))
def test_library_matches_digitised_curve(panel, curve):
    mw, p = _curve(panel, curve)
    model = np.asarray(Moss2013PrimarySR().get_prob(
        mw, style=panel, vs30=VS30[curve]))
    d = np.abs(model - p)
    assert d.size >= 50
    assert np.median(d) < MEDIAN_ABS, (panel, curve, np.median(d))
    assert d.max() < MAX_ABS, (panel, curve, d.max())


@pytest.mark.parametrize("panel,curve", sorted(LEGEND))
def test_plotted_curve_is_the_legend_equation(panel, curve):
    """A logistic refit of each digitised curve recovers its legend
    coefficients, so the figure and its printed equations agree."""
    a, b = LEGEND[(panel, curve)]
    mw, p = _curve(panel, curve)
    b_fit, a_fit = np.polyfit(mw, np.log(p / (1.0 - p)), 1)
    assert a_fit == pytest.approx(a, rel=0.02)
    assert b_fit == pytest.approx(b, rel=0.02)
