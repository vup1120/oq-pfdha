"""Phase G / v4 test 7: per-site arithmetic on the map-mode fixture demo.

Runs the three sibling map-mode configs::

    openquake/fdha/test/fixtures/examples_archive/logic_tree_validation/map_mode/
    {single_bilinear_map, single_elliptical_map, blend_50_50_map}/job.ini

and asserts that on every active grid site ``i`` and every D0 level ``j``::

    lambda_blend(i, j) == 0.5 * lambda_A(i, j) + 0.5 * lambda_B(i, j)

to within 1e-12. This is the map-mode analogue of the single-site
two-branch weighted sum test.
"""
from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest

from openquake.fdha.logic_tree.driver import FdhaLogicTree


DEMO_ROOT = (
    Path(__file__).resolve().parents[5]
    / "openquake/fdha/test/fixtures/examples_archive/logic_tree_validation/map_mode"
)

RUNS = {
    "A": DEMO_ROOT / "single_bilinear_map",
    "B": DEMO_ROOT / "single_elliptical_map",
    "C": DEMO_ROOT / "blend_50_50_map",
}


def _rates_mean_path(run_dir: Path) -> Path:
    return run_dir / "out" / "aggregate" / "rates_mean.h5"


def _read_rates_mean(path: Path):
    with h5py.File(path, "r") as f:
        rates = np.asarray(f["rates_mean"][()], dtype=float)
        d0 = np.asarray(f["d0"][()], dtype=float)
        lons = np.asarray(f["site_lons"][()], dtype=float)
        lats = np.asarray(f["site_lats"][()], dtype=float)
    return rates, d0, lons, lats


@pytest.fixture(scope="module")
def run_three_scenarios():
    """Execute any scenario whose rates_mean.h5 is missing."""
    if not DEMO_ROOT.is_dir():
        pytest.skip(f"Map-mode demo root not found: {DEMO_ROOT}")

    for label, run_dir in RUNS.items():
        out = _rates_mean_path(run_dir)
        if out.is_file():
            continue
        ini = run_dir / "job.ini"
        if not ini.is_file():
            pytest.skip(f"Map-mode INI missing: {ini}")
        lt = FdhaLogicTree.from_ini(ini)
        lt.run(outdir=run_dir / "out")
        assert out.is_file(), f"driver failed to write {out}"


def test_map_per_site_arithmetic(run_three_scenarios):
    rA, d0A, lonsA, latsA = _read_rates_mean(_rates_mean_path(RUNS["A"]))
    rB, d0B, lonsB, latsB = _read_rates_mean(_rates_mean_path(RUNS["B"]))
    rC, d0C, lonsC, latsC = _read_rates_mean(_rates_mean_path(RUNS["C"]))

    # Grid consistency preconditions.
    assert rA.shape == rB.shape == rC.shape
    np.testing.assert_array_equal(d0A, d0B)
    np.testing.assert_array_equal(d0A, d0C)
    np.testing.assert_allclose(lonsA, lonsB)
    np.testing.assert_allclose(lonsA, lonsC)
    np.testing.assert_allclose(latsA, latsB)
    np.testing.assert_allclose(latsA, latsC)

    # On-fault sites must have non-zero hazard so the test is non-trivial.
    assert rA.max() > 0.0, "Scenario A produced zero hazard; demo is degenerate."
    assert rB.max() > 0.0, "Scenario B produced zero hazard; demo is degenerate."
    assert rC.max() > 0.0, "Scenario C produced zero hazard; demo is degenerate."

    expected = 0.5 * rA + 0.5 * rB
    diff = np.abs(rC - expected)
    max_abs = float(diff.max())
    assert max_abs <= 1e-12, (
        f"Per-site arithmetic failed: max|blend - 0.5*(A+B)| = {max_abs:.3e} > 1e-12"
    )
