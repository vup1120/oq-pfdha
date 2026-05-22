"""Phase G / v4 test 8: inversion consistency for the map-mode demo.

For the ``single_bilinear_map`` scenario, read ``rates_mean.h5`` and
manually invert each site's mean rate curve at the configured return
period using :func:`get_map_from_curves`, then compare to the values
stored in ``displacement_map_mean.csv``. The two must agree to 1e-10 m
on every site.
"""
from __future__ import annotations

import csv
from pathlib import Path

import h5py
import numpy as np
import pytest

from openquake.fdha.calc.utils.interpolation import get_map_from_curves
from openquake.fdha.logic_tree.driver import FdhaLogicTree


DEMO_DIR = (
    Path(__file__).resolve().parents[5]
    / "openquake/fdha/test/fixtures/examples_archive/logic_tree_validation"
    / "map_mode/single_bilinear_map"
)


def _read_rates_mean(path: Path):
    with h5py.File(path, "r") as f:
        rates = np.asarray(f["rates_mean"][()], dtype=float)
        d0 = np.asarray(f["d0"][()], dtype=float)
        lons = np.asarray(f["site_lons"][()], dtype=float)
        lats = np.asarray(f["site_lats"][()], dtype=float)
    return rates, d0, lons, lats


def _read_map_csv(path: Path):
    rp = None
    rows: list[tuple[int, float, float, int, float]] = []
    with path.open() as f:
        # Parse the '# return_period = T' header first.
        lines = f.readlines()
    for line in lines:
        if line.startswith("#"):
            if "return_period" in line:
                rp = float(line.split("=")[1].strip())
            continue
        break
    non_comment = [ln for ln in lines if not ln.startswith("#")]
    reader = csv.reader(non_comment)
    header = next(reader)
    for row in reader:
        rows.append(
            (int(row[0]), float(row[1]), float(row[2]), int(row[3]), float(row[4]))
        )
    return rp, header, rows


@pytest.fixture(scope="module")
def ensure_demo_run():
    if not DEMO_DIR.is_dir():
        pytest.skip(f"Map-mode demo dir missing: {DEMO_DIR}")
    rates_h5 = DEMO_DIR / "out" / "aggregate" / "rates_mean.h5"
    map_csv = DEMO_DIR / "out" / "aggregate" / "displacement_map_mean.csv"
    if rates_h5.is_file() and map_csv.is_file():
        return rates_h5, map_csv
    ini = DEMO_DIR / "job.ini"
    if not ini.is_file():
        pytest.skip(f"Map-mode INI missing: {ini}")
    lt = FdhaLogicTree.from_ini(ini)
    lt.run(outdir=DEMO_DIR / "out")
    assert rates_h5.is_file() and map_csv.is_file()
    return rates_h5, map_csv


def test_map_inversion_consistency(ensure_demo_run):
    rates_h5, map_csv = ensure_demo_run
    rates, d0, lons, lats = _read_rates_mean(rates_h5)
    rp, header, rows = _read_map_csv(map_csv)

    assert rp is not None and rp > 0, "return_period header missing in displacement map CSV"
    csv_lons = np.array([r[1] for r in rows])
    csv_lats = np.array([r[2] for r in rows])
    csv_displ = np.array([r[4] for r in rows])

    # Must share the exact site ordering with rates_mean.h5.
    np.testing.assert_allclose(csv_lons, lons)
    np.testing.assert_allclose(csv_lats, lats)

    pex = 1.0 / rp
    manual = get_map_from_curves(d0, rates, pex)
    diff = np.abs(manual - csv_displ)
    max_abs = float(diff.max())
    assert max_abs <= 1e-10, (
        f"Inversion mismatch: max|manual - csv| = {max_abs:.3e} > 1e-10"
    )

    # And at least one site should have non-zero displacement so the check
    # is not vacuous.
    assert float(csv_displ.max()) > 0.0, (
        "All sites report zero displacement; demo is degenerate for inversion."
    )
