# -*- coding: utf-8 -*-
"""Pytest wrapper for the IAEA PFDHA exercise benchmark.

Runs every quantitative manifest entry end-to-end through the canonical
logic-tree driver and asserts the computed hazard curve against the
digitized team curve of the IAEA exercise paper (see README.md /
REFERENCE.md for provenance and the meaning of ``post_factor``).

Run from the repo root::

    pytest openquake/fdha/test/benchmark/IAEA/ -m benchmark -v

Individual cases::

    pytest openquake/fdha/test/benchmark/IAEA/ -m benchmark -k norcia -v
"""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import numpy as np
import pytest

os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from manifest import MANIFEST, Entry  # noqa: E402

QUANTITATIVE = [e for e in MANIFEST if e.assert_max_relerr is not None]


def _load_reference(entry: Entry):
    with (HERE / "reference" / entry.figure_csv).open() as f:
        rows = list(csv.DictReader(f))
    disp_m = np.array([float(r["disp_cm"]) for r in rows]) / 100.0
    ref = np.array([float(r[entry.column]) for r in rows])
    return disp_m, ref


@pytest.mark.benchmark
@pytest.mark.parametrize(
    "entry", QUANTITATIVE, ids=[f"{e.case}-{e.job}" for e in QUANTITATIVE])
def test_iaea_exercise_curve(entry: Entry, tmp_path):
    from openquake.fdha.logic_tree.driver import FdhaLogicTree

    ini = HERE / entry.case / f"job_{entry.job}.ini"
    result = FdhaLogicTree.from_ini(ini).run(outdir=tmp_path / "out")

    d0 = np.asarray(result.d0, dtype=float)
    rates = np.asarray(result.mean_rates, dtype=float)
    if rates.ndim == 2:
        rates = rates[0]
    assert np.any(rates > 0), "model chain produced an all-zero hazard curve"

    ref_d, ref = _load_reference(entry)
    np.testing.assert_allclose(d0, ref_d, rtol=1e-9)
    computed = rates * entry.post_factor

    ok = (ref > 0) & (computed > 0)
    if entry.assert_dmax_m is not None:
        ok &= ref_d <= entry.assert_dmax_m + 1e-12
    assert ok.any(), "no comparable points between computed and reference"

    rel = np.abs(computed[ok] / ref[ok] - 1.0)
    max_rel = float(rel.max())
    worst = ref_d[ok][int(rel.argmax())]
    assert max_rel <= entry.assert_max_relerr, (
        f"{entry.case}/{entry.job} vs {entry.column} in {entry.figure_csv}: "
        f"max |rel err| = {max_rel:.1%} at d = {worst:g} m exceeds "
        f"{entry.assert_max_relerr:.0%}"
    )
