# -*- coding: utf-8 -*-
"""Tests for the IAEA-exercise sensitivity-case demonstration jobs.

The exercise published no reference curves for the sensitivity cases, so
these jobs are validated by

1. structural checks — positive curve head, non-increasing hazard curve;
2. regression snapshots (``reference_snapshots/``, rtol 1e-4) — update
   with ``python run_sensitivity.py --update-snapshots`` after an
   intentional model change and review the diff;
3. physical-consistency relations across jobs:
   - Norcia sens2 (r = 2.4 km) carries more hazard than the base-case
     distributed site (r = 7.6 km) for the same source and chain;
   - Norcia sens3 (MVFS + NFS at the base site) carries more hazard than
     the base-case MVFS alone;
   - Kumamoto sens3 head equals rate x Takao P1p(5.8) (principal on-fault
     site with the full trace rupturing);
   - Kumamoto sens4 (r = 10 km) head is below the base-case distributed
     head (r = 5.22 km) despite summing four sources.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from sensitivity import SENSITIVITY_JOBS  # noqa: E402
from run_sensitivity import load_snapshot, run_job  # noqa: E402

pytestmark = pytest.mark.benchmark

_RESULTS = {}


def _get(case, job):
    tag = f"{case}_{job}"
    if tag not in _RESULTS:
        _RESULTS[tag] = run_job(case, job, HERE / "out" / tag)
    return _RESULTS[tag]


@pytest.mark.parametrize("sj", SENSITIVITY_JOBS,
                         ids=[f"{s.case}-{s.job}" for s in SENSITIVITY_JOBS])
def test_structure_and_snapshot(sj):
    d0, rates = _get(sj.case, sj.job)
    assert rates[0] > 0, "zero curve head — model chain failed"
    assert np.all(np.diff(rates) <= rates[:-1] * 1e-12 + 1e-30), \
        "hazard curve must be non-increasing"

    snap_d, snap_r = load_snapshot(f"{sj.case}_{sj.job}")
    assert snap_d is not None, \
        "missing snapshot — run run_sensitivity.py --update-snapshots"
    np.testing.assert_allclose(d0, snap_d, rtol=1e-9)
    np.testing.assert_allclose(rates, snap_r, rtol=1e-4, atol=1e-30)


def test_norcia_sens2_closer_site_higher_hazard():
    _, base = run_job("norcia", "distributed_Y03", HERE / "out" / "n_base")
    _, sens2 = _get("norcia", "distributed_Y03_sens2")
    assert sens2[0] > base[0]


def test_norcia_sens3_added_source_higher_hazard():
    _, base = run_job("norcia", "distributed_Y03", HERE / "out" / "n_base")
    _, sens3 = _get("norcia", "distributed_Y03_sens3")
    assert sens3[0] > base[0]


def test_kumamoto_sens3_head_is_rate_times_p1p():
    from openquake.fdha.primary_surf_rup import Takao2013PrimarySR

    _, rates = _get("kumamoto", "principal_T13_sens3")
    expected = 23.3e-5 * float(Takao2013PrimarySR().get_prob(5.8))
    assert rates[0] == pytest.approx(expected, rel=0.05)


def test_kumamoto_sens4_farther_site_lower_than_basecase():
    # same P11 chain and the same four sources: the r = 10 km site (sens4,
    # middle branches) must carry less hazard than the r = 5.22 km
    # base-case site (mean-collapsed branches, similar total rate)
    _, basecase = run_job("kumamoto", "distributed_P11_basecase",
                          HERE / "out" / "k_basecase")
    _, sens4 = _get("kumamoto", "distributed_P11_sens4")
    assert sens4[0] < basecase[0]
