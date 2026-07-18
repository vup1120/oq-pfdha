# -*- coding: utf-8 -*-
"""
ECS validation gate against the R (mgcv) oracle fixtures.

Gate 1 checks that our ``data4ecs`` selection matches R's ``counts.csv``
(disp=27, rup=320, replicated=3133) before anything downstream is
meaningful. The other validation layers live elsewhere: the GC2 engine is
validated in ``test_ecs_gc2_engine.py`` (oq MultiLine tests plus the
analytic Spudich & Chiou 2015 single-segment reduction) and the tp-spline
against the mgcv ``spl_*`` oracle fixtures in ``test_ecs_spline.py``.
"""
import os

import pandas as pd
import pytest

from openquake.fdha.calc.utils import ecs

pytestmark = pytest.mark.unit

_FIX = os.path.join(os.path.dirname(__file__), "..", "fixtures", "ecs", "calingiri")


def _need(*names):
    paths = [os.path.join(_FIX, n) for n in names]
    missing = [n for n, p in zip(names, paths) if not os.path.exists(p)]
    if missing:
        pytest.skip(f"R oracle fixture(s) not committed yet: {missing}")
    return paths if len(paths) > 1 else paths[0]


# --------------------------------------------------------------------------- #
# Gate 1 -- selection / counts
# --------------------------------------------------------------------------- #
def test_gate1_counts_match_oracle():
    counts_p, md, rp = _need("counts.csv", "flatfile_measurements.csv",
                             "flatfile_ruptures.csv")
    counts = pd.read_csv(counts_p)
    disp, rup = pd.read_csv(md), pd.read_csv(rp)
    data4ecs, wd, wr = ecs.assemble_data4ecs(disp, rup)
    wt = ecs.weight_ecs_data(data4ecs)
    assert len(wd) == int(counts.n_disp_kept[0]), "disp selection diverged"
    assert len(wr) == int(counts.n_rup_kept[0]), "rup selection diverged"
    assert len(data4ecs) == int(counts.n_data4ecs[0])
    assert len(wt) == int(counts.n_replicated[0]), "replication diverged"
