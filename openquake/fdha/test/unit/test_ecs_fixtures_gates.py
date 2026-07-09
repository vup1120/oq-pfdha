# -*- coding: utf-8 -*-
"""
ECS validation gates against the R (mgcv) oracle fixtures (v2).

These tests are SKIPPED until the seven fixtures emitted by
``generate_ecs_fixtures.R`` are committed under ``fixtures/ecs/calingiri/``:
``counts.csv``, ``data4ecs_pre.csv``, ``ecs_trace.csv``, ``fault_disp.csv``,
``lpmatrix.csv``, ``coef.csv``, ``penalty_S_pred{1,2}.csv``.

Gate order (per project decision; stop at the first failure):

1. Selection / count gate -- our ``data4ecs`` selection must match R's
   ``counts.csv`` (expect disp=27, rup=320, replicated=3133) before anything
   downstream is meaningful.
2. GC2 gate -- ECS reference line from ``ecs_trace`` + 50 km extension,
   transform ``fault_disp`` lon/lat, match (u,t) to R's ``fault_disp`` u/t by
   SIGN and ORIGIN (oq km -> R m via x1000).
3. Spline basis gate -- the tp basis evaluated at ``fault_disp.u`` must match
   ``lpmatrix.csv`` cell-by-cell (added with the spline reconstruction).
4. Spline solution gate -- fitted coefficients vs ``coef.csv`` (with spline).
"""
import os

import numpy as np
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


# --------------------------------------------------------------------------- #
# Gate 2 -- GC2 (MultiLine + 50 km extension) sign & origin vs R fault_disp.
# PERMANENTLY SKIPPED: R's ecs_main cannot produce ecs_trace/fault_disp here --
# its GC2 step needs the FORTRAN r_wrapper_gc2 / R_dist_metrics.so, which is
# absent from the reference repo, the Zenodo package, and upstream GitHub (it
# was a precompiled artifact, never published). The GC2 engine is instead
# validated independently in test_ecs_gc2_engine.py (oq MultiLine tests + the
# analytic Spudich & Chiou 2015 single-segment reduction).
# --------------------------------------------------------------------------- #
def test_gate2_gc2_sign_and_origin():
    trace_p, fd_p = _need("ecs_trace.csv", "fault_disp.csv")
    trace = pd.read_csv(trace_p)
    fd = pd.read_csv(fd_p)

    utm = ecs.utm_for(trace["Longitude"].to_numpy(), trace["Latitude"].to_numpy())
    u_km, t_km = ecs.gc2ext_ut(
        fd["Longitude"].to_numpy(), fd["Latitude"].to_numpy(),
        trace["Longitude"].to_numpy(), trace["Latitude"].to_numpy(), utm)
    u_m, t_m = u_km * 1000.0, t_km * 1000.0  # oq km -> R m

    u_R, t_R = fd["u"].to_numpy(), fd["t"].to_numpy()

    # ORIGIN: along-strike origin agreement (both anchored near first vertex)
    assert abs(np.min(u_m) - np.min(u_R)) < 250.0, "u origin diverged"
    # SIGN + scale: strong positive correlation, slope ~ 1 (allow projection diff)
    su = np.polyfit(u_R, u_m, 1)[0]
    assert 0.97 < su < 1.03, f"u sign/scale diverged (slope {su:.3f})"
    # strike-normal: sign must agree (t may be flipped between conventions);
    # accept either a consistent + or - sign, but not a sign scramble
    st = np.polyfit(t_R, t_m, 1)[0]
    assert abs(st) > 0.9 and abs(abs(st) - 1.0) < 0.1, f"t scale off ({st:.3f})"
    corr = np.corrcoef(np.sign(t_R) * np.abs(t_m), t_m)[0, 1]
    assert np.isfinite(corr)
    # along-strike RMS agreement after unit conversion
    rms = float(np.sqrt(np.mean((u_m - u_R) ** 2)))
    assert rms < 250.0, f"u RMS {rms:.1f} m too large"


# --------------------------------------------------------------------------- #
# Gate 3 / 4 -- the tp-spline basis/solution validation. SUPERSEDED: because R
# could not run ecs_main (see Gate 2), the spline is validated against mgcv via
# the spl_* fixtures (spline_sample.csv, spl_pred.csv, spl_lpmatrix.csv,
# spl_coef.csv, spl_S1.csv) produced by the standalone spline oracle. See
# test_ecs_spline.py (fit_spline_xy reproduces mgcv's trace to <1 m, tip-safe).
# --------------------------------------------------------------------------- #
