# -*- coding: utf-8 -*-
"""
Validate the pure-Python penalized thin-plate spline (reconstruction of mgcv's
``s(u)`` tp basis + ``mvn`` family) and the end-to-end ``ecs_main`` pipeline.

Fidelity policy (project decision 2026-06-30): the spline is faithful to the
*method* (Wood 2003 tprs + mvn-coupled penalized GLS, mgcv ``scale.penalty``
normalisation) but not bit-identical to mgcv's internal ``repara``. The
committed tolerance below is set just above the observed convention residual so
a genuine basis/penalty bug fails, while mgcv's internal-normalisation residual
passes. The residual does NOT spike at the rupture tips (x/L-sensitive); the
max sits in the interior, and the propagated x/L deviation is < 1e-3.
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
    if not all(os.path.exists(p) for p in paths):
        pytest.skip("spline fixtures not present")
    return paths


# --------------------------------------------------------------------------- #
# fit_spline_xy reproduces mgcv's fitted ECS trace (within the committed
# convention tolerance), and is tip-safe.
# --------------------------------------------------------------------------- #
def test_fit_spline_reproduces_mgcv_trace_tip_safe():
    samp_p, pred_p = _need("spline_sample.csv", "spl_pred.csv")
    samp = pd.read_csv(samp_p)
    pred = pd.read_csv(pred_p)
    u = samp["u"].to_numpy()
    fault_len = u.max() - u.min()
    sp = 0.05 / fault_len
    predict, _ = ecs.fit_spline_xy(u, samp["x"].to_numpy(), samp["y"].to_numpy(), sp)
    px, py = predict(pred["u"].to_numpy())
    off = np.sqrt((px - pred["x"].to_numpy()) ** 2 + (py - pred["y"].to_numpy()) ** 2)
    ug = pred["u"].to_numpy()
    xl = (ug - ug.min()) / (ug.max() - ug.min())

    # committed convention tolerance: just above the observed ~0.84 m residual
    # (whitened natural parameterisation + mgcv scale.penalty). A real basis or
    # penalty bug pushes well past this; mgcv's internal-repara residual passes.
    assert off.max() < 1.1, f"trace offset {off.max():.2f} m exceeds convention tolerance"
    # propagated x/L deviation stays below the project bar (<3e-4)
    assert off.max() / fault_len < 3e-4
    # TIP-SAFE: tip offsets must not exceed the interior maximum (x/L is
    # tip-sensitive; a tip spike would signal a GC2/end-extension bug)
    tip = (xl < 0.1) | (xl > 0.9)
    assert off[tip].max() <= off[~tip].max() + 1e-9


def test_fit_spline_tracks_the_data():
    # the penalized fit should follow the data (residual well below data spread).
    samp_p, = _need("spline_sample.csv")
    samp = pd.read_csv(samp_p)
    u = samp["u"].to_numpy(); x = samp["x"].to_numpy(); y = samp["y"].to_numpy()
    sp = 0.05 / (u.max() - u.min())
    predict, beta = ecs.fit_spline_xy(u, x, y, sp)
    px, py = predict(u)
    # fit should track the data (residual far smaller than the data spread)
    assert np.std(x - px) < 0.5 * np.std(x)
    assert np.std(y - py) < 0.5 * np.std(y)


# --------------------------------------------------------------------------- #
# end-to-end ecs_main: converges, x/L in [0,1]
# --------------------------------------------------------------------------- #
@pytest.fixture
def calingiri_weighted():
    md = os.path.join(_FIX, "flatfile_measurements.csv")
    rp = os.path.join(_FIX, "flatfile_ruptures.csv")
    if not (os.path.exists(md) and os.path.exists(rp)):
        pytest.skip("Calingiri flatfiles not present")
    disp, rup = pd.read_csv(md), pd.read_csv(rp)
    data4ecs, _, _ = ecs.assemble_data4ecs(disp, rup)
    return ecs.weight_ecs_data(data4ecs), disp


def test_ecs_main_converges(calingiri_weighted):
    wt, _disp = calingiri_weighted
    res = ecs.ecs_main(wt, lambda_p=0.05, ecs_du=100.0, flt_max_ds=50.0, start="PCA")
    assert res.n_iter >= 1
    assert res.flt_ds < 50.0                       # converged below threshold
    assert len(res.u) >= 2 and res.u.max() > res.u.min()
    # ECS nodes lie near the rupture (sane lon/lat near Calingiri)
    assert 116.4 < np.mean(res.lon) < 116.5
    assert -31.2 < np.mean(res.lat) < -31.0


def test_ecs_main_xl_in_unit_interval(calingiri_weighted):
    wt, disp = calingiri_weighted
    res = ecs.ecs_main(wt, lambda_p=0.05, start="PCA")
    xl, L = res.x_l(disp["longitude_degrees"].to_numpy(),
                    disp["latitude_degrees"].to_numpy())
    assert L > 0
    assert np.all(xl >= 0.0) and np.all(xl <= 1.0)
    assert np.ptp(xl) > 0.5                         # points span a good range of the fault


def test_ecs_main_mrs_start_runs(calingiri_weighted):
    # forward multi-fault case uses mean-rupture-strike start (rupture-only).
    wt, _ = calingiri_weighted
    res = ecs.ecs_main(wt, lambda_p=0.05, start="MRS", rup_df=wt)
    assert res.u.max() > res.u.min()
