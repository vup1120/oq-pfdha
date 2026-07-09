# -*- coding: utf-8 -*-
"""
Unit tests for the deterministic (fixture-independent) ECS layers in
``openquake.fdha.calc.utils.ecs``: projection, geometry, curvature, and the
data-assembly / weighting that ports CalcEventCoordinateSystem.R.

These check the pieces that can be validated by construction / analytic
geometry, independently of the R ``lpmatrix``/``penalty_S`` spline fixtures.
"""
import os

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose

from openquake.fdha.calc.utils import ecs

pytestmark = pytest.mark.unit

_FIX = os.path.join(
    os.path.dirname(__file__), "..", "fixtures", "ecs", "calingiri"
)


# --------------------------------------------------------------------------- #
# Projection (pyproj UTM, must reproduce R's EPSG:32750 for Calingiri)
# --------------------------------------------------------------------------- #
def test_utm_zone_calingiri():
    # Calingiri ~116.47E -> UTM zone 50 (epsg_for_analysis=32750 == 50S)
    assert ecs.longlat2utm_zone(116.47) == 50


def test_utm_roundtrip_is_exact():
    lon = np.array([116.47, 116.48, 116.46])
    lat = np.array([-31.11, -31.12, -31.10])
    proj = ecs.utm_for(lon, lat)
    assert proj.zone == 50
    x, y = proj.to_xy(lon, lat)
    lon2, lat2 = proj.to_lonlat(x, y)
    assert_allclose(lon2, lon, atol=1e-9)
    assert_allclose(lat2, lat, atol=1e-9)


# --------------------------------------------------------------------------- #
# Curvature (analytic circle / straight line)
# --------------------------------------------------------------------------- #
def test_curvature_circle():
    R = 1000.0
    t = np.linspace(0, np.pi, 40)
    u = R * t  # arc length parametrisation
    cur = ecs.compute_curvature(u, R * np.cos(t), R * np.sin(t))
    assert_allclose(cur[1:-1], 1.0 / R, rtol=2e-3)
    assert np.isnan(cur[-1])  # faithful R endpoint quirk


def test_curvature_straight_line_is_zero():
    x = np.linspace(0, 5000, 30)
    cur = ecs.compute_curvature(x.copy(), x, 2.0 * x + 1.0)
    assert np.nanmax(np.abs(cur[1:-1])) < 1e-9


# --------------------------------------------------------------------------- #
# Geometry ports
# --------------------------------------------------------------------------- #
def test_rup_length_and_strike_xy():
    # horizontal segment of length 300 along +x -> strike 0
    x = np.array([0.0, 100.0, 300.0])
    y = np.zeros(3)
    assert ecs.rup_length_xy(x, y) == pytest.approx(300.0)
    assert ecs.rup_avg_strike_xy(x, y) == pytest.approx(0.0, abs=1e-12)
    # 45-degree segment
    assert ecs.rup_avg_strike_xy(
        np.array([0.0, 1.0]), np.array([0.0, 1.0])
    ) == pytest.approx(np.pi / 4)


# --------------------------------------------------------------------------- #
# Data assembly + weighting on the real Calingiri fixtures
# --------------------------------------------------------------------------- #
@pytest.fixture
def calingiri():
    md = os.path.join(_FIX, "flatfile_measurements.csv")
    rp = os.path.join(_FIX, "flatfile_ruptures.csv")
    if not (os.path.exists(md) and os.path.exists(rp)):
        pytest.skip("Calingiri flatfile fixtures not present")
    return pd.read_csv(md), pd.read_csv(rp)


def test_assemble_keeps_only_rank2consider_and_drops_outliers(calingiri):
    disp, rup = calingiri
    data4ecs, wd, wr = ecs.assemble_data4ecs(disp, rup)
    # only Principal disp rows, non-outlier (>= -900), non-NaN weight survive
    exp_disp = (
        disp["rank"].isin(ecs.RANK2CONSIDER)
        & (disp[ecs.FIELD_DISP_WT] >= -900.0)
    ).sum()
    assert len(wd) == exp_disp
    # only Principal rupture vertices survive
    assert len(wr) == int(rup["rank"].isin(ecs.RANK2CONSIDER).sum())
    # integer replication weights, all >= 1
    assert wd.min() >= 1 and wr.min() >= 1
    assert data4ecs[["Longitude", "Latitude", "wt"]].notna().all().all()


def test_weight_ecs_data_replicates_by_weight(calingiri):
    disp, rup = calingiri
    data4ecs, _, _ = ecs.assemble_data4ecs(disp, rup)
    wt = ecs.weight_ecs_data(data4ecs)
    # total rows == sum of (capped, >=1) integer weights
    n_rep = np.clip(np.round(data4ecs["wt"] / data4ecs["wt"].min()),
                    1, ecs.R_THRES_MAX).astype(int)
    assert len(wt) == int(n_rep.sum())
    assert (wt["wt"] == 1).all()


def test_forward_case_uses_rupture_only(calingiri):
    disp, rup = calingiri
    data4ecs, wd, wr = ecs.assemble_data4ecs(disp, rup, use_disp=False)
    assert len(wd) == 0
    assert len(wr) == int(rup["rank"].isin(ecs.RANK2CONSIDER).sum())
    assert "RUP_ID" in data4ecs.columns


# --------------------------------------------------------------------------- #
# Start solution (PCA / MRS) + initial coarse ECS
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("theta_deg", [0.0, 30.0, -50.0, 80.0])
def test_pca_recovers_principal_strike(theta_deg):
    rng = np.random.default_rng(0)
    a = rng.normal(0, 10.0, 400)
    b = rng.normal(0, 0.5, 400)
    th = np.deg2rad(theta_deg)
    c, s = np.cos(th), np.sin(th)
    x, y = c * a - s * b, s * a + c * b
    rot = ecs.pca_strike(x, y, np.ones_like(x))
    # -rot is the principal-axis angle; compare mod pi to theta
    err = np.angle(np.exp(1j * 2 * (-rot - th))) / 2
    assert abs(np.rad2deg(err)) < 1.0


def test_initial_ecs_is_collinear_and_spans_data():
    rng = np.random.default_rng(1)
    x = rng.normal(0, 10, 200)
    y = 0.3 * x + rng.normal(0, 0.5, 200)
    rot = ecs.pca_strike(x, y, np.ones_like(x))
    ex, ey = ecs.initial_ecs(x - x.mean(), y - y.mean(), rot, n_pt=10)
    assert len(ex) == 10
    A = np.column_stack([ex, np.ones_like(ex)])
    res = ey - A @ np.linalg.lstsq(A, ey, rcond=None)[0]
    assert np.max(np.abs(res)) < 1e-6


def test_mrs_strike_on_calingiri_is_finite(calingiri):
    _disp, rup = calingiri
    rup = rup.rename(columns={"longitude_degrees": "Longitude",
                              "latitude_degrees": "Latitude"})
    rup = rup[rup["rank"].isin(ecs.RANK2CONSIDER)]
    utm = ecs.utm_for(rup["Longitude"].to_numpy(), rup["Latitude"].to_numpy())
    rot = ecs.mrs_strike(rup, utm)
    assert np.isfinite(rot)
    assert -np.pi <= rot <= np.pi
