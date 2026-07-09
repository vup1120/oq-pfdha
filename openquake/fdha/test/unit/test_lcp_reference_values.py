# -*- coding: utf-8 -*-
"""
Reference-value validation of the LCP Dijkstra core against the Thomas
Ridgecrest outputs (see test/fixtures/lcp/README.md).

The committed `LCP_path_{PL,S2}.txt` were produced by the original notebook
with `skimage.graph.route_through_array(geometric=True, fully_connected=True)`
on the committed cost rasters; our scipy replication must reproduce them
node-for-node (verified exact at port time).
"""
import os

import numpy as np
import pytest

pytest.importorskip("PIL")
from PIL import Image  # noqa: E402  (test-only dependency)

from openquake.fdha.calc.utils import lcp

pytestmark = pytest.mark.unit

FIXDIR = os.path.join(os.path.dirname(__file__), os.pardir,
                      "fixtures", "lcp", "ridgecrest")

# GeoTIFF tags (test-only minimal reader; runtime lcp.py has no raster I/O)
MODEL_PIXEL_SCALE = 33550
MODEL_TIEPOINT = 33922


def load_geotiff(path):
    im = Image.open(path)
    arr = np.array(im, dtype=float)
    sx, sy, _sz = im.tag_v2[MODEL_PIXEL_SCALE]
    i, j, _k, x, y, _z = im.tag_v2[MODEL_TIEPOINT]
    gt = lcp.GeoTransform(origin_x=x - i * sx, origin_y=y + j * sy,
                          pixel_dx=sx, pixel_dy=-sy)
    return arr, gt


def _run_case(case):
    cost, gt = load_geotiff(os.path.join(
        FIXDIR, f"Images_4_LCP_Ridge_main_Ridgecrest_{case}_Cost_Raster_4_LCP.tif"))
    start = np.loadtxt(os.path.join(
        FIXDIR, f"Ridgecrest_{case}_LCP_start_point.txt"), delimiter="\t")
    stop = np.loadtxt(os.path.join(
        FIXDIR, f"Ridgecrest_{case}_LCP_end_points.txt"), delimiter="\t")
    ref = np.loadtxt(os.path.join(FIXDIR, f"LCP_path_{case}.txt"))
    path_xy, path_rc, total = lcp.route_through_raster(cost, gt, start, stop)
    return path_xy, path_rc, total, ref, cost, gt


@pytest.mark.parametrize("case", ["PL", "S2"])
def test_lcp_matches_thomas_reference_path(case):
    """Node-for-node match of the routed path (UTM pixel-corner coords)."""
    path_xy, _rc, _total, ref, _cost, _gt = _run_case(case)
    assert len(path_xy) == len(ref)
    # reference txt is written with %f (6 decimals); coords are exact
    # pixel-corner products so the match is exact up to that formatting
    np.testing.assert_allclose(path_xy, ref[:, :2], rtol=0, atol=1e-6)


@pytest.mark.parametrize("case", ["PL", "S2"])
def test_lcp_matches_thomas_reference_lonlat(case):
    """The lon/lat columns of the reference output (EPSG:32611 inverse)."""
    pyproj = pytest.importorskip("pyproj")
    path_xy, _rc, _total, ref, _cost, _gt = _run_case(case)
    tr = pyproj.Transformer.from_crs("epsg:32611", "epsg:4326", always_xy=True)
    lon, lat = tr.transform(path_xy[:, 0], path_xy[:, 1])
    np.testing.assert_allclose(lon, ref[:, 2], rtol=0, atol=5e-6)
    np.testing.assert_allclose(lat, ref[:, 3], rtol=0, atol=5e-6)


@pytest.mark.parametrize("case", ["PL", "S2"])
def test_lcp_cost_is_optimal_vs_reference(case):
    """Our accumulated cost equals the reference path's cost under the
    MCP_Geometric edge metric (both are the optimum)."""
    _xy, _rc, total, ref, cost, gt = _run_case(case)
    cols = np.round((ref[:, 0] - gt.origin_x) / gt.pixel_dx).astype(int)
    rows = np.round((ref[:, 1] - gt.origin_y) / gt.pixel_dy).astype(int)
    steps = np.diff(np.column_stack([rows, cols]), axis=0).astype(float)
    # the reference path moves one 8-connected step at a time
    assert np.abs(steps).max() == 1
    c = cost[rows, cols]
    ref_cost = float(np.sum(0.5 * (c[:-1] + c[1:]) * np.hypot(*steps.T)))
    assert total == pytest.approx(ref_cost, rel=1e-12)
