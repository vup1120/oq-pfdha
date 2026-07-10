# -*- coding: utf-8 -*-
"""
Unit tests for the pure-Python LCP module (synthetic geometries, no raster
I/O). The Dijkstra core is additionally validated node-for-node against the
Thomas Ridgecrest outputs in test_lcp_reference_values.py.
"""
import numpy as np
import pytest

from openquake.fdha.calc.utils import lcp

pytestmark = pytest.mark.unit


# --------------------------------------------------------------------------- #
# route_through_grid: MCP_Geometric edge metric, hand-computed optima
# --------------------------------------------------------------------------- #
def test_route_uniform_grid_diagonal():
    # uniform cost: straight diagonal is optimal, cost = n_steps * v * sqrt(2)
    cost = np.full((3, 3), 2.0)
    path, total = lcp.route_through_grid(cost, (0, 0), (2, 2))
    np.testing.assert_array_equal(path, [[0, 0], [1, 1], [2, 2]])
    assert total == pytest.approx(2 * 2.0 * np.sqrt(2.0))


def test_route_cheap_channel():
    # expensive background, cost-1 channel along row 2: path follows the channel
    cost = np.full((5, 5), 100.0)
    cost[2, :] = 1.0
    path, total = lcp.route_through_grid(cost, (2, 0), (2, 4))
    np.testing.assert_array_equal(path, [[2, 0], [2, 1], [2, 2], [2, 3], [2, 4]])
    assert total == pytest.approx(4.0)


def test_route_geometric_prefers_cheap_diagonal():
    # (0,0)->(0,2): straight through the cost-5 cell = (1+5)/2+(5+1)/2 = 6;
    # dipping through row 1 = 2 * (1+1)/2 * sqrt(2) ~ 2.83 -> the dip wins
    cost = np.array([[1.0, 5.0, 1.0],
                     [1.0, 1.0, 1.0]])
    path, total = lcp.route_through_grid(cost, (0, 0), (0, 2))
    np.testing.assert_array_equal(path, [[0, 0], [1, 1], [0, 2]])
    assert total == pytest.approx(2.0 * np.sqrt(2.0))


def test_route_nonfinite_cells_are_barriers():
    cost = np.ones((3, 3))
    cost[:2, 1] = np.nan                      # wall with a gap at the bottom
    path, _total = lcp.route_through_grid(cost, (0, 0), (0, 2))
    rc = {tuple(p) for p in path}
    assert not rc & {(0, 1), (1, 1)}          # never crosses the wall
    # fully walled-off -> no path
    cost[:, 1] = np.inf
    with pytest.raises(ValueError, match="no traversable path"):
        lcp.route_through_grid(cost, (0, 0), (0, 2))


def test_route_rejects_out_of_grid_endpoints():
    with pytest.raises(ValueError, match="outside"):
        lcp.route_through_grid(np.ones((3, 3)), (0, 0), (3, 0))


# --------------------------------------------------------------------------- #
# GeoTransform: notebook MapCoord2Pixel / Pixel2Map conventions
# --------------------------------------------------------------------------- #
def test_geotransform_conventions():
    gt = lcp.GeoTransform(origin_x=1000.0, origin_y=5000.0,
                          pixel_dx=20.0, pixel_dy=-20.0)
    # int() truncation: anywhere inside pixel (3, 2) maps to (3, 2)
    assert gt.map_to_pixel(1079.9, 4941.0) == (3, 2)
    # pixel_to_map returns the pixel *corner* (no half-pixel shift)
    x, y = gt.pixel_to_map(3, 2)
    assert (x, y) == (1060.0, 4960.0)
    # roundtrip: corner maps back to the same pixel
    assert gt.map_to_pixel(x, y) == (3, 2)


# --------------------------------------------------------------------------- #
# rasterize_traces_xy: extent rules + burning
# --------------------------------------------------------------------------- #
def test_rasterize_extent_and_burn():
    px = 10.0
    trace = np.array([[100.0, 50.0], [160.0, 50.0]])       # horizontal, 60 m
    start, stop = (95.0, 45.0), (170.0, 55.0)
    cost, gt = lcp.rasterize_traces_xy([trace], px, start, stop,
                                       cost_fault=1.0, cost_background=100.0)
    # extent = union of trace bbox and endpoints, +- one pixel buffer
    assert gt.origin_x == 95.0 - px and gt.origin_y == 55.0 + px
    assert gt.pixel_dx == px and gt.pixel_dy == -px
    # endpoints fall inside the grid
    for pt in (start, stop):
        c, r = gt.map_to_pixel(*pt)
        assert 0 <= r < cost.shape[0] and 0 <= c < cost.shape[1]
    # exactly the pixels under the segment are burned, all in one row
    rows, cols = np.nonzero(cost == 1.0)
    assert len(rows) > 0 and np.all(rows == rows[0])
    xs, _ = gt.pixel_to_map(cols, rows)
    assert xs.min() >= 100.0 - px and xs.max() <= 160.0
    # everything else is background
    assert np.all(np.isin(cost, [1.0, 100.0]))


def test_rasterize_burn_is_continuous_along_diagonal():
    # what matters for the router: the burned trace has no gaps, so the LCP
    # from one end to the other never has to step on a background pixel
    px = 10.0
    trace = np.array([[0.0, 0.0], [100.0, 100.0]])
    cost, gt = lcp.rasterize_traces_xy([trace], px, (0.0, 0.0), (100.0, 100.0))
    path_xy, path_rc, _total = lcp.route_through_raster(
        cost, gt, (0.0, 0.0), (100.0, 100.0))
    assert np.all(cost[path_rc[:, 0], path_rc[:, 1]] == lcp.COST_FAULT)


# --------------------------------------------------------------------------- #
# shortcut smoothing (string pulling)
# --------------------------------------------------------------------------- #
def test_shortcut_straightens_uniform_background():
    # on uniform cost every monotone staircase is an equal-cost tie; the
    # shortcut must collapse the arbitrary Dijkstra pick to the straight chord
    cost = np.full((30, 90), 100.0)
    path, _ = lcp.route_through_grid(cost, (25, 3), (5, 80))
    out = lcp.shortcut_path(cost, path)
    np.testing.assert_array_equal(out, [[25, 3], [5, 80]])


def test_shortcut_does_not_cut_expensive_corners():
    # L-shaped cheap corridor: the chord across the corner samples cost-100
    # cells, so the shortcut must keep the corner node
    cost = np.full((40, 40), 100.0)
    cost[35, 2:38] = 1.0                     # horizontal leg
    cost[2:36, 37] = 1.0                     # vertical leg
    path, _ = lcp.route_through_grid(cost, (35, 2), (2, 37))
    out = lcp.shortcut_path(cost, path)
    assert len(out) >= 3                     # corner survives
    # all retained nodes sit on the cheap corridor
    assert np.all(cost[out[:, 0], out[:, 1]] == 1.0)
    # no multi-pixel shortcut segment cuts across the expensive background:
    # its interior samples (same rounding as shortcut_path) stay on cost 1
    # (single grid moves, e.g. the raw path's own diagonal corner step, are
    # legal 8-connected moves between cheap cells, not shortcuts)
    for a, b in zip(out[:-1].astype(float), out[1:].astype(float)):
        length = np.hypot(*(b - a))
        if length <= np.sqrt(2.0):
            continue
        frac = np.linspace(0.0, 1.0, max(int(np.ceil(length * 4)), 1) + 1)
        rr = np.round(a[0] + frac * (b[0] - a[0])).astype(int)
        cc = np.round(a[1] + frac * (b[1] - a[1])).astype(int)
        assert np.all(cost[rr, cc] == 1.0)


def test_shortcut_never_increases_cost_and_keeps_endpoints():
    rng = np.random.default_rng(42)
    cost = rng.uniform(1.0, 100.0, size=(40, 60))
    path, total = lcp.route_through_grid(cost, (2, 3), (35, 55))
    out = lcp.shortcut_path(cost, path)
    np.testing.assert_array_equal(out[0], path[0])
    np.testing.assert_array_equal(out[-1], path[-1])
    # retained nodes are a subsequence of the raw path
    as_set = {tuple(p) for p in path}
    assert all(tuple(p) in as_set for p in out)
    # continuous-metric cost of the shortcut path <= raw grid path cost
    def poly_cost(rc):
        c = 0.0
        for i in range(len(rc) - 1):
            a, b = rc[i].astype(float), rc[i + 1].astype(float)
            length = float(np.hypot(*(b - a)))
            n_s = max(int(np.ceil(length * 8)), 1)
            frac = np.linspace(0, 1, n_s + 1)
            rr = np.round(a[0] + frac * (b[0] - a[0])).astype(int)
            cc = np.round(a[1] + frac * (b[1] - a[1])).astype(int)
            c += float(np.mean(cost[rr, cc])) * length
        return c
    assert poly_cost(out) <= poly_cost(path) * (1.0 + 1e-6)


def test_lcp_from_traces_smooth_bridge_is_straight():
    # with smooth=True the inter-section bridge is a single straight segment;
    # the default (smooth=False, original-algorithm faithful) keeps the raw
    # staircase with many more nodes
    res_s = lcp.lcp_from_traces(_two_section_traces(), pixel_size=50.0,
                                smooth=True)
    res_r = lcp.lcp_from_traces(_two_section_traces(), pixel_size=50.0)
    assert len(res_s.lon) < len(res_r.lon)
    # smoothed path deviates < 1.5 px from the ideal line (S tip -> gap
    # bridge -> N tip); the raw staircase bridge deviates by much more
    x, y = res_s.utm.to_xy(res_s.lon, res_s.lat)
    pxy = np.column_stack([x, y])
    (lo1, la1), (lo2, la2) = _two_section_traces()
    t1 = np.column_stack(res_s.utm.to_xy(lo1, la1))
    t2 = np.column_stack(res_s.utm.to_xy(lo2, la2))
    ideal = np.vstack([t1, t2])              # traces + implied gap chord
    for p in pxy:
        assert _min_dist_to_polyline(p, ideal) < 1.5 * res_s.pixel_size


# --------------------------------------------------------------------------- #
# path decimation + endpoint selection
# --------------------------------------------------------------------------- #
def test_decimate_collinear():
    path = np.array([[0, 0], [0, 1], [0, 2], [1, 3], [2, 4], [2, 5]])
    out = lcp._decimate_collinear(path)
    np.testing.assert_array_equal(out, [[0, 0], [0, 2], [2, 4], [2, 5]])
    # geometry endpoints always kept; short paths unchanged
    np.testing.assert_array_equal(lcp._decimate_collinear(path[:2]), path[:2])


def test_farthest_endpoint_pair_orientation():
    # two N-S en-echelon traces; farthest endpoints are the extreme S and N
    # tips, oriented so the start is the southern one (strike ~ +y)
    t1 = np.column_stack([np.zeros(5), np.linspace(0.0, 4000.0, 5)])
    t2 = np.column_stack([np.full(5, 300.0), np.linspace(4500.0, 9000.0, 5)])
    a, b = lcp._farthest_endpoint_pair([t1, t2])
    np.testing.assert_allclose(a, [0.0, 0.0])
    np.testing.assert_allclose(b, [300.0, 9000.0])


# --------------------------------------------------------------------------- #
# lcp_from_traces: forward multi-fault entry point
# --------------------------------------------------------------------------- #
def _two_section_traces():
    # two roughly colinear ~N-S sections separated by a stepover (lon/lat)
    lat1 = np.linspace(-31.13, -31.09, 8)
    lat2 = np.linspace(-31.085, -31.05, 8)
    return [(np.full(8, 116.470), lat1), (np.full(8, 116.472), lat2)]


def test_lcp_from_traces_stepover():
    res = lcp.lcp_from_traces(_two_section_traces(), pixel_size=50.0)
    assert isinstance(res, lcp.LcpResult)
    # u is the cumulative arc length, strictly increasing from 0
    assert res.u[0] == 0.0 and np.all(np.diff(res.u) > 0)
    # one continuous line spanning the whole rupture (~8.9 km incl. stepover)
    span = res.u[-1]
    assert 8000.0 < span < 11000.0
    # the path stays on/near the traces (within a couple of pixels) and
    # bridges the inter-section gap
    x, y = res.utm.to_xy(res.lon, res.lat)
    pxy = np.column_stack([x, y])
    for lon_t, lat_t in _two_section_traces():
        tx, ty = res.utm.to_xy(lon_t, lat_t)
        for vx, vy in zip(tx, ty):
            d = _min_dist_to_polyline(np.array([vx, vy]), pxy)
            assert d < 2.5 * res.pixel_size


def test_lcp_x_l_ordering_and_range():
    res = lcp.lcp_from_traces(_two_section_traces(), pixel_size=50.0)
    lons = np.array([116.4705, 116.471, 116.4715])
    lats = np.array([-31.125, -31.0875, -31.055])       # S -> stepover -> N
    xl, L = res.x_l(lons, lats)
    assert L == pytest.approx(res.u[-1])
    assert np.all((xl >= 0.0) & (xl <= 1.0))
    # strike-consistent orientation: x/L grows from S to N (start = S tip)
    assert xl[0] < xl[1] < xl[2]
    assert xl[0] < 0.15 and xl[2] > 0.85


def test_lcp_from_traces_input_validation():
    with pytest.raises(ValueError, match=">=1 trace"):
        lcp.lcp_from_traces([(np.array([116.0]), np.array([-31.0]))])


def _min_dist_to_polyline(p, poly):
    dmin = np.inf
    for i in range(len(poly) - 1):
        a, b = poly[i], poly[i + 1]
        ab = b - a
        denom = float(ab @ ab)
        t = 0.0 if denom == 0.0 else float(np.clip((p - a) @ ab / denom, 0, 1))
        dmin = min(dmin, float(np.linalg.norm(p - (a + t * ab))))
    return dmin
