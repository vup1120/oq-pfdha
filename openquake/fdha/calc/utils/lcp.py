# -*- coding: utf-8 -*-
"""
Least Cost Path (LCP) representative reference-line construction.

Pure-Python port of the Thomas/Milliner LCP workflow
(``pyLCP_Field_Maps_v3.ipynb`` / ``pyLCP_Cost2Fault.ipynb``; see
test/usercase/LCP/Thomas_LCP) used to build a single representative
principal-rupture reference line for a multi-section / multi-fault rupture,
from which the normalized along-strike coordinate ``x/L`` is derived. It is a
companion to the ECS builder (:mod:`openquake.fdha.calc.utils.ecs`): both
produce a reference line with the same ``x_l`` interface and are selected via
the ``reference_line_method`` config switch ('ecs' | 'lcp').

Reference workflow (notebooks):

1. Rasterize the rupture traces onto a UTM grid (``shp_to_raster``): extent =
   union of the trace bounding box and the LCP start/end points, buffered by
   one pixel; GDAL corner-anchored geotransform.
2. Turn the rupture-map raster into a cost raster
   (``Calc_Cost_raster_from_rupture_map``): fault pixels cost 1, background
   cost 100 (defaults).
3. Run a geometric least-cost path between the start and end points
   (``skimage.graph.route_through_array(..., geometric=True,
   fully_connected=True)``): 8-connected moves, edge cost = mean of the two
   pixel costs x Euclidean step length.
4. Convert the path pixel indices back to UTM via the geotransform
   (``Pixel2Map``: map = origin + index * pixel, i.e. pixel corners).

Runtime constraints (same as ecs.py): pure Python - numpy / scipy / pyproj
only. No GDAL, no scikit-image, no rasterio. The Dijkstra core replicates
skimage ``MCP_Geometric`` via ``scipy.sparse.csgraph``; it is validated
against the committed Thomas Ridgecrest outputs (cost rasters + picked
endpoints + reference paths) in test/fixtures/lcp/ridgecrest. GC2 for x/L is
taken from oq-engine ``MultiLine`` (via :func:`ecs.gc2ext_ut`), NOT
reimplemented.

One documented optional post-step beyond the notebooks: the forward builder
(:func:`lcp_from_traces`) can apply a cost-checked shortcut pass
(:func:`shortcut_path`) to the raw grid path. Across a uniform-cost gap the
8-connected optimum is a degenerate tie (every monotone staircase costs the
same, and all overcharge the straight chord by up to ~8% - the grid-metric
artifact; the original workflow leaves whatever arbitrary staircase the
router returns), so the pass collapses such bridges to the straight line
while provably never entering more expensive cells. The DEFAULT is
``smooth=False`` - exactly the original ``route_through_array`` behaviour
(project decision 2026-07-02); pass ``smooth=True`` to straighten tie
bridges.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra

from openquake.fdha.calc.utils import ecs as _ecs

# --- defaults, from pyLCP_Field_Maps_v3.ipynb __main__ ------------------------
COST_FAULT = 1.0          # Cost_val_4_fault
COST_BACKGROUND = 100.0   # Cost_val_4_nonfault
PIXEL_SIZE_M = 100.0      # forward-model default (notebook field case uses 20 m)


# =============================================================================
# Geotransform helpers -- GDAL corner-anchored affine, faithful to the notebook
# =============================================================================
@dataclass(frozen=True)
class GeoTransform:
    """GDAL-style north-up geotransform (no rotation terms).

    ``origin_x/origin_y`` is the map coordinate of the *top-left corner* of
    pixel (0, 0); ``pixel_dy`` is negative for top-down rasters.
    """
    origin_x: float
    origin_y: float
    pixel_dx: float
    pixel_dy: float

    def map_to_pixel(self, x, y):
        """Map coords -> integer pixel (col, row). Ports ``MapCoord2Pixel``
        (int() truncation, exactly as the notebook)."""
        col = int((x - self.origin_x) / self.pixel_dx)
        row = int((y - self.origin_y) / self.pixel_dy)
        return col, row

    def pixel_to_map(self, col, row):
        """Pixel (col, row) -> map coords of the pixel *corner*. Ports
        ``Pixel2Map`` (origin + index * pixel, no half-pixel shift)."""
        x = self.origin_x + np.asarray(col, dtype=float) * self.pixel_dx
        y = self.origin_y + np.asarray(row, dtype=float) * self.pixel_dy
        return x, y


# =============================================================================
# Geometric least-cost path -- replicates skimage MCP_Geometric on a grid
# =============================================================================
def route_through_grid(cost: np.ndarray, start_rc, stop_rc):
    """Least-cost path across ``cost`` between two (row, col) cells.

    Replicates ``skimage.graph.route_through_array(cost, start, stop,
    geometric=True, fully_connected=True)`` with scipy: the grid is an
    8-connected graph whose edge weight between adjacent cells a, b is
    ``(cost[a] + cost[b]) / 2 * step`` with ``step`` = 1 (orthogonal) or
    sqrt(2) (diagonal). Cells with non-finite cost are impassable.

    :param cost: 2-D array of non-negative traversal costs.
    :param start_rc: (row, col) of the start cell.
    :param stop_rc: (row, col) of the end cell.
    :returns: ``(path_rc, total_cost)`` where ``path_rc`` is an (n, 2) int
        array of (row, col) from start to stop inclusive.
    """
    cost = np.asarray(cost, dtype=float)
    nr, nc = cost.shape
    for name, (r, c) in (("start", start_rc), ("stop", stop_rc)):
        if not (0 <= r < nr and 0 <= c < nc):
            raise ValueError(f"{name} cell {(r, c)} outside {cost.shape} grid")

    n = nr * nc
    passable = np.isfinite(cost)
    idx = np.arange(n).reshape(nr, nc)

    # one edge per unordered neighbour pair: E, S, SE, SW offsets
    rows_l, cols_l, wts_l = [], [], []
    for dr, dc in ((0, 1), (1, 0), (1, 1), (1, -1)):
        r0, r1 = max(0, -dr), min(nr, nr - dr)
        c0, c1 = max(0, -dc), min(nc, nc - dc)
        a = idx[r0:r1, c0:c1]
        b = idx[r0 + dr:r1 + dr, c0 + dc:c1 + dc]
        ok = passable[r0:r1, c0:c1] & passable[r0 + dr:r1 + dr, c0 + dc:c1 + dc]
        step = np.sqrt(float(dr * dr + dc * dc))
        w = 0.5 * (cost[r0:r1, c0:c1] + cost[r0 + dr:r1 + dr, c0 + dc:c1 + dc]) * step
        rows_l.append(a[ok].ravel())
        cols_l.append(b[ok].ravel())
        wts_l.append(w[ok].ravel())
    graph = coo_matrix(
        (np.concatenate(wts_l), (np.concatenate(rows_l), np.concatenate(cols_l))),
        shape=(n, n)).tocsr()

    start = int(start_rc[0]) * nc + int(start_rc[1])
    stop = int(stop_rc[0]) * nc + int(stop_rc[1])
    dist, pred = dijkstra(graph, directed=False, indices=start,
                          return_predecessors=True)
    if not np.isfinite(dist[stop]):
        raise ValueError("no traversable path between start and stop cells")

    path = [stop]
    while path[-1] != start:
        path.append(int(pred[path[-1]]))
    path = np.asarray(path[::-1], dtype=int)
    return np.column_stack([path // nc, path % nc]), float(dist[stop])


def route_through_raster(cost: np.ndarray, gt: GeoTransform,
                         start_xy, stop_xy):
    """LCP between two *map-coordinate* points across a georeferenced raster.

    Ports the notebook ``Est_LCP_Path`` + ``Pixel2Map`` post-processing:
    endpoints are truncated to pixel indices, the path is routed on the grid,
    and the path is returned as pixel-corner map coordinates.

    :returns: ``(path_xy, path_rc, total_cost)``.
    """
    c0, r0 = gt.map_to_pixel(start_xy[0], start_xy[1])
    c1, r1 = gt.map_to_pixel(stop_xy[0], stop_xy[1])
    path_rc, total = route_through_grid(cost, (r0, c0), (r1, c1))
    x, y = gt.pixel_to_map(path_rc[:, 1], path_rc[:, 0])
    return np.column_stack([x, y]), path_rc, total


# =============================================================================
# Trace rasterization -- ports shp_to_raster + Calc_Cost_raster_from_rupture_map
# =============================================================================
def rasterize_traces_xy(traces_xy, pixel_size: float, start_xy, stop_xy,
                        cost_fault: float = COST_FAULT,
                        cost_background: float = COST_BACKGROUND):
    """Burn projected (m) polyline traces into a cost raster.

    Extent handling ports ``shp_to_raster``: the raster covers the union of
    the trace bounding box and the start/end points, buffered by one pixel;
    the geotransform is corner-anchored at (left, top) with negative pixel_dy.
    Burning marks every pixel within half a pixel-diagonal sampling step of a
    segment (dense sampling at pixel/2 spacing - 8-connected coverage like
    GDAL's Bresenham line burning); marked pixels get ``cost_fault``, the rest
    ``cost_background`` (ports ``Calc_Cost_raster_from_rupture_map``).

    :param traces_xy: iterable of (n_i, 2) arrays of projected trace vertices.
    :returns: ``(cost, gt)``.
    """
    allv = np.vstack([np.asarray(t, dtype=float) for t in traces_xy])
    left = min(allv[:, 0].min(), start_xy[0], stop_xy[0]) - pixel_size
    right = max(allv[:, 0].max(), start_xy[0], stop_xy[0]) + pixel_size
    bot = min(allv[:, 1].min(), start_xy[1], stop_xy[1]) - pixel_size
    top = max(allv[:, 1].max(), start_xy[1], stop_xy[1]) + pixel_size
    # int() truncation as in shp_to_raster; +1 guards the degenerate case where
    # truncation would leave the far edge (or a 1-pixel raster) uncovered
    nc = int((right - left) / pixel_size) + 1
    nr = int((top - bot) / pixel_size) + 1
    gt = GeoTransform(left, top, pixel_size, -pixel_size)

    cost = np.full((nr, nc), cost_background, dtype=float)
    step = 0.5 * pixel_size
    for t in traces_xy:
        t = np.asarray(t, dtype=float)
        for i in range(len(t) - 1):
            seg = t[i + 1] - t[i]
            n_s = max(int(np.ceil(np.hypot(*seg) / step)), 1)
            frac = np.linspace(0.0, 1.0, n_s + 1)
            xs = t[i, 0] + frac * seg[0]
            ys = t[i, 1] + frac * seg[1]
            cols = ((xs - gt.origin_x) / gt.pixel_dx).astype(int)
            rows = ((ys - gt.origin_y) / gt.pixel_dy).astype(int)
            ok = (rows >= 0) & (rows < nr) & (cols >= 0) & (cols < nc)
            cost[rows[ok], cols[ok]] = cost_fault
    return cost, gt


# =============================================================================
# Forward-model entry point -- LCP reference line from section traces
# =============================================================================
@dataclass
class LcpResult:
    """Final LCP reference line. Same x/L interface as :class:`ecs.EcsResult`."""
    lon: np.ndarray           # LCP node longitudes
    lat: np.ndarray           # LCP node latitudes
    u: np.ndarray             # cumulative along-path coordinate of nodes (m)
    utm: "_ecs._Utm"
    total_cost: float         # accumulated traversal cost of the path
    pixel_size: float         # cost-raster resolution (m)

    def x_l(self, lon, lat):
        """Normalized along-strike x/L in [0,1] for target points via GC2 on
        the LCP line (identical convention to ``EcsResult.x_l``)."""
        u_km, _t = _ecs.gc2ext_ut(lon, lat, self.lon, self.lat, self.utm)
        u_m = u_km * 1000.0
        umin, umax = self.u.min(), self.u.max()
        L = umax - umin
        if L <= 0:
            return np.zeros(len(np.atleast_1d(lon))), 0.0
        return np.clip((u_m - umin) / L, 0.0, 1.0), L


def _farthest_endpoint_pair(traces_xy):
    """Start/end for the forward case: the farthest-apart pair among the
    section-trace endpoint vertices (no picked points exist in a source
    model). Oriented along the length^2-weighted mean rupture strike so the
    x/L direction matches the ECS ('MRS') convention."""
    ends = np.vstack([np.asarray(t, dtype=float)[[0, -1]] for t in traces_xy])
    d2 = ((ends[:, None, :] - ends[None, :, :]) ** 2).sum(axis=2)
    i, j = np.unravel_index(np.argmax(d2), d2.shape)
    a, b = ends[i], ends[j]

    angs, lens = [], []
    for t in traces_xy:
        t = np.asarray(t, dtype=float)
        angs.append(_ecs.rup_avg_strike_xy(t[:, 0], t[:, 1]))
        lens.append(_ecs.rup_length_xy(t[:, 0], t[:, 1]))
    strike = float(np.average(np.asarray(angs), weights=np.asarray(lens) ** 2))
    sdir = np.array([np.cos(strike), np.sin(strike)])
    if float(np.dot(b - a, sdir)) < 0.0:
        a, b = b, a
    return a, b


def shortcut_path(cost: np.ndarray, path_rc: np.ndarray,
                  samples_per_px: int = 4, max_span: int = 256) -> np.ndarray:
    """Straighten a grid LCP by cost-checked shortcuts ("string pulling").

    The 8-connected grid metric overcharges travel in directions between the
    8 move angles (up to ~8%), so on locally uniform cost the raw path is an
    arbitrary equal-cost staircase (e.g. across an inter-section gap). This
    post-step finds, by dynamic programming over the path nodes, the
    minimum-cost polyline through a *subsequence* of them, where every
    candidate chord is charged the same continuous metric: integrated raster
    cost sampled every 1/``samples_per_px`` pixel (nearest cell) times chord
    length. Uniform-background bridges collapse to the straight chord, while
    chords that would cut across expensive background price themselves out,
    so the path keeps hugging the cheap (fault) corridor at pixel fidelity.
    Charging *all* candidates with one metric (rather than comparing chords
    against the inflated staircase cost) prevents corner-shaving at
    trace/bridge junctions. The result never costs more than the input path
    under the continuous metric.

    :param path_rc: (n, 2) int array of (row, col) grid path nodes.
    :param max_span: DP window - a single shortcut may skip at most this many
        input nodes (bounds the O(n * max_span) chord evaluations).
    :returns: (m, 2) subset of ``path_rc`` nodes (m <= n), endpoints kept.
    """
    n = len(path_rc)
    if n <= 2:
        return path_rc
    rc = path_rc.astype(float)

    def chord_cost(i, j):
        length = float(np.hypot(*(rc[j] - rc[i])))
        n_s = max(int(np.ceil(length * samples_per_px)), 1)
        frac = np.linspace(0.0, 1.0, n_s + 1)
        rr = np.round(rc[i, 0] + frac * (rc[j, 0] - rc[i, 0])).astype(int)
        cc = np.round(rc[i, 1] + frac * (rc[j, 1] - rc[i, 1])).astype(int)
        return float(np.mean(cost[rr, cc])) * length

    best = np.full(n, np.inf)
    pred = np.zeros(n, dtype=int)
    best[0] = 0.0
    for j in range(1, n):
        for i in range(max(0, j - max_span), j):
            c = best[i] + chord_cost(i, j)
            if c < best[j]:
                best[j] = c
                pred[j] = i
    keep = [n - 1]
    while keep[-1] != 0:
        keep.append(int(pred[keep[-1]]))
    return path_rc[keep[::-1]]


def _decimate_collinear(path_rc: np.ndarray) -> np.ndarray:
    """Drop interior nodes of straight pixel runs (exact geometry-preserving;
    keeps GC2 over the path cheap)."""
    if len(path_rc) <= 2:
        return path_rc
    d = np.diff(path_rc, axis=0)
    turn = np.any(d[1:] != d[:-1], axis=1)
    keep = np.concatenate([[True], turn, [True]])
    return path_rc[keep]


def lcp_from_traces(traces, pixel_size: float = PIXEL_SIZE_M,
                    cost_fault: float = COST_FAULT,
                    cost_background: float = COST_BACKGROUND,
                    smooth: bool = False) -> LcpResult:
    """Build an LCP reference line from a multi-section rupture's top-edge traces.

    Forward multi-fault case: the section traces are projected to UTM,
    rasterized into a cost grid (fault pixels cheap, background expensive),
    and the least-cost path is routed between the two farthest-apart section
    endpoints - it follows the traces and bridges inter-section gaps/stepovers
    with straight jumps, yielding one continuous representative line.

    :param traces: iterable of ``(lon, lat)`` arrays, one per section top edge.
    :param pixel_size: cost-raster resolution (m); finer = a more faithful
        path at higher cost (notebook field case uses 20 m).
    :param smooth: optionally apply the cost-checked shortcut pass
        (:func:`shortcut_path`) to the raw grid path. On the uniform-cost
        background between sections the raw 8-connected path is an arbitrary
        member of an equal-cost staircase tie set; the shortcut collapses
        those bridges to the (continuously cheaper) straight chord while
        keeping the on-trace pixels. Default False = the raw
        ``route_through_array``-faithful staircase, exactly like the
        original workflow.
    """
    traces = [(np.asarray(lo, dtype=float), np.asarray(la, dtype=float))
              for lo, la in traces]
    traces = [(lo, la) for lo, la in traces if len(lo) >= 2]
    if len(traces) < 1:
        raise ValueError("lcp_from_traces needs >=1 trace with >=2 vertices")
    all_lon = np.concatenate([lo for lo, _ in traces])
    all_lat = np.concatenate([la for _, la in traces])
    utm = _ecs.utm_for(all_lon, all_lat)
    traces_xy = []
    for lo, la in traces:
        x, y = utm.to_xy(lo, la)
        traces_xy.append(np.column_stack([x, y]))

    start_xy, stop_xy = _farthest_endpoint_pair(traces_xy)
    cost, gt = rasterize_traces_xy(traces_xy, pixel_size, start_xy, stop_xy,
                                   cost_fault, cost_background)
    _, path_rc, total = route_through_raster(cost, gt, start_xy, stop_xy)
    path_rc = _decimate_collinear(path_rc)   # exact; shrinks the DP below
    if smooth:
        path_rc = shortcut_path(cost, path_rc)
    x, y = gt.pixel_to_map(path_rc[:, 1], path_rc[:, 0])
    lon, lat = utm.to_lonlat(x, y)
    u = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(x), np.diff(y)))])
    return LcpResult(lon=np.asarray(lon), lat=np.asarray(lat), u=u, utm=utm,
                     total_cost=total, pixel_size=pixel_size)
