# -*- coding: utf-8 -*-
import numpy as np
from typing import List, Tuple, Optional, Any, TYPE_CHECKING
try:
    from openquake.hazardlib.site import Site, SiteCollection  # type: ignore
    from openquake.hazardlib.geo import Point  # type: ignore
except Exception:  # pragma: no cover - optional dependency for distance-only utilities
    Site = None  # type: ignore[assignment]
    SiteCollection = None  # type: ignore[assignment]
    Point = None  # type: ignore[assignment]


# -------- Multi-section (ECS) detection --------
def _section_traces(surface: Any):
    """Return per-section top-edge ``(lon, lat)`` traces for a multi-section
    rupture surface, or ``None`` if the surface is single-strand.

    Multi-section ruptures (e.g. ``multiFaultSource``) arrive as a MultiSurface
    whose ``.surfaces`` each expose a top-of-rupture line (``.tor``). For these
    the single-trace ``_extract_fault_trace_from_mesh`` is invalid (it picks one
    mesh row across disjoint sections), so x/L is computed via the ECS instead.
    """
    subs = getattr(surface, "surfaces", None)
    if not subs or len(subs) <= 1:
        return None
    traces = []
    for s in subs:
        tor = getattr(s, "tor", None)
        try:
            coo = np.asarray(tor.coo if tor is not None else s._get_tor(),
                             dtype=float)
        except Exception:
            return None
        if coo.ndim == 2 and coo.shape[0] >= 2 and coo.shape[1] >= 2:
            traces.append((coo[:, 0], coo[:, 1]))
    return traces if len(traces) >= 2 else None


# -------- Trace extraction from surface.mesh (no get_fault_trace dependency) --------
def _extract_fault_trace_from_mesh(surface: Any) -> np.ndarray:
    """
    Extract the near-surface fault trace (lon, lat) polyline from an OQ surface mesh.
    Prefer a row/column with depths ~ 0; if none, pick the shallowest row by mean depth.
    Returns an ndarray of shape (N, 2).
    """
    # Prefer mesh when available; otherwise fallback to get_fault_trace()
    if not hasattr(surface, "mesh") or surface.mesh is None:
        # Expect get_fault_trace() to return an iterable of (lon, lat[, depth])
        if hasattr(surface, "get_fault_trace"):
            try:
                tr = np.asarray(surface.get_fault_trace(), dtype=float)
                # Ensure shape (N, 2)
                if tr.ndim == 2 and tr.shape[1] >= 2:
                    return tr[:, :2]
            except Exception:
                pass
        # As parity with reference runs, avoid introducing synthetic segments; let caller handle errors
        raise AttributeError("surface has no mesh and get_fault_trace() did not return a valid polyline")
    mesh = surface.mesh
    lons = np.asarray(mesh.lons)
    lats = np.asarray(mesh.lats)
    deps = np.asarray(mesh.depths)

    zero_mask = np.isclose(deps, 0.0)
    if np.any(zero_mask):
        row_counts = zero_mask.sum(axis=1)
        col_counts = zero_mask.sum(axis=0)
        if row_counts.max() >= col_counts.max():
            i = int(np.argmax(row_counts))
            xs = lons[i, :]
            ys = lats[i, :]
        else:
            j = int(np.argmax(col_counts))
            xs = lons[:, j]
            ys = lats[:, j]
    else:
        mean_dep_row = deps.mean(axis=1)
        i = int(np.argmin(mean_dep_row))
        xs = lons[i, :]
        ys = lats[i, :]

    coords = np.column_stack([xs, ys])

    # remove consecutive exact duplicates (bitwise-identical mesh nodes only;
    # do NOT use np.allclose here — its default rtol=1e-5 collapses fine-mesh
    # traces where adjacent nodes differ by < rtol*|lon|, e.g. 0.02 km spacing
    # at lon ~120° produces node differences of ~2e-4° < 1e-5*120 = 1.2e-3°)
    if len(coords) >= 2:
        keep = [0]
        for k in range(1, len(coords)):
            if not np.array_equal(coords[k], coords[k - 1]):
                keep.append(k)
        coords = coords[keep]

    if coords.shape[0] < 2:
        raise ValueError(
            f"Extracted fault trace has only {coords.shape[0]} distinct point(s). "
            f"Cannot compute distances. Check the surface mesh."
        )

    return coords


# ----------------------------- Polyline distance helpers (xy plane) -------------------------
def _min_distance_point_to_polyline_xy(pxy: np.ndarray, poly_xy: np.ndarray) -> float:
    """Shortest distance from point *pxy* to polyline *poly_xy* in the xy plane (km)."""
    d_min = float("inf")
    for i in range(len(poly_xy) - 1):
        p1 = poly_xy[i]
        p2 = poly_xy[i + 1]
        seg = p2 - p1
        seg_len2 = float(np.dot(seg, seg))
        if seg_len2 == 0.0:
            d = float(np.linalg.norm(pxy - p1))
        else:
            t = float(np.clip(np.dot(pxy - p1, seg) / seg_len2, 0.0, 1.0))
            closest = p1 + t * seg
            d = float(np.linalg.norm(pxy - closest))
        if d < d_min:
            d_min = d
    return d_min


def _horizontal_distance_to_trace_km(site_lonlat: np.ndarray, trace_lonlat: np.ndarray) -> float:
    """Horizontal distance (km) from a single site to a trace polyline.

    Uses the same local equirectangular projection as
    ``project_point_onto_trace_km``.
    """
    tr = np.asarray(trace_lonlat, dtype=float)
    if tr.shape[0] < 2:
        return 0.0
    tr_lons = unwrap_longitudes(tr[:, 0])
    tr_lats = tr[:, 1]
    lon0 = float(tr_lons[0])
    lat0 = float(np.mean(tr_lats))
    tx, ty = to_local_equirectangular_km(tr_lons, tr_lats, lon0=lon0, lat0=lat0)
    txy = np.column_stack([tx, ty])
    sx, sy = to_local_equirectangular_km(site_lonlat[0], site_lonlat[1], lon0=lon0, lat0=lat0)
    return _min_distance_point_to_polyline_xy(np.array([float(sx), float(sy)]), txy)


# ----------------------------- Base calculators (point sites) -----------------------------
class RuptureDistanceCalculator:
    """
    Distance utilities for a single-site SiteCollection (first site used).

    Multi-site safety
    -----------------
    The scalar methods ``calculate_site_to_trace_distance`` and
    ``calculate_x_l_ratio`` operate on the first site only (via
    ``_first_site_lonlat``). They are SAFE for multi-site jobs because the
    main calculation path never calls them: ``FDHAContextMaker`` always routes
    through ``VectorizedRuptureDistanceCalculator`` (see
    ``FDHAContextMaker._get_distance_calculator``), which provides the
    ``*_distances`` / ``*_ratios`` plural methods that return per-site arrays.
    This base class is retained only as the shared projection/trace-extraction
    foundation for the vectorized subclass.

    Notes
    -----
    - All along-trace distances and total lengths are computed in kilometers
      using a local equirectangular projection centered at the mean trace
      latitude. This avoids degree/km mixing and ensures consistent units.
    - Crossing the International Date Line (±180°) is handled via longitude
      unwrapping/wrapping in the projection step.
    - Returns a tuple (x_over_L, L_km) where x_over_L ∈ [0,1] and L_km ≥ 0.
    """
    def __init__(self, sitecol, rup_surface):
        self.sitecol = sitecol
        self.rup_surface = rup_surface
        # Multi-section ruptures: build the ECS representative reference line and
        # use it as the trace for both r and x/L. Single-strand falls back to the
        # existing single-trace extraction.
        self._ecs = None
        traces = _section_traces(rup_surface)
        if traces is not None:
            try:
                from openquake.fdha.calc.utils import ecs as _ecs_mod
                self._ecs = _ecs_mod.ecs_from_traces(traces)
                self.trace_points = np.column_stack(
                    [self._ecs.lon, self._ecs.lat])
            except Exception:
                self._ecs = None
        if self._ecs is None:
            self.trace_points = _extract_fault_trace_from_mesh(self.rup_surface)

    def calculate_site_to_trace_distance(self) -> float:
        """Return r (km) = horizontal distance from site to surface trace polyline.

        .. deprecated::
            This scalar (first-site-only) method is NOT used by the production
            hazard pipeline — ``FDHAContextMaker`` always uses
            ``VectorizedRuptureDistanceCalculator.calculate_site_to_trace_distances``
            (plural, per-site array).  This method is retained for unit tests
            only; new code must use the vectorized calculator.

        Unlike OQ ``get_min_distance`` (Rrup to the 3-D mesh), this computes
        the shortest distance in a local equirectangular km frame from the
        site to the 2-D fault trace extracted from the mesh top edge.  This
        decouples the result from ``rupture_mesh_spacing`` and avoids the
        massive ``cdist`` allocation that ``get_min_distance`` requires for
        fine meshes on long faults.
        """
        site_lonlat = np.array(_first_site_lonlat(self.sitecol), dtype=float)
        return _horizontal_distance_to_trace_km(site_lonlat, self.trace_points)

    def calculate_x_l_ratio(self) -> Tuple[float, float]:
        """
        Project the site onto the trace and return (x_over_L, L_km).

        .. deprecated::
            This scalar (first-site-only) method is NOT used by the production
            hazard pipeline — ``FDHAContextMaker`` always uses
            ``VectorizedRuptureDistanceCalculator.calculate_x_l_ratios``
            (plural, per-site array).  This method is retained for unit tests
            only; new code must use the vectorized calculator.

        Distances are computed in kilometers using a local equirectangular
        projection centered at the mean trace latitude. Robust to IDL crossing.
        """
        # Robustly extract first site's lon/lat from SiteCollection
        site_lonlat = np.array(_first_site_lonlat(self.sitecol), dtype=float)
        x_km, L_km = project_point_onto_trace_km(site_lonlat, self.trace_points)
        x_over_L = 0.0 if L_km <= 0.0 else float(np.clip(x_km / L_km, 0.0, 1.0))
        return x_over_L, L_km


class VectorizedRuptureDistanceCalculator(RuptureDistanceCalculator):
    """
    Vectorized variant for multiple sites in a SiteCollection.
    """
    def __init__(
        self,
        sitecol: SiteCollection,
        rup_surface: Any
    ) -> None:
        super().__init__(sitecol, rup_surface)
        self.sites = [site.location for site in sitecol]

    def calculate_site_to_trace_distances(self) -> np.ndarray:
        """Return r (km) = horizontal distance to surface trace for all sites.

        Uses a local equirectangular projection (same frame as
        ``calculate_x_l_ratios``) instead of OQ ``get_min_distance`` (Rrup),
        avoiding the O(n_mesh_nodes × n_sites) ``cdist`` allocation.
        """
        tr = np.asarray(self.trace_points, dtype=float)
        if tr.shape[0] < 2:
            return np.zeros(len(self.sites))

        # Set up projection frame (same as calculate_x_l_ratios)
        tr_lons = unwrap_longitudes(tr[:, 0])
        tr_lats = tr[:, 1]
        lon0 = float(tr_lons[0])
        lat0 = float(np.mean(tr_lats))
        tx, ty = to_local_equirectangular_km(tr_lons, tr_lats, lon0=lon0, lat0=lat0)
        txy = np.column_stack([tx, ty])

        dists = np.empty(len(self.sites))
        for idx, site in enumerate(self.sites):
            sx, sy = to_local_equirectangular_km(
                site.longitude, site.latitude, lon0=lon0, lat0=lat0
            )
            dists[idx] = _min_distance_point_to_polyline_xy(
                np.array([float(sx), float(sy)]), txy
            )
        return dists

    def calculate_x_l_ratios(self) -> Tuple[np.ndarray, float]:
        """
        Return (x_over_L for each site, L_km) using orthogonal projection
        onto the polyline with distances in kilometers.

        Uses a local equirectangular projection centered at the mean trace
        latitude; robust to IDL crossing.

        Multi-section ruptures route through the ECS reference line (GC2 along
        the representative principal-rupture trace); single-strand uses the
        orthogonal-projection path below.
        """
        if self._ecs is not None:
            lons = np.array([s.longitude for s in self.sites], dtype=float)
            lats = np.array([s.latitude for s in self.sites], dtype=float)
            xl, L_m = self._ecs.x_l(lons, lats)
            return np.asarray(xl, dtype=float), L_m / 1000.0  # km, like the trace path

        tr = np.asarray(self.trace_points, dtype=float)
        if tr.shape[0] < 2:
            return np.zeros(len(self.sites)), 0.0

        # Projection frame based on trace
        tr_lons = unwrap_longitudes(tr[:, 0])
        tr_lats = tr[:, 1]
        lon0 = float(tr_lons[0])
        lat0 = float(np.mean(tr_lats))
        tx, ty = to_local_equirectangular_km(tr_lons, tr_lats, lon0=lon0, lat0=lat0)
        txy = np.column_stack([tx, ty])
        segs = txy[1:] - txy[:-1]
        seg_lens = np.linalg.norm(segs, axis=1)
        L_km = float(np.sum(seg_lens))

        def _project_point_xy(pxy: np.ndarray, poly_xy: np.ndarray, total_length: float) -> float:
            x_best = 0.0
            d_min = float("inf")
            cumul = 0.0
            for i in range(len(poly_xy) - 1):
                p1 = poly_xy[i]
                p2 = poly_xy[i + 1]
                seg = p2 - p1
                seg_len2 = float(np.dot(seg, seg))
                seg_len = float(np.sqrt(seg_len2))
                
                if seg_len2 == 0.0:
                    d = float(np.linalg.norm(pxy - p1))
                    if d < d_min:
                        d_min = d
                        x_best = cumul
                    continue
                
                t = float(np.clip(np.dot(pxy - p1, seg) / seg_len2, 0.0, 1.0))
                closest = p1 + t * seg
                d = float(np.linalg.norm(pxy - closest))
                if d < d_min:
                    d_min = d
                    x_best = cumul + t * seg_len
                    # Ensure x_best doesn't exceed total length due to floating point precision
                    x_best = min(x_best, total_length)
                cumul += seg_len
            return x_best

        ratios = []
        for site in self.sites:
            sx, sy = to_local_equirectangular_km(site.longitude, site.latitude, lon0=lon0, lat0=lat0)
            x_proj_km = _project_point_xy(np.array([float(sx), float(sy)]), txy, L_km)
            # x_proj_km is now guaranteed to be <= L_km
            ratio = x_proj_km / L_km if L_km > 0 else 0.0
            ratios.append(ratio)
        
        ratios_arr = np.array(ratios, dtype=float)
        
        # DIAGNOSTIC: Check for values outside [0, 1] before any clipping
        out_of_range = (ratios_arr < 0.0) | (ratios_arr > 1.0)
        if np.any(out_of_range):
            import warnings
            n_invalid = int(np.sum(out_of_range))
            invalid_vals = ratios_arr[out_of_range]
            warnings.warn(
                f"VectorizedRuptureDistanceCalculator.calculate_x_l_ratios: "
                f"{n_invalid} ratio(s) outside [0, 1]: {invalid_vals[:5]}... "
                f"This should not happen - check projection algorithm. "
                f"L_km={L_km:.4f}, trace_points={len(self.trace_points)}",
                RuntimeWarning
            )
        
        # Clip to [0, 1] to match single-site version behavior
        return np.clip(ratios_arr, 0.0, 1.0), L_km




def project_point_onto_trace_deg(site: np.ndarray, trace: np.ndarray) -> tuple:
    """
    DEGREE-SPACE: for internal debug only; do not use in production calculators.

    Project a point onto a polyline; return (x, L) measured in degrees.
    The ratio x/L is dimensionless. This ignores Earth curvature and should
    not be used for production hazard calculations.
    """
    x_best = 0.0
    d_min = float("inf")
    L = 0.0
    cumul = 0.0
    for i in range(len(trace) - 1):
        p1 = trace[i]
        p2 = trace[i + 1]
        seg = p2 - p1
        seg_len = np.linalg.norm(seg)
        L += seg_len
        if seg_len == 0:
            continue
        t = np.clip(np.dot(site - p1, seg) / (seg_len ** 2), 0.0, 1.0)
        closest = p1 + t * seg
        d = np.linalg.norm(site - closest)
        if d < d_min:
            d_min = d
            x_best = cumul + t * seg_len
        cumul += seg_len
    return x_best, L


# ----------------------------- New unified km-based helpers -----------------------------
def _first_site_lonlat(sitecol) -> Tuple[float, float]:
    """Extract (lon, lat) for the first site from a SiteCollection or compatible iterable.

    .. deprecated::
        First-site-only helper used solely by the deprecated scalar methods of
        ``RuptureDistanceCalculator``.  The production hazard pipeline never
        calls this — it uses the per-site arrays from
        ``VectorizedRuptureDistanceCalculator``.
    """
    try:
        site0 = sitecol[0]
        # hazardlib Site has .location with .longitude/.latitude
        if hasattr(site0, "location"):
            return float(site0.location.longitude), float(site0.location.latitude)
        # Fallback: tuple-like or record
        if hasattr(site0, "longitude") and hasattr(site0, "latitude"):
            return float(site0.longitude), float(site0.latitude)
    except Exception:
        pass
    # As a last resort, try to coerce
    site0 = list(sitecol)[0]
    return float(site0.location.longitude), float(site0.location.latitude)
def unwrap_longitudes(lons: np.ndarray) -> np.ndarray:
    """
    Unwrap longitudes to avoid jumps greater than 180°, robust to IDL crossings.

    Parameters
    ----------
    lons : np.ndarray
        Array of longitudes in degrees.

    Returns
    -------
    np.ndarray
        Unwrapped longitudes in degrees with minimized jumps between
        consecutive points.
    """
    lons = np.asarray(lons, dtype=float)
    if lons.size == 0:
        return lons
    r = np.deg2rad(lons)
    r_unw = np.unwrap(r, discont=np.deg2rad(180.0))
    return np.rad2deg(r_unw)


def to_local_equirectangular_km(lon, lat, lon0: float, lat0: float):
    """
    Project geographic coordinates to a local equirectangular frame in km.

    Parameters
    ----------
    lon, lat : array-like or float
        Longitudes and latitudes in degrees.
    lon0, lat0 : float
        Projection center in degrees. Typically lon0 is the first trace vertex
        and lat0 is the mean trace latitude.

    Returns
    -------
    (x_km, y_km) : tuple[np.ndarray, np.ndarray] or (float, float)
        Coordinates in kilometers in the local frame.
    """
    R_KM = 6371.0088
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    lon0 = float(lon0)
    lat0 = float(lat0)

    lon_rad = np.deg2rad(lon)
    lat_rad = np.deg2rad(lat)
    lon0_rad = np.deg2rad(lon0)
    lat0_rad = np.deg2rad(lat0)

    # Wrap delta-lon to [-pi, pi] to handle IDL crossing consistently
    dlon = lon_rad - lon0_rad
    dlon = np.arctan2(np.sin(dlon), np.cos(dlon))
    x_km = R_KM * dlon * np.cos(lat0_rad)
    y_km = R_KM * (lat_rad - lat0_rad)
    return x_km, y_km


def project_point_onto_trace_km(site_lonlat: np.ndarray, trace_lonlat: np.ndarray) -> tuple:
    """
    Project a point onto a polyline; return (x_km, L_km), both in kilometers.

    Uses a local equirectangular projection centered at the mean trace
    latitude and first-trace-vertex longitude. Robust to IDL crossing.

    Returns (0.0, 0.0) if the trace has fewer than 2 points.
    """
    tr = np.asarray(trace_lonlat, dtype=float)
    if tr.shape[0] < 2:
        return 0.0, 0.0

    tr_lons = unwrap_longitudes(tr[:, 0])
    tr_lats = tr[:, 1]
    lon0 = float(tr_lons[0])
    lat0 = float(np.mean(tr_lats))

    # Trace to local km
    tx, ty = to_local_equirectangular_km(tr_lons, tr_lats, lon0=lon0, lat0=lat0)
    txy = np.column_stack([tx, ty])
    segs = txy[1:] - txy[:-1]
    seg_lens = np.linalg.norm(segs, axis=1)
    L_km = float(np.sum(seg_lens))
    if L_km <= 0.0:
        return 0.0, 0.0

    # Site to local km
    sx, sy = to_local_equirectangular_km(site_lonlat[0], site_lonlat[1], lon0=lon0, lat0=lat0)
    pxy = np.array([float(sx), float(sy)])

    # Project
    x_best = 0.0
    d_min = float("inf")
    cumul = 0.0
    for i in range(len(txy) - 1):
        p1 = txy[i]
        p2 = txy[i + 1]
        seg = p2 - p1
        seg_len2 = float(np.dot(seg, seg))
        if seg_len2 == 0.0:
            d = float(np.linalg.norm(pxy - p1))
            if d < d_min:
                d_min = d
                x_best = cumul
            continue
        t = float(np.clip(np.dot(pxy - p1, seg) / seg_len2, 0.0, 1.0))
        closest = p1 + t * seg
        d = float(np.linalg.norm(pxy - closest))
        if d < d_min:
            d_min = d
            x_best = cumul + t * float(np.sqrt(seg_len2))
        cumul += float(np.sqrt(seg_len2))

    return x_best, L_km


#def trace_total_length(trace: np.ndarray) -> float:
#    return float(sum(np.linalg.norm(trace[i + 1] - trace[i]) for i in range(len(trace) - 1)))

def trace_total_length(trace: np.ndarray) -> float:
    """
    Compute the total length of a fault surface trace in kilometers.

    Parameters
    ----------
    trace : np.ndarray, shape (N, 2) or (N, >=2)
        Each row is [lon, lat] in degrees. Any extra columns (e.g., depth) are ignored.

    Returns
    -------
    float
        Total length in kilometers.
    """
    pts = np.asarray(trace, dtype=float)
    if pts.ndim != 2 or pts.shape[0] < 2:
        return 0.0

    # Use only lon/lat columns
    lon = np.deg2rad(pts[:, 0])
    lat = np.deg2rad(pts[:, 1])

    dlon = np.diff(lon)
    dlat = np.diff(lat)

    # Haversine great-circle distance for each segment (in radians)
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat[:-1]) * np.cos(lat[1:]) * np.sin(dlon / 2.0) ** 2
    a = np.clip(a, 0.0, 1.0)  # guard against floating-point drift
    c = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))

    R_KM = 6371.0088  # mean Earth radius in km
    segment_km = R_KM * c

    return float(np.sum(segment_km))


def cumulative_length_upto_segment_deg(trace: np.ndarray, seg_index: int) -> float:
    """Cumulative length from the start up to (but not including) segment 'seg_index'."""
    if seg_index <= 0:
        return 0.0
    return float(sum(np.linalg.norm(trace[i + 1] - trace[i]) for i in range(seg_index)))