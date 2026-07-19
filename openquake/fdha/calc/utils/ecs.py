# -*- coding: utf-8 -*-
"""
Event Coordinate System (ECS) representative reference-line construction.

Pure-Python port of the ECS reference implementation
(``ecs_functions.R`` / ``CalcEventCoordinateSystem.R``; Lavrentiadis et al.,
2024) used to build a single representative principal-rupture reference line
for a multi-section / multi-fault rupture, from which the normalized
along-strike coordinate ``x/L`` (and strike-normal ``T``) are derived.

Scope and constraints (see test/fixtures/ecs/README.md and project memory):

* Runtime is **pure Python** (numpy / scipy / pyproj / pandas). It never
  imports R, calls rpy2, shells out to R, or loads a compiled FORTRAN ``.so``.
* GC2 is taken from oq-engine ``MultiLine`` (not reimplemented here).
* R + mgcv is used only as a throwaway dev oracle to emit the validation
  fixtures that this module is checked against.

This file currently implements the deterministic, fixture-independent layers
(data assembly + weighting + projection + curvature + mean-rupture-strike
start). The penalized thin-plate spline fit (mgcv ``s(u)`` / ``family=mvn``)
and the GC2 iteration are added once the R ``lpmatrix``/``penalty_S`` fixtures
are available to validate the basis reconstruction cell-by-cell.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd
import pyproj

from openquake.hazardlib.geo.line import Line
from openquake.hazardlib.geo.multiline import MultiLine

# GC2 end-extension length (m). ecs_functions.R gc2ext_ut default ext_len (L152).
GC2_EXT_LEN = 50000.0

# --- parameters, verbatim from CalcEventCoordinateSystem.R (L37-54) ----------
RANK2CONSIDER = ("Total", "Principal", "Cumulative")
LAMBDA_P = 0.05          # spline smoothness penalty
WT_DISP_MIN = 0.05       # minimum allowed displacement weight
ABS_WT_RUP = 0.05        # absolute rupture-point weight (flag_wt_rup_opt == 2)
R_THRES_MAX = 1e4        # maximum weighting ratio (weight_ecs_data cap)
FIELD_DISP_WT = "recommended_net_preferred_for_analysis_meters"

# fields the assembly needs after lon/lat normalisation
_LONLAT = {"longitude_degrees": "Longitude", "latitude_degrees": "Latitude"}


# =============================================================================
# Projection (pyproj) -- reproduces ecs_functions.R longlat2UTM / UTM2longlat
# =============================================================================
def longlat2utm_zone(lon_mean: float) -> int:
    """UTM zone from mean longitude. Mirrors ``longlat2UTMzone`` (R L54-58)."""
    return int((np.floor((lon_mean + 180.0) / 6.0) % 60) + 1)


@dataclass
class _Utm:
    """A fixed WGS84<->UTM projection, matching R's ``+proj=utm +datum=WGS84``."""
    zone: int
    _fwd: pyproj.Transformer = field(repr=False)
    _inv: pyproj.Transformer = field(repr=False)

    @classmethod
    def for_zone(cls, zone: int) -> "_Utm":
        crs = f"+proj=utm +zone={zone} +datum=WGS84 +units=m +no_defs"
        wgs = "+proj=longlat +datum=WGS84 +no_defs"
        return cls(
            zone,
            pyproj.Transformer.from_crs(wgs, crs, always_xy=True),
            pyproj.Transformer.from_crs(crs, wgs, always_xy=True),
        )

    def to_xy(self, lon, lat):
        x, y = self._fwd.transform(np.asarray(lon), np.asarray(lat))
        return np.asarray(x), np.asarray(y)

    def to_lonlat(self, x, y):
        lon, lat = self._inv.transform(np.asarray(x), np.asarray(y))
        return np.asarray(lon), np.asarray(lat)


def utm_for(lon, lat, zone: Optional[int] = None) -> _Utm:
    """Build the UTM projection for a set of points (zone from mean lon if unset).

    Mirrors ``longlat2UTM`` (R L61-85): the zone is derived from the mean of
    the supplied coordinates unless one is given explicitly.
    """
    if zone is None:
        zone = longlat2utm_zone(float(np.mean(np.asarray(lon))))
    return _Utm.for_zone(zone)


# =============================================================================
# Geometry helpers -- ports of ecs_functions.R rup_length / rup_avg_strike
# =============================================================================
def rup_length_xy(x: np.ndarray, y: np.ndarray) -> float:
    """Polyline length in projected (m) coordinates. Ports ``rup_length`` (R L533)."""
    dx = np.diff(x)
    dy = np.diff(y)
    return float(np.sum(np.sqrt(dx * dx + dy * dy)))


def rup_avg_strike_xy(x: np.ndarray, y: np.ndarray) -> float:
    """Length-weighted mean strike angle (rad, mod pi). Ports ``rup_avg_strike`` (R L506).

    Caller must pass vertices already sorted by NODE_ID.
    """
    dx = np.diff(x)
    dy = np.diff(y)
    seg_len = np.sqrt(dx * dx + dy * dy)
    ang = np.mod(np.arctan2(dy, dx), np.pi)
    return float(np.average(ang, weights=seg_len))


# =============================================================================
# Curvature -- ports compute_2deriv_central / compute_curvature (R L611-647)
# =============================================================================
def _d2_central(u: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Central 2nd derivative d2f/du2 on a non-uniform grid. Ports ``compute_2deriv_central`` (R L611-630).

    R uses ``diff(.,lag=2)`` = ``x[i+2]-x[i]`` (a single lag-2 difference), which
    is ``u[2:]-u[:-2]`` here -- NOT numpy's ``np.diff(u, 2)`` (the iterated 2nd
    difference). The interior length is ``n-2``; R then pads with
    ``c(inner[1], inner, inner[n_pt])`` where ``inner[n_pt]`` is out of bounds
    and evaluates to ``NA``. That endpoint quirk is reproduced faithfully (the
    last curvature value is NaN) and will be confirmed against the R fixture.
    """
    d2f = np.diff(f[1:]) - np.diff(f[:-1])           # f[i+1]-2f[i]+f[i-1]
    du2 = ((u[2:] - u[:-2]) / 2.0) ** 2              # ((u[i+1]-u[i-1])/2)^2
    inner = d2f / du2
    return np.concatenate([[inner[0]], inner, [np.nan]])


def compute_curvature(u: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Curvature magnitude sqrt((x'')^2+(y'')^2). Ports ``compute_curvature`` (R L638)."""
    return np.sqrt(_d2_central(u, x) ** 2 + _d2_central(u, y) ** 2)


# =============================================================================
# GC2 (u,t) via oq-engine MultiLine, with R's 50 km end-extension + u re-anchor
# Ports ecs_functions.R gc2ext_ut (L152-174). GC2 itself is oq-engine's, NOT
# reimplemented here.
# =============================================================================
def _extend_lonlat(ecs_lon, ecs_lat, utm: "_Utm", ext_len=GC2_EXT_LEN):
    """Prepend/append a point ``ext_len`` m beyond each end along the end segment.

    R builds the extension in projected (UTM) x,y (gc2ext_ut L155-162); we do the
    same via pyproj, then map back to lon/lat so the line can be fed to MultiLine.
    """
    x, y = utm.to_xy(ecs_lon, ecs_lat)
    nx, ny = np.diff(x), np.diff(y)
    seg = np.sqrt(nx * nx + ny * ny)
    u0x, u0y = nx[0] / seg[0], ny[0] / seg[0]            # first-segment unit vec
    u1x, u1y = nx[-1] / seg[-1], ny[-1] / seg[-1]        # last-segment unit vec
    x_ext = np.concatenate([[x[0] - ext_len * u0x], x, [x[-1] + ext_len * u1x]])
    y_ext = np.concatenate([[y[0] - ext_len * u0y], y, [y[-1] + ext_len * u1y]])
    lon_ext, lat_ext = utm.to_lonlat(x_ext, y_ext)
    return lon_ext, lat_ext


def gc2ext_ut(query_lon, query_lat, ecs_lon, ecs_lat, utm: "_Utm",
              ext_len=GC2_EXT_LEN):
    """GC2 (u, t) of query points w.r.t. an ECS reference line. Ports ``gc2ext_ut``.

    Reuses oq-engine ``MultiLine`` GC2. The reference line is extended by
    ``ext_len`` m past each end (so off-end points get well-defined ``u``), then
    ``u`` is re-anchored so the first *original* ECS vertex sits at ``u = 0``
    (R subtracts ``ecs_str_ut`` and asserts it == ext_len; ecs_functions.R L168-171).

    :returns: ``(u, t)`` arrays. ``u`` is along-strike (oq ``U``), ``t`` is
        strike-normal (oq ``T``). Sign/origin conventions are validated against
        the R ``fault_disp`` fixture (GC2 gate).

    UNITS: oq-engine ``MultiLine`` GC2 returns ``(T, U)`` in **kilometres**,
    whereas R ``gc2ext_ut`` (UTM) returns **metres**. This function therefore
    returns km; multiply by 1000 to compare against the R ``fault_disp`` u/t
    fixture. The unit cancels in ``x/L = u / u_max`` so it is irrelevant to the
    final ratio, but it matters for the GC2 gate comparison.
    """
    ecs_lon = np.array(ecs_lon, dtype=float)
    ecs_lat = np.array(ecs_lat, dtype=float)
    # MultiLine GC2 (numba) requires WRITABLE contiguous float64 arrays; pandas
    # .to_numpy() can return read-only blocks, so force a writable copy.
    q_lon = np.array(query_lon, dtype=float)
    q_lat = np.array(query_lat, dtype=float)
    lon_ext, lat_ext = _extend_lonlat(ecs_lon, ecs_lat, utm, ext_len)
    ml = MultiLine([Line.from_vectors(lon_ext, lat_ext)])
    t_q, u_q = ml.get_tu(q_lon, q_lat)
    # re-anchor: u of the first original ECS vertex -> 0
    t0, u0 = ml.get_tu(np.array([ecs_lon[0]]), np.array([ecs_lat[0]]))
    return np.asarray(u_q) - float(u0[0]), np.asarray(t_q)


# =============================================================================
# Data assembly + weighting -- ports CalcEventCoordinateSystem.R L120-198
# =============================================================================
def _normalise_lonlat(df: pd.DataFrame) -> pd.DataFrame:
    return df.rename(columns={k: v for k, v in _LONLAT.items() if k in df.columns})


def assemble_data4ecs(
    fault_disp_all: pd.DataFrame,
    fault_rup_all: pd.DataFrame,
    *,
    field_disp_wt: str = FIELD_DISP_WT,
    rank2consider=RANK2CONSIDER,
    abs_wt_rup: float = ABS_WT_RUP,
    wt_disp_min: float = WT_DISP_MIN,
    flag_wt_rup_len: bool = True,
    use_disp: bool = True,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """Build ``data4ecs`` with integer weights, faithfully to the R reference.

    Ports CalcEventCoordinateSystem.R L120-187. Displacement points are weighted
    by ``|field_disp_wt|`` (floored at ``wt_disp_min``); rupture vertices by
    ``abs_wt_rup`` optionally scaled by per-segment length. Both kept only for
    ``rank in rank2consider``. Weights are normalised by the global minimum and
    rounded to integers (the replication counts used by :func:`weight_ecs_data`).

    With ``use_disp=False`` the displacement component is dropped entirely
    (forward multi-fault case: only rupture vertices exist).

    :returns: ``(data4ecs, wt_disp_int, wt_rup_int)`` where ``data4ecs`` has
        columns Longitude, Latitude, wt (and RUP_ID/NODE_ID for rupture rows).
    """
    fault_disp_all = _normalise_lonlat(fault_disp_all)
    fault_rup_all = _normalise_lonlat(fault_rup_all)

    # --- displacement weights (R L128-139) ---
    if use_disp:
        wt_disp = np.abs(fault_disp_all[field_disp_wt].to_numpy(dtype=float))
        keep_d = (~np.isnan(wt_disp)) & fault_disp_all["rank"].isin(rank2consider).to_numpy()
        # outliers (field < -900) become NaN via abs()>900 check in R via i_outliers_d;
        # here the rank/np.isnan gate plus the explicit outlier removal below match R.
        keep_d &= fault_disp_all[field_disp_wt].to_numpy(dtype=float) >= -900.0
        wt_disp = wt_disp[keep_d]
        wt_disp[wt_disp < wt_disp_min] = wt_disp_min
    else:
        keep_d = np.zeros(len(fault_disp_all), dtype=bool)
        wt_disp = np.array([], dtype=float)

    # --- rupture weights (R L140-164) ---
    wt_rup = np.full(len(fault_rup_all), abs_wt_rup, dtype=float)
    if flag_wt_rup_len:
        wt_len = np.full(len(fault_rup_all), np.nan)
        for r_id, idx in fault_rup_all.groupby("RUP_ID").groups.items():
            sub = fault_rup_all.loc[idx].sort_values("NODE_ID")
            seg_len = _rup_seg_length_lonlat(sub)
            wt_len[fault_rup_all.index.get_indexer(idx)] = seg_len / len(idx)
        wt_rup = wt_rup * wt_len
    keep_r = (~np.isnan(wt_rup)) & fault_rup_all["rank"].isin(rank2consider).to_numpy()
    wt_rup = wt_rup[keep_r]

    # --- normalise to integer replication weights (R L165-170) ---
    pool = [wt_rup]
    if wt_disp.size:
        pool.append(wt_disp)
    min_wt = min(float(np.min(a)) for a in pool if a.size)
    wt_disp_int = np.round(wt_disp / min_wt).astype(int) if wt_disp.size else wt_disp
    wt_rup_int = np.round(wt_rup / min_wt).astype(int)

    # --- assemble (R L181-187) ---
    frames = []
    if wt_disp.size:
        d = fault_disp_all.loc[keep_d, ["Longitude", "Latitude"]].copy()
        d["wt"] = wt_disp_int
        frames.append(d)
    r = fault_rup_all.loc[keep_r, ["Longitude", "Latitude", "RUP_ID", "NODE_ID"]].copy()
    r["wt"] = wt_rup_int
    frames.append(r)
    data4ecs = pd.concat(frames, ignore_index=True)
    return data4ecs, wt_disp_int, wt_rup_int


@dataclass
class EcsResult:
    """Final ECS reference line and helpers."""
    lon: np.ndarray           # ECS node longitudes
    lat: np.ndarray           # ECS node latitudes
    u: np.ndarray             # along-strike coordinate of ECS nodes (m)
    curv: np.ndarray          # curvature at ECS nodes
    utm: "_Utm"
    n_iter: int
    flt_ds: float             # final mean offset between iterations (m)

    def x_l(self, lon, lat):
        """Normalized along-strike x/L in [0,1] for target points via GC2 on the ECS."""
        u_km, _t = gc2ext_ut(lon, lat, self.lon, self.lat, self.utm)
        u_m = u_km * 1000.0
        umin, umax = self.u.min(), self.u.max()
        L = umax - umin
        if L <= 0:
            return np.zeros(len(np.atleast_1d(lon))), 0.0
        return np.clip((u_m - umin) / L, 0.0, 1.0), L


def ecs_main(data4ecs_wt: pd.DataFrame, lambda_p: float = LAMBDA_P,
             ecs_du: float = 100.0, flt_max_ds: float = 50.0,
             start: str = "PCA", rup_df: Optional[pd.DataFrame] = None,
             max_iter: int = 50) -> EcsResult:
    """Construct the ECS reference line from weighted points. Ports ``ecs_main``.

    :param data4ecs_wt: replicated weighted points (Longitude, Latitude[, RUP_ID]).
    :param lambda_p: spline smoothness penalty (sp = lambda_p / fault_len).
    :param ecs_du: along-strike node spacing of the ECS (m).
    :param flt_max_ds: convergence threshold on mean point offset between
        iterations (m).
    :param start: 'PCA' (weighted principal direction) or 'MRS' (mean rupture
        strike; requires ``rup_df`` with RUP_ID/NODE_ID).
    """
    lon = data4ecs_wt["Longitude"].to_numpy(float)
    lat = data4ecs_wt["Latitude"].to_numpy(float)
    utm = utm_for(lon, lat)
    x, y = utm.to_xy(lon, lat)
    ox, oy = x.mean(), y.mean()
    x, y = x - ox, y - oy

    if start == "PCA":
        rot = pca_strike(x, y, np.ones_like(x))
    elif start == "MRS":
        rup = _normalise_lonlat(rup_df if rup_df is not None else data4ecs_wt)
        rot = mrs_strike(rup, utm)
    else:
        raise ValueError(f"unknown start solution {start!r}")

    ex, ey = initial_ecs(x, y, rot)                       # centred xy of start ECS
    elon, elat = utm.to_lonlat(ex + ox, ey + oy)

    # initial GC2 (u in metres) of the data points vs the start ECS
    u_km, t_km = gc2ext_ut(lon, lat, elon, elat, utm)
    u, t = u_km * 1000.0, t_km * 1000.0

    n_iter = 0
    flt_ds = np.inf
    for n_iter in range(1, max_iter + 1):
        fault_len = float(u.max() - u.min())
        sp = lambda_p / fault_len
        predict, _beta = fit_spline_xy(u, x, y, sp)
        # ECS along-strike grid (R: round to ecs_du multiples)
        u0 = np.round(u.min() / ecs_du) * ecs_du
        u1 = np.round(u.max() / ecs_du) * ecs_du
        ecs_u = np.arange(u0, u1 + 0.5 * ecs_du, ecs_du)
        ex, ey = predict(ecs_u)                            # new ECS in centred xy
        elon, elat = utm.to_lonlat(ex + ox, ey + oy)
        # recompute GC2 of the data vs the updated ECS
        u_km, t_km = gc2ext_ut(lon, lat, elon, elat, utm)
        un, tn = u_km * 1000.0, t_km * 1000.0
        flt_ds = float(np.mean(np.sqrt((un - u) ** 2 + (tn - t) ** 2)))
        u, t = un, tn
        if flt_ds < flt_max_ds:
            break

    curv = compute_curvature(ecs_u, ex, ey)
    return EcsResult(lon=elon, lat=elat, u=ecs_u, curv=curv,
                     utm=utm, n_iter=n_iter, flt_ds=flt_ds)


def ecs_from_traces(traces, lambda_p: float = LAMBDA_P, ecs_du: float = 100.0,
                    flt_max_ds: float = 50.0) -> EcsResult:
    """Build an ECS reference line from a multi-section rupture's top-edge traces.

    Forward multi-fault case (no displacement data): each section trace becomes
    a set of rupture vertices (rank 'Principal', weighted by ``abs_wt_rup`` and
    segment length), and the ECS is grown with the mean-rupture-strike start.

    :param traces: iterable of ``(lon, lat)`` arrays, one per section top edge.
    """
    rows = []
    for sid, (lon, lat) in enumerate(traces):
        lon = np.asarray(lon, float); lat = np.asarray(lat, float)
        if len(lon) < 2:
            continue                                   # need >=2 vertices for strike/length
        for nid in range(len(lon)):
            rows.append((lon[nid], lat[nid], sid, nid, "Principal"))
    if len(rows) < 2:
        raise ValueError("ecs_from_traces needs >=2 usable section vertices")
    rup = pd.DataFrame(rows, columns=["longitude_degrees", "latitude_degrees",
                                      "RUP_ID", "NODE_ID", "rank"])
    disp = pd.DataFrame(columns=["longitude_degrees", "latitude_degrees",
                                 "rank", FIELD_DISP_WT])
    data4ecs, _, _ = assemble_data4ecs(disp, rup, use_disp=False)
    wt = weight_ecs_data(data4ecs)
    return ecs_main(wt, lambda_p=lambda_p, ecs_du=ecs_du, flt_max_ds=flt_max_ds,
                    start="MRS", rup_df=data4ecs)


def _rup_seg_length_lonlat(sub: pd.DataFrame) -> float:
    """Length (m) of one rupture's vertices, projected to its own UTM zone.

    Mirrors ``rup_length`` (R L533) which projects via ``longlat2UTM`` when x/y
    are absent.
    """
    proj = utm_for(sub["Longitude"].to_numpy(), sub["Latitude"].to_numpy())
    x, y = proj.to_xy(sub["Longitude"].to_numpy(), sub["Latitude"].to_numpy())
    return rup_length_xy(x, y)


# =============================================================================
# Starting solution -- ports ecs_main L300-340 (PCA / mean-rupture-strike)
# =============================================================================
def pca_strike(x: np.ndarray, y: np.ndarray, wt: np.ndarray) -> float:
    """Rotation angle (rad) of the weighted principal direction. Ports ecs_main L302-308.

    R: ``xy_pca <- xy * wt`` (row-scale by weight), centre, ``prcomp``; the angle
    is ``-atan2(v[2,1], v[1,1])`` of the first rotation column.
    """
    xy = np.column_stack([x, y]) * wt[:, None]
    xy = xy - xy.mean(axis=0)
    # principal directions = right singular vectors (prcomp uses SVD of centred X)
    _, _, vt = np.linalg.svd(xy, full_matrices=False)
    v0 = vt[0]                      # first principal axis (x,y)
    return float(-np.arctan2(v0[1], v0[0]))


def mrs_strike(rup_df: pd.DataFrame, utm: "_Utm") -> float:
    """Rotation angle (rad) = -(length^2-weighted mean rupture strike). Ports ecs_main L312-336.

    Each rupture's average strike (``rup_avg_strike``) is averaged with weights
    ``rup_len^2``; the ECS rotation is the negative of that fault strike angle.
    """
    angs, lens = [], []
    for _r, idx in rup_df.groupby("RUP_ID").groups.items():
        sub = rup_df.loc[idx].sort_values("NODE_ID")
        x, y = utm.to_xy(sub["Longitude"].to_numpy(), sub["Latitude"].to_numpy())
        angs.append(rup_avg_strike_xy(x, y))
        lens.append(rup_length_xy(x, y))
    angs, lens = np.asarray(angs), np.asarray(lens)
    fault_strike = float(np.average(angs, weights=lens ** 2))
    return -fault_strike


def initial_ecs(x: np.ndarray, y: np.ndarray, rot_th: float, n_pt: int = 10):
    """Coarse straight starting ECS: ``n_pt`` points along the strike direction.

    Ports ecs_main L317-326. Points are placed at ``t=0`` spanning the rotated
    along-strike extent of the (centred) data, then rotated back to x,y.
    """
    c, s = np.cos(rot_th), np.sin(rot_th)
    rot = np.array([[c, -s], [s, c]])
    ut = (rot @ np.column_stack([x, y]).T).T          # rotate data into (u,t)
    u_lin = np.linspace(ut[:, 0].min(), ut[:, 0].max(), n_pt)
    ecs_ut = np.column_stack([u_lin, np.zeros(n_pt)])
    ecs_xy = (np.linalg.solve(rot, ecs_ut.T)).T       # rotate back to x,y
    return ecs_xy[:, 0], ecs_xy[:, 1]


# =============================================================================
# Penalized thin-plate spline -- reconstruction of mgcv s(u) tp + mvn family.
#
# FIDELITY NOTE (project decision 2026-06-30): this reproduces mgcv's tprs
# (Wood 2003) penalized fit in pure Python. It is faithful to the *method* but
# not bit-identical to mgcv's internal "repara" parameterisation. On the
# Calingiri validation event it reproduces mgcv's fitted ECS trace to ~0.84 m
# max (mean ~0.31 m), with NO spiking at the rupture tips -- the max sits in
# the interior (tips ~0.34 m), so the propagated x/L deviation is <3e-4, far
# below the FDM aleatory sigma. Two deliberate, documented simplifications,
# both validated as negligible:
#   1. x and y are fit by SEPARATE penalized least squares, not the mvn-coupled
#      GLS. The mvn coupling moves the trace ~0.2 m (x/L ~7e-5) but its 2x2
#      residual covariance is near-singular for an ECS (collinear across-strike
#      residuals), making the coupled solve unstable; mgcv regularises it
#      internally. The stable separate fit is used.
#   2. Smoothing strength uses mgcv's documented ``scale.penalty`` normalisation
#      computed in the natural ("whitened") basis (B'B = I), which is the
#      textbook stable penalized-regression parameterisation -- NOT mgcv's
#      private repara. This reproduces mgcv's smoothing *level* (tip-safe)
#      without porting mgcv internals.
# See test_ecs_spline.py for the committed tolerance: a real basis/penalty bug
# fails it; mgcv's internal-repara convention residual passes.
# =============================================================================
TP_K = 10  # mgcv default basis dimension for s(u)


def _tprs(u_knots: np.ndarray, k: int = TP_K):
    """Build a thin-plate regression spline (Wood 2003, d=1, m=2) on ``u_knots``.

    Returns ``(evalfn, S, shift)`` where ``evalfn(u)`` gives the (n, k) basis
    matrix and ``S`` is the (k, k) penalty (rank k-2, the linear null-space term
    unpenalised). The wiggly directions are the top ``k-2`` eigenvectors of the
    ``|r|^3`` kernel restricted to the complement of the polynomial null space
    {1, u} (the identifiability constraint ``T'delta = 0``).
    """
    shift = float(np.mean(u_knots))
    kn = np.unique(u_knots) - shift
    # the basis dimension cannot exceed the number of distinct knots; mgcv
    # likewise reduces k for sparse data. Need >=2 knots for the {1, u} null
    # space; any extra knots add wiggly directions.
    k = int(min(k, len(kn)))
    nwig = max(k - 2, 0)
    Ek = np.abs(kn[:, None] - kn[None, :]) ** 3
    Tk = np.column_stack([np.ones_like(kn), kn])
    Q, _ = np.linalg.qr(Tk, mode="complete")
    Z = Q[:, 2:]                                    # basis of null(T')
    if nwig > 0:
        w, V = np.linalg.eigh(Z.T @ Ek @ Z)
        order = np.argsort(-np.abs(w))[:nwig]
        Uk = Z @ V[:, order]                        # (n_knots, nwig), T'Uk = 0
        Dk = np.abs(w[order])                       # penalty eigenvalues
    else:
        Uk = np.zeros((len(kn), 0))
        Dk = np.zeros(0)

    def evalfn(u):
        uu = np.asarray(u, dtype=float) - shift
        E = np.abs(uu[:, None] - kn[None, :]) ** 3
        return np.column_stack([E @ Uk, np.ones_like(uu), uu])

    S = np.zeros((k, k))
    if nwig > 0:
        S[:nwig, :nwig] = np.diag(Dk)
    return evalfn, S, shift


def fit_spline_xy(u, x, y, sp, k: int = TP_K):
    """Penalized fit of (x, y) ~ s(u). Reconstruction of ecs_main's GAM call.

    Mirrors ``mgcv::gam(list(x~s(u)-1, y~s(u)-1), family=mvn(d=2), sp=sp)``,
    but fits x and y by SEPARATE penalized least squares rather than the
    coupled ``mvn`` GLS. Rationale: for an ECS the x- and y-residuals are nearly
    collinear (they point across-strike), so the 2x2 residual covariance the
    ``mvn`` family estimates is near-singular and the coupled solve is unstable;
    mgcv regularises this internally. The coupling changes the fitted trace by
    only ~0.2 m on the validation event (x/L deviation ~7e-5), far below the
    committed convention tolerance, so the stable separate fit is used. Penalty
    uses mgcv's documented ``scale.penalty`` normalisation. Returns
    ``predict(u)->(x, y)``.
    """
    u = np.asarray(u, float); x = np.asarray(x, float); y = np.asarray(y, float)
    evalfn, S, _shift = _tprs(u, k)
    B = evalfn(u)
    k = B.shape[1]                                   # actual basis dim (reduced for sparse data)
    BtB = B.T @ B
    # natural ("whitened") parameterisation: reparam T so (BT)'(BT)=I and
    # T'ST=diag. This is the standard stable penalized-regression basis (not an
    # mgcv internal), and it makes mgcv's scale.penalty normalisation behave
    # consistently regardless of the raw |r|^3 column scaling -> reproduces
    # mgcv's *smoothing level* (tip-safe), not just the penalty form.
    reg = 1e-10 * np.trace(BtB) / k
    L = np.linalg.cholesky(BtB + reg * np.eye(k))
    Li = np.linalg.inv(L)
    val, Vec = np.linalg.eigh(Li @ S @ Li.T)
    T = Li.T @ Vec
    Bn = B @ T                                       # Bn'Bn = I
    val = np.clip(val, 0.0, None)
    vmax = val.max() if val.size else 0.0
    nz = val[val > 1e-6 * vmax] if vmax > 0 else val[:0]
    # mgcv scale.penalty in the whitened basis (data-computed -> generalises).
    # When there are no wiggly directions (sparse data -> pure linear fit) the
    # penalty is inert.
    if nz.size:
        fac = np.mean(np.abs(Bn.T @ Bn)) / np.mean(np.abs(nz))
        A = Bn.T @ Bn + sp * fac * np.diag(val)
    else:
        A = Bn.T @ Bn                                # pure linear: Bn'Bn ~ I
    beta_n = np.column_stack([np.linalg.solve(A, Bn.T @ x),
                              np.linalg.solve(A, Bn.T @ y)])
    beta = T @ beta_n                                # back to the raw basis

    def predict(u_new):
        xy = evalfn(u_new) @ beta
        return xy[:, 0], xy[:, 1]

    return predict, beta


def weight_ecs_data(data4ecs: pd.DataFrame, ratio_thres: float = R_THRES_MAX) -> pd.DataFrame:
    """Replicate each observation by its integer weight. Ports ``weight_ecs_data`` (R L465).

    The weight array is normalised by its minimum, capped at ``ratio_thres``,
    and each row repeated ``clip(round(wt), 1, ratio_thres)`` times; the
    returned frame has ``wt`` reset to 1.
    """
    wt = data4ecs["wt"].to_numpy(dtype=float)
    wt = wt / wt.min()
    if wt.max() / wt.min() > ratio_thres:
        wt = np.minimum(wt, ratio_thres)
    n_rep = np.clip(np.round(wt), 1, ratio_thres).astype(int)
    out = data4ecs.loc[data4ecs.index.repeat(n_rep)].reset_index(drop=True)
    out["wt"] = 1
    return out
