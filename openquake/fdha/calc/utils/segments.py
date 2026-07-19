# -*- coding: utf-8 -*-
"""
Segmentation-direct rupture distances for multi-section (multiFaultSource)
ruptures - NO smoothed representative reference line.

Third option of the ``reference_line_method`` config switch
(``'ecs' | 'lcp' | 'segments'``). Unlike the ECS (penalized-spline) and LCP
(least-cost-path) builders, which construct one smooth representative trace
and measure both ``r`` and ``x/L`` against it, this mode takes the rupture
segmentation *directly*, in the spirit of Visini et al. (2025), whose
distributed-faulting regressions are calibrated on distances to the actual
(segmented) principal rupture:

* ``r``  = minimum horizontal distance to ANY section top trace. Crucially,
  inter-section gaps/stepovers are NOT bridged: a site in the middle of a
  gap gets its true distance to the nearest rupture segment, whereas a
  smoothed reference line passing through the gap would (wrongly for this
  model family) assign it ``r ~ 0``.
* ``x/L`` and ``L`` come from GC2 over the raw section traces via oq-engine
  ``MultiLine`` (NOT reimplemented) - the standard engine treatment of
  multi-segment ruptures. ``L`` is the along-strike ``u``-span of the section
  vertices; ``x/L`` is the site ``u`` normalised by that span, clipped to
  [0, 1].

With a multiFaultSource, the probability of a rupture jumping a gap
(``P(gap)``) is already carried by the *source model*: each section
combination is a distinct rupture with its own occurrence probability. The
displacement-model side must therefore not smooth across gaps a second time.
"""
from __future__ import annotations

import numpy as np

from openquake.hazardlib.geo.line import Line
from openquake.hazardlib.geo.multiline import MultiLine

# local-frame Earth radius shared with rupture_distance (safe to import at
# module level: rupture_distance only imports this module lazily, inside
# _build_reference_line)
from openquake.fdha.calc.utils.rupture_distance import R_KM


class SegmentsResult:
    """Raw rupture segmentation used directly as the distance reference.

    Shares the ``x_l(lon, lat) -> (xl, L_m)`` interface of ``EcsResult`` /
    ``LcpResult`` (so ``calculate_x_l_ratios`` needs no special casing) and
    additionally exposes ``r_km(lon, lat)`` - the per-section minimum distance
    that replaces the polyline-to-smoothed-trace ``r`` (the presence of
    ``r_km`` is what routes the distance computation, by duck typing, in
    ``VectorizedRuptureDistanceCalculator``).
    """

    def __init__(self, traces):
        """:param traces: iterable of ``(lon, lat)`` arrays, one per section
        top edge (>= 2 vertices each)."""
        self.traces = [(np.array(lo, dtype=float), np.array(la, dtype=float))
                       for lo, la in traces]
        self.traces = [(lo, la) for lo, la in self.traces if len(lo) >= 2]
        if not self.traces:
            raise ValueError("SegmentsResult needs >=1 trace with >=2 vertices")
        # concatenated vertices: used for plot overlays and as the u anchor
        self.lon = np.concatenate([lo for lo, _ in self.traces])
        self.lat = np.concatenate([la for _, la in self.traces])
        # GC2 (T, U) frame over the raw section polylines (oq-engine MultiLine)
        self._ml = MultiLine(
            [Line.from_vectors(lo.copy(), la.copy()) for lo, la in self.traces])
        _t, u_v = self._ml.get_tu(self.lon.copy(), self.lat.copy())
        self._u_min = float(np.min(u_v))
        self._u_max = float(np.max(u_v))

    @property
    def length_km(self) -> float:
        """Total along-strike length L (km): GC2 u-span of the section
        vertices (gaps contribute their along-strike extent, overlapping
        en-echelon strands are not double counted)."""
        return self._u_max - self._u_min

    def x_l(self, lon, lat):
        """Normalized along-strike x/L in [0, 1] via raw-segmentation GC2.

        :returns: ``(xl, L_m)`` with L in metres (same convention as
            ``EcsResult.x_l`` / ``LcpResult.x_l``).
        """
        q_lon = np.array(np.atleast_1d(lon), dtype=float)
        q_lat = np.array(np.atleast_1d(lat), dtype=float)
        L_km = self.length_km
        if L_km <= 0.0:
            return np.zeros(len(q_lon)), 0.0
        _t, u = self._ml.get_tu(q_lon, q_lat)
        xl = np.clip((np.asarray(u) - self._u_min) / L_km, 0.0, 1.0)
        return xl, L_km * 1000.0

    def r_km(self, lon, lat) -> np.ndarray:
        """Minimum horizontal distance (km) from each site to the nearest
        section top trace. Computed per section (gaps are NOT bridged) in a
        local equirectangular frame centred on that section."""
        q_lon = np.array(np.atleast_1d(lon), dtype=float)
        q_lat = np.array(np.atleast_1d(lat), dtype=float)
        best = np.full(len(q_lon), np.inf)
        for t_lon, t_lat in self.traces:
            best = np.minimum(best, _dist_to_polyline_km(
                q_lon, q_lat, t_lon, t_lat))
        return best


def _dist_to_polyline_km(q_lon, q_lat, t_lon, t_lat) -> np.ndarray:
    """Vectorized min distance (km) from sites to one lon/lat polyline,
    in a local equirectangular frame centred on the polyline."""
    lon0 = float(t_lon[0])
    lat0 = float(np.mean(t_lat))
    coslat = np.cos(np.radians(lat0))

    def to_xy(lon, lat):
        dlon = np.radians(np.asarray(lon) - lon0)
        dlon = np.arctan2(np.sin(dlon), np.cos(dlon))     # IDL-safe
        return R_KM * dlon * coslat, R_KM * np.radians(np.asarray(lat) - lat0)

    tx, ty = to_xy(t_lon, t_lat)
    qx, qy = to_xy(q_lon, q_lat)
    # segments (S, 2); sites (N, 2) -> pointwise segment distance (N, S)
    ax, ay = tx[:-1], ty[:-1]
    bx, by = tx[1:], ty[1:]
    dx, dy = bx - ax, by - ay
    den = dx * dx + dy * dy                                # (S,)
    apx = qx[:, None] - ax[None, :]
    apy = qy[:, None] - ay[None, :]
    with np.errstate(invalid="ignore", divide="ignore"):
        t = (apx * dx[None, :] + apy * dy[None, :]) / den[None, :]
    t = np.where(den[None, :] > 0.0, np.clip(t, 0.0, 1.0), 0.0)
    ex = apx - t * dx[None, :]
    ey = apy - t * dy[None, :]
    return np.sqrt((ex * ex + ey * ey).min(axis=1))


def segments_from_traces(traces) -> SegmentsResult:
    """Build the segmentation-direct distance reference from a multi-section
    rupture's top-edge traces (companion to ``ecs.ecs_from_traces`` /
    ``lcp.lcp_from_traces``)."""
    return SegmentsResult(traces)
