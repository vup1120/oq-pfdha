# -*- coding: utf-8 -*-
"""
FDHA Context classes following OpenQuake ContextMaker pattern.

This module provides unified context arrays for FDHA calculations,
ensuring consistency between hazard curve and hazard map computations.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Dict, Any, Tuple
import logging

logger = logging.getLogger(__name__)


def classify_style(rake: float) -> str:
    """
    Classify faulting style from rake angle (degrees).

    This is the canonical rake→style mapping used by FDHA.

    Returns one of: 'normal', 'reverse', 'strike-slip'.
    """
    rake = float(rake)
    if -150.0 <= rake <= -30.0:
        return 'normal'
    if 30.0 <= rake <= 150.0:
        return 'reverse'
    return 'strike-slip'


@dataclass
class FDHAContext:
    """
    Unified context for FDHA calculations.
    
    All arrays have shape (N,) where N = number of site-rupture pairs.
    This follows the OpenQuake ContextMaker pattern where contexts are
    "flattened" arrays representing all site-rupture combinations.
    
    Attributes:
        sids: Site indices, shape (N,)
        mag: Earthquake magnitude, shape (N,)
        rake: Rake angle in degrees, shape (N,)
        dip: Dip angle in degrees, shape (N,)
        ztor: Depth to top of rupture in km, shape (N,)
        occurrence_rate: Annual occurrence rate, shape (N,)
        vs30: Shear wave velocity in m/s, shape (N,)
        r: Minimum horizontal distance to trace in km, shape (N,)
        rx: Signed perpendicular distance in km (+ = hanging wall), shape (N,)
        x_L: Normalized along-strike position [0,1], shape (N,)
        L: Fault trace length in km, shape (N,)
    """
    # Site identification
    sids: np.ndarray
    
    # Rupture parameters (broadcast to all sites)
    mag: np.ndarray
    rake: np.ndarray
    dip: np.ndarray
    ztor: np.ndarray
    occurrence_rate: np.ndarray
    
    # Site parameters
    vs30: np.ndarray
    
    # Distance parameters (all in km)
    r: np.ndarray
    rx: np.ndarray
    x_L: np.ndarray
    L: np.ndarray
    
    # Site coordinates (for Visini model)
    lons: np.ndarray = field(default=None)
    lats: np.ndarray = field(default=None)

    # Per-reference-line metrics for multi-section ruptures:
    # {method: {'r': (N,), 'x_L': (N,), 'L': (N,)}} for each method some
    # configured model declared via MULTIFAULT_REFERENCE_LINE (the FDHA
    # analogue of hazardlib computing every distance type in the union of
    # the GMPEs' REQUIRES_DISTANCES). None for single-strand ruptures,
    # where every model sees the same trace-based r/x_L/L.
    ref_metrics: Optional[Dict[str, Dict[str, np.ndarray]]] = field(
        default=None, repr=False)

    # Derived properties (computed on demand)
    _style: Optional[np.ndarray] = field(default=None, repr=False)
    _hw_fw: Optional[np.ndarray] = field(default=None, repr=False)
    _near_far_threshold_km: float = field(default=0.2, repr=False)
    
    def __len__(self) -> int:
        """Return number of site-rupture pairs."""
        return len(self.sids)
    
    def __post_init__(self):
        """Convert and validate arrays."""
        # Ensure proper dtypes
        self.sids = np.asarray(self.sids, dtype=np.uint32)
        
        for name in ['mag', 'rake', 'dip', 'ztor', 'occurrence_rate',
                     'vs30', 'r', 'rx', 'x_L', 'L']:
            arr = getattr(self, name)
            if arr is not None:
                setattr(self, name, np.asarray(arr, dtype=np.float64))
        
        # Handle optional coordinate arrays
        if self.lons is not None:
            self.lons = np.asarray(self.lons, dtype=np.float64)
        if self.lats is not None:
            self.lats = np.asarray(self.lats, dtype=np.float64)
        
        # Validate shapes
        N = len(self.sids)
        for name in ['mag', 'rake', 'dip', 'ztor', 'occurrence_rate',
                     'vs30', 'r', 'rx', 'x_L', 'L']:
            arr = getattr(self, name)
            if arr is not None and len(arr) != N:
                raise ValueError(f"{name} has length {len(arr)}, expected {N}")
    
    @property
    def style(self) -> np.ndarray:
        """
        Faulting style derived from rake angle.
        
        Returns:
            Array of strings: 'normal', 'reverse', or 'strike-slip'
        """
        if self._style is None:
            self._style = np.asarray([classify_style(r) for r in self.rake], dtype='U12')
        return self._style
    
    @property
    def hw_fw(self) -> np.ndarray:
        """
        Hanging wall ('HW') or footwall ('FW') indicator.
        
        Based on sign of rx (perpendicular distance).
        Positive rx = hanging wall side.
        """
        if self._hw_fw is None:
            self._hw_fw = np.where(self.rx >= 0, 'HW', 'FW')
        return self._hw_fw
    
    @property
    def near_far(self) -> np.ndarray:
        """
        Near-field ('near') or far-field ('far') label for each site.

        A site is 'near' when its distance to the trace ``r`` is
        <= ``near_far_threshold_km`` (default 0.2 km), otherwise 'far'.
        This regime label is consumed *inside* the secondary (Visini) SR
        computation for the along-strike Monte Carlo (Rank 2). It is NOT the
        principal/distributed split: that decision uses the separate, smaller
        ``r_threshold_km`` (default 0.1 km) via ``get_principal_mask``.
        """
        return np.where(self.r <= self._near_far_threshold_km, 'near', 'far')
    
    @property
    def r_m(self) -> np.ndarray:
        """Distance to trace in meters."""
        return self.r * 1000.0
    
    @property
    def rx_m(self) -> np.ndarray:
        """Signed perpendicular distance in meters."""
        return self.rx * 1000.0
    
    @property
    def L_m(self) -> np.ndarray:
        """Fault length in meters."""
        return self.L * 1000.0
    
    @property
    def site_coords(self) -> Optional[np.ndarray]:
        """
        Site coordinates as (N, 2) array of (lon, lat) pairs.
        
        Returns None if coordinates not available.
        """
        if self.lons is not None and self.lats is not None:
            return np.column_stack([self.lons, self.lats])
        return None
    
    def filter(self, mask: np.ndarray) -> 'FDHAContext':
        """
        Return filtered context where mask is True.
        
        Args:
            mask: Boolean array of shape (N,)
            
        Returns:
            New FDHAContext with only masked elements
        """
        return FDHAContext(
            sids=self.sids[mask],
            mag=self.mag[mask],
            rake=self.rake[mask],
            dip=self.dip[mask],
            ztor=self.ztor[mask],
            occurrence_rate=self.occurrence_rate[mask],
            vs30=self.vs30[mask],
            r=self.r[mask],
            rx=self.rx[mask],
            x_L=self.x_L[mask],
            L=self.L[mask],
            lons=self.lons[mask] if self.lons is not None else None,
            lats=self.lats[mask] if self.lats is not None else None,
            ref_metrics={
                method: {key: arr[mask] for key, arr in metrics.items()}
                for method, metrics in self.ref_metrics.items()
            } if self.ref_metrics is not None else None,
            _near_far_threshold_km=self._near_far_threshold_km,
        )

    def metrics_for(self, method: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Return the (r, x_L, L) arrays for a model's declared reference-line
        method (its ``MULTIFAULT_REFERENCE_LINE`` class attribute).

        For single-strand ruptures — or a method that was not required by
        any configured model — this falls back to the canonical ctx arrays,
        which for single-strand sources are the trace-based values every
        method would produce anyway.
        """
        if self.ref_metrics and method in self.ref_metrics:
            m = self.ref_metrics[method]
            return m['r'], m['x_L'], m['L']
        return self.r, self.x_L, self.L

    def get_principal_mask(self, r_threshold_km: float) -> np.ndarray:
        """
        Get boolean mask for sites within principal rupture zone.
        
        Args:
            r_threshold_km: Distance threshold in km
            
        Returns:
            Boolean array where True = principal zone
        """
        return self.r <= r_threshold_km
    
    def get_distributed_mask(self, r_threshold_km: float) -> np.ndarray:
        """
        Get boolean mask for sites in distributed rupture zone.
        
        Args:
            r_threshold_km: Distance threshold in km
            
        Returns:
            Boolean array where True = distributed zone
        """
        return self.r > r_threshold_km


class FDHAContextMaker:
    """
    Factory for creating FDHA contexts from ruptures and sites.
    
    Follows the OpenQuake ContextMaker pattern with:
    - Distance filtering by maximum_distance
    - Caching of distance calculators by surface geometry
    - Pre-extraction of site arrays for efficiency
    
    Example:
        cmaker = FDHAContextMaker(sitecol, fdha_params, max_distance=50.0)
        
        for rup in source.iter_ruptures():
            ctx = cmaker.get_ctx(rup)
            if ctx is None:
                continue  # All sites too far
            
            # Use ctx for hazard calculation
            P_sr = adapter.compute_primary_sr(ctx, red_cfg)
    """
    
    # Default tolerance for surface rupture detection
    SURFACE_DEPTH_TOLERANCE_KM = 0.5
    
    def __init__(
        self,
        sitecol,
        fdha_params: dict[str, any],
        maximum_distance: float = 50.0
    ):
        """
        Initialize context maker.
        
        Args:
            sitecol: OpenQuake SiteCollection
            fdha_params: Dict of FDHA distance thresholds (km) with keys:
                - 'r_threshold_km' (default 0.1): the principal/distributed
                  split. Sites with r <= threshold are handled by the
                  *primary* (on-trace) SR x FD models; sites beyond it by the
                  *secondary* (distributed) models. This chooses which model
                  family applies and is enforced during hazard-curve
                  integration (see ``get_principal_mask`` /
                  ``get_distributed_mask`` and ``hazard.py``).
                - 'near_far_threshold_km' (default 0.2): the near/far regime
                  split used *inside* the secondary (Visini) computation to
                  label each site 'near' or 'far' for the along-strike Monte
                  Carlo (SR Rank 2). It tunes behaviour within the secondary
                  model rather than selecting the model family (see the
                  ``near_far`` property).
                The two thresholds are independent and have different defaults
                (0.1 vs 0.2 km); do not conflate them.
            maximum_distance: Maximum source-site distance in km
        """
        self.sitecol = sitecol
        self.maximum_distance = maximum_distance
        self.r_threshold_km = fdha_params.get('r_threshold_km', 0.1)
        self.near_far_threshold_km = fdha_params.get('near_far_threshold_km', 0.2)
        # Reference-line methods required by the configured models on
        # multi-section ruptures (union of the models' declared
        # MULTIFAULT_REFERENCE_LINE attributes, collected by the calculator —
        # the FDHA analogue of hazardlib's union of REQUIRES_DISTANCES).
        # 'segments' is always added for multi-section ruptures in get_ctx:
        # the canonical ctx.r (principal/distributed mask, near/far label)
        # is the distance to the nearest surface-reaching section.
        self.multifault_reference_lines = tuple(
            fdha_params.get('multifault_reference_lines', ('lcp',)))
        
        # Distance calculator cache: hash -> calculator
        self._dist_cache: Dict[int, Any] = {}
        self._cache_hits = 0
        self._cache_misses = 0
        
        # Pre-extract site arrays from SiteCollection
        self._n_sites = len(sitecol)
        
        # Site IDs
        if hasattr(sitecol, 'sids'):
            self._sids = np.asarray(sitecol.sids, dtype=np.uint32)
        else:
            self._sids = np.arange(self._n_sites, dtype=np.uint32)
        
        # Coordinates - OpenQuake SiteCollection has .lons and .lats arrays
        if hasattr(sitecol, 'lons'):
            self._lons = np.asarray(sitecol.lons, dtype=np.float64)
            self._lats = np.asarray(sitecol.lats, dtype=np.float64)
        elif hasattr(sitecol, 'mesh'):
            self._lons = np.asarray(sitecol.mesh.lons, dtype=np.float64)
            self._lats = np.asarray(sitecol.mesh.lats, dtype=np.float64)
        else:
            # Fallback: extract from array or sites
            self._lons, self._lats = self._extract_coords_fallback(sitecol)
        
        # Vs30
        if hasattr(sitecol, 'vs30'):
            self._vs30 = np.asarray(sitecol.vs30, dtype=np.float64)
        else:
            self._vs30 = np.full(self._n_sites, 760.0, dtype=np.float64)
        
        logger.debug(f"FDHAContextMaker initialized with {self._n_sites} sites")
    
    def _extract_coords_fallback(self, sitecol) -> Tuple[np.ndarray, np.ndarray]:
        """Extract coordinates using fallback methods."""
        lons = []
        lats = []
        
        for i in range(len(sitecol)):
            site = sitecol[i]
            if hasattr(site, 'location'):
                lons.append(site.location.longitude)
                lats.append(site.location.latitude)
            elif hasattr(site, '__getitem__'):
                # Structured array record
                if 'lon' in site.dtype.names:
                    lons.append(float(site['lon']))
                    lats.append(float(site['lat']))
                else:
                    raise ValueError(f"Cannot extract coordinates from site: {site}")
            else:
                raise ValueError(f"Unknown site type: {type(site)}")
        
        return np.array(lons, dtype=np.float64), np.array(lats, dtype=np.float64)
    
    def _get_surface_hash(self, surface) -> int:
        """
        Create hash key from surface geometry.

        Multi-section surfaces (multiFaultSource ruptures) are keyed by their
        per-section identity — ``suid`` when the engine set one, else the
        section top-trace endpoints — so ruptures made of different section
        combinations never collide and the expensive concatenated
        ``MultiSurface.mesh`` property is not touched. Single surfaces keep
        the mesh first/last-point key.
        """
        subs = getattr(surface, 'surfaces', None)
        if subs:
            parts = []
            for s in subs:
                suid = getattr(s, 'suid', None)
                if suid is not None:
                    parts.append(('suid', suid))
                    continue
                tor = getattr(s, 'tor', None)
                try:
                    coo = np.asarray(
                        tor.coo if tor is not None else s._get_tor(),
                        dtype=float)
                    parts.append((
                        coo.shape[0],
                        round(float(coo[0, 0]), 6), round(float(coo[0, 1]), 6),
                        round(float(coo[-1, 0]), 6), round(float(coo[-1, 1]), 6),
                    ))
                except Exception:
                    parts.append(('id', id(s)))
            return hash(tuple(parts))

        mesh = surface.mesh
        lons = mesh.lons.flatten()
        lats = mesh.lats.flatten()

        if len(lons) == 0:
            return 0

        key_parts = (
            mesh.lons.shape,
            round(float(lons[0]), 6),
            round(float(lons[-1]), 6),
            round(float(lats[0]), 6),
            round(float(lats[-1]), 6),
        )
        return hash(key_parts)

    def _all_sites_far(self, surface) -> bool:
        """
        Cheap bounding-box screen: True when every site is beyond
        ``maximum_distance`` from the surface's bounding box, so ``get_ctx``
        can bail out BEFORE building distance calculators / reference lines
        (which are costly for multi-section ruptures).

        Conservative by construction: the box separation is a lower bound of
        the true site-surface distance, so a rupture is only screened out
        when it is provably too far; False just means "compute properly".
        """
        from openquake.fdha.calc.utils.rupture_distance import _sections_info

        try:
            info = _sections_info(surface)
            if info is not None:
                lons = np.concatenate([lo for lo, _la, _d in info])
                lats = np.concatenate([la for _lo, la, _d in info])
            else:
                mesh = getattr(surface, 'mesh', None)
                if mesh is None:
                    return False
                lons = np.asarray(mesh.lons, dtype=float).ravel()
                lats = np.asarray(mesh.lats, dtype=float).ravel()
        except Exception:
            return False
        ok = np.isfinite(lons) & np.isfinite(lats)
        if not np.any(ok):
            return False
        lons, lats = lons[ok], lats[ok]

        # bounding-box separation in degrees
        dlon = max(0.0, float(self._lons.min()) - float(lons.max()),
                   float(lons.min()) - float(self._lons.max()))
        dlat = max(0.0, float(self._lats.min()) - float(lats.max()),
                   float(lats.min()) - float(self._lats.max()))
        if dlon == 0.0 and dlat == 0.0:
            return False
        # lower-bound km conversion: use the largest |lat| in play so the
        # per-degree longitude length is never overestimated
        max_abs_lat = min(89.0, max(
            float(np.max(np.abs(lats))), float(np.max(np.abs(self._lats)))))
        dx_km = dlon * 111.32 * np.cos(np.radians(max_abs_lat))
        dy_km = dlat * 110.57
        return bool(np.hypot(dx_km, dy_km) > self.maximum_distance)
    
    def _get_distance_calculator(self, surface, reference_line_method: str = "ecs"):
        """
        Get or create cached distance calculator for surface.

        Args:
            surface: Rupture surface
            reference_line_method: reference-line treatment for multi-section
                ruptures ('ecs' | 'lcp' | 'segments'); ignored by
                single-strand surfaces, whose trace is used directly.

        Returns:
            VectorizedRuptureDistanceCalculator instance
        """
        from openquake.fdha.calc.utils.rupture_distance import (
            VectorizedRuptureDistanceCalculator
        )

        # Use vectorized calculator (with caching)
        key = (self._get_surface_hash(surface), reference_line_method)

        if key in self._dist_cache:
            self._cache_hits += 1
            return self._dist_cache[key]

        self._cache_misses += 1
        calc = VectorizedRuptureDistanceCalculator(
            self.sitecol, surface,
            reference_line_method=reference_line_method)
        self._dist_cache[key] = calc
        return calc
    
    def is_surface_rupturing(self, rupture, tolerance_km: float = None) -> bool:
        """
        Check if rupture reaches surface within tolerance.
        
        Args:
            rupture: OpenQuake rupture object
            tolerance_km: Depth tolerance (default: SURFACE_DEPTH_TOLERANCE_KM)
            
        Returns:
            True if rupture is surface-rupturing
        """
        if tolerance_km is None:
            tolerance_km = self.SURFACE_DEPTH_TOLERANCE_KM
        
        depths = rupture.surface.mesh.depths
        if depths is None or depths.size == 0:
            # No depth information - conservatively include
            return True
        
        min_depth = np.nanmin(depths)
        return min_depth <= tolerance_km
    
    def get_ctx(self, rupture,
                investigation_time: Optional[float] = None
                ) -> Optional[FDHAContext]:
        """
        Create FDHA context for a rupture.

        Args:
            rupture: OpenQuake rupture object — parametric (has an
                ``occurrence_rate``) or non-parametric (has a PMF in
                ``probs_occur``, e.g. multiFaultSource ruptures).
            investigation_time: Time span (years) covered by the PMF of a
                non-parametric rupture; its Poisson-equivalent annual rate is
                ``-ln(P(0 events)) / investigation_time``. Ignored for
                parametric ruptures; defaults to 1.0 when needed but absent.

        Returns:
            FDHAContext or None if all sites are beyond maximum_distance
        """
        # Cheap bounding-box screen BEFORE any distance-calculator /
        # reference-line construction (those are costly for multi-section
        # ruptures and would pollute the cache with far ruptures).
        if self._all_sites_far(rupture.surface):
            return None

        from openquake.fdha.calc.utils.rupture_distance import _sections_info
        sections = _sections_info(rupture.surface)
        is_multisection = sections is not None and len(sections) >= 2

        N = self._n_sites
        ref_metrics = None
        if is_multisection:
            # Union of the reference-line treatments declared by the
            # configured models (MULTIFAULT_REFERENCE_LINE — the FDHA
            # analogue of hazardlib computing the union of the GMPEs'
            # REQUIRES_DISTANCES), plus 'segments' for the canonical ctx.r:
            # the principal/distributed mask and the near/far label need the
            # distance to the nearest surface-reaching section, because a
            # site inside an inter-section gap is not on the principal
            # rupture (a smoothed ECS/LCP bridge would give it r ~ 0).
            needed = set(self.multifault_reference_lines) | {'segments'}
            ref_metrics = {}
            for method in sorted(needed):
                dc = self._get_distance_calculator(rupture.surface, method)
                r_m = np.asarray(
                    dc.calculate_site_to_trace_distances(), dtype=np.float64)
                x_L_m, L_m = dc.calculate_x_l_ratios()
                ref_metrics[method] = {
                    'r': r_m,
                    'x_L': np.asarray(x_L_m, dtype=np.float64),
                    'L': np.full(N, L_m, dtype=np.float64),
                }
            r = ref_metrics['segments']['r']
            # Canonical x_L/L: prefer a smoothed representative line when one
            # was required (most models declare 'lcp'); model adapters select
            # their own method via ctx.metrics_for() regardless.
            canonical = next(m for m in ('lcp', 'ecs', 'segments')
                             if m in ref_metrics)
            x_L = ref_metrics[canonical]['x_L']
            L = ref_metrics[canonical]['L']
            rx = rupture.surface.get_rx_distance(self.sitecol)
        else:
            # Single-strand (or single-section) surface: one trace, one set
            # of metrics shared by every model.
            dist_calc = self._get_distance_calculator(rupture.surface)
            r = dist_calc.calculate_site_to_trace_distances()
            x_L, L = dist_calc.calculate_x_l_ratios()
            if dist_calc.trace_is_original:
                # Signed distance from the exact NRML trace: keeps |rx| == r
                # and the HW/FW sign independent of rupture_mesh_spacing.
                # hazardlib's get_rx_distance uses the resampled mesh top
                # edge, which drifts off the true trace at coarse spacing and
                # can flip the side for near-fault sites.
                rx = dist_calc.calculate_signed_site_to_trace_distances()
            else:
                rx = rupture.surface.get_rx_distance(self.sitecol)

        # Early exit if all sites too far
        if np.all(r > self.maximum_distance):
            return None

        # Extract rupture parameters
        mag = rupture.mag
        rake = getattr(rupture, 'rake', 0.0)
        # Prefer the dip declared in the NRML (attached at parse time, see
        # parsing._attach_original_traces): SimpleFaultSurface.get_dip()
        # averages apparent mesh-cell dips, which on a wiggly trace is
        # biased steep and drifts with rupture_mesh_spacing.
        dip = getattr(rupture.surface, 'original_dip', None)
        if dip is None:
            dip = rupture.surface.get_dip()

        # Depth to top of rupture
        ztor = getattr(rupture.surface, 'ztor', None)
        if ztor is None:
            depths = rupture.surface.mesh.depths
            ztor = float(np.nanmin(depths)) if depths is not None and depths.size > 0 else 0.0

        # Occurrence rate: parametric ruptures carry it directly;
        # non-parametric (PMF) ruptures get the Poisson-equivalent
        # annual rate from the probability of zero events.
        occurrence_rate = getattr(rupture, 'occurrence_rate', None)
        if occurrence_rate is None:
            probs = getattr(rupture, 'probs_occur', None)
            if probs is not None and len(probs) > 0:
                t = 1.0 if investigation_time is None else float(investigation_time)
                occurrence_rate = -np.log(float(probs[0])) / t
            else:
                occurrence_rate = np.nan

        # Build context with all sites
        ctx = FDHAContext(
            sids=self._sids.copy(),
            mag=np.full(N, mag, dtype=np.float64),
            rake=np.full(N, rake, dtype=np.float64),
            dip=np.full(N, dip, dtype=np.float64),
            ztor=np.full(N, ztor, dtype=np.float64),
            occurrence_rate=np.full(N, occurrence_rate, dtype=np.float64),
            vs30=self._vs30.copy(),
            r=np.asarray(r, dtype=np.float64),
            rx=np.asarray(rx, dtype=np.float64),
            x_L=np.asarray(x_L, dtype=np.float64),
            L=np.full(N, L, dtype=np.float64) if np.isscalar(L) else np.asarray(L, dtype=np.float64),
            lons=self._lons.copy(),
            lats=self._lats.copy(),
            ref_metrics=ref_metrics,
            _near_far_threshold_km=self.near_far_threshold_km,
        )

        # Filter by maximum distance
        mask = r <= self.maximum_distance
        if not np.all(mask):
            ctx = ctx.filter(mask)

        return ctx
    
    def get_cache_stats(self) -> dict[str, any]:
        """
        Return cache hit/miss statistics.
        
        Returns:
            Dict with 'hits', 'misses', 'hit_rate', 'cached_surfaces'
        """
        total = self._cache_hits + self._cache_misses
        return {
            'hits': self._cache_hits,
            'misses': self._cache_misses,
            'hit_rate': self._cache_hits / total if total > 0 else 0.0,
            'cached_surfaces': len(self._dist_cache),
        }
    
    def clear_cache(self):
        """Clear the distance calculator cache."""
        self._dist_cache.clear()
        self._cache_hits = 0
        self._cache_misses = 0