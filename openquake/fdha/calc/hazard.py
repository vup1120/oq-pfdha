# -*- coding: utf-8 -*-
"""
Unified FDHA hazard calculation module.

This module provides a single implementation for both hazard curve and
hazard map calculations.
"""

import numpy as np
from tqdm import tqdm
import logging
from typing import Dict, Any, Optional, Tuple, TYPE_CHECKING

from openquake.fdha.calc.contexts import FDHAContext, FDHAContextMaker
from openquake.hazardlib.site import SiteCollection
from openquake.fdha.calc.calculators import BaseFaultRuptureCalculator

logger = logging.getLogger(__name__)


class ApplicabilityTracker:
    """Track sites evaluated outside a distributed FD model's declared
    applicability range and emit ONE ``logging.warning`` per model per run.

    Distributed displacement models declare their calibrated distance range
    as the ``APPLICABILITY_RANGE`` class attribute (see
    ``primary_surf_displ.base.BaseSecondarySurfDispl``), expressed in the
    model's OWN distance metric (``MULTIFAULT_REFERENCE_LINE``). The tracker
    accumulates, across all ruptures of a run, the ids of sites whose
    distributed contribution was computed at a distance outside that range
    (i.e. an extrapolation of the regression), then reports the offending
    site count once at the end of the run.

    Sites where the distributed term carries zero weight are NOT counted:
    on the sigma = 0 (complementary) W_p path the distributed component is
    masked inside ``|r| <= r_threshold_km``, so e.g. an on-trace site below
    Visini's 5 m data floor is not an extrapolation -- the model is never
    used there.
    """

    def __init__(self, r_threshold_km: float, r_sigma_km: float):
        self._r_threshold_km = float(r_threshold_km)
        self._r_sigma_km = float(r_sigma_km)
        self._offending: Dict[str, set] = {}
        self._sources: Dict[str, str] = {}

    def observe(self, model: Any, ctx: 'FDHAContext') -> None:
        """Record ``ctx`` sites outside ``model``'s declared range."""
        if model is None:
            return
        rng = getattr(model, 'APPLICABILITY_RANGE', None)
        if not rng:
            return
        method = getattr(model, 'MULTIFAULT_REFERENCE_LINE', 'lcp')
        r_sel, _x_L, _L = ctx.metrics_for(method)
        r = np.abs(np.asarray(r_sel, dtype=np.float64))

        outside = np.zeros(r.shape, dtype=bool)
        if 'r_min_km' in rng:
            outside |= r < float(rng['r_min_km'])
        if 'r_max_km' in rng:
            outside |= r > float(rng['r_max_km'])
        if 'r_max_hw_km' in rng or 'r_max_fw_km' in rng:
            # Tool-wide wall convention (cf. Visini2025SecondaryFD):
            # rx < 0 = footwall, rx >= 0 = hanging wall.
            fw = np.asarray(ctx.rx, dtype=np.float64) < 0.0
            if 'r_max_hw_km' in rng:
                outside |= (~fw) & (r > float(rng['r_max_hw_km']))
            if 'r_max_fw_km' in rng:
                outside |= fw & (r > float(rng['r_max_fw_km']))

        # Only count sites where the distributed term actually contributes:
        # sigma = 0 -> complementary split masks distributed inside the
        # boxcar (G = 1 - W_p = 0 there); sigma > 0 -> additive, G = 1
        # everywhere (docs/design/rupture_location_uncertainty.md, D1).
        if self._r_sigma_km == 0.0:
            outside &= np.abs(np.asarray(ctx.r, dtype=np.float64)) \
                > self._r_threshold_km

        if not outside.any():
            return
        name = model.__class__.__name__
        sids = np.asarray(ctx.sids)[outside]
        self._offending.setdefault(name, set()).update(
            int(s) for s in sids)
        self._sources.setdefault(name, str(rng.get('source', '')))

    def emit(self) -> None:
        """Emit the once-per-model warnings (call after the rupture loop)."""
        for name in sorted(self._offending):
            n = len(self._offending[name])
            src = self._sources.get(name, '')
            logger.warning(
                "%s: %d site(s) were evaluated outside the model's declared "
                "applicability range (%s). The distributed regression is "
                "extrapolating beyond its calibration data there; results "
                "are still computed unchanged.",
                name, n, src or 'declared range',
            )


def calculate_fdha_hazard(
    calculator: 'BaseFaultRuptureCalculator',
    sitecol: Optional[SiteCollection] = None,
    show_progress: bool = True
) -> Dict[str, Any]:
    """
    Unified FDHA hazard calculation for both curves and maps.
    
    This function replaces the separate hazard_curve.py and hazard_map_calculator.py
    implementations with a single, consistent calculation.
    
    Args:
        calculator: Configured FDHA calculator with models and parameters
        sitecol: Optional SiteCollection override (for hazard maps with grid)
        show_progress: Show tqdm progress bar
        
    Returns:
        Dictionary of numpy arrays:
        - 'imls': displacement levels (m), shape (n_displ,)
        - 'poes': ANNUAL EXCEEDANCE RATES (historical key name, not
          probabilities), shape (n_sites, n_displ)
        - 'rate_principal': principal contribution, shape (n_sites, n_displ)
        - 'rate_distributed': distributed contribution, shape (n_sites, n_displ)
        - 'site_lons', 'site_lats': site coordinates, shape (n_sites,)
        - 'n_sites', 'n_displ': ints
    """
    # Use provided sitecol or calculator's sitecol
    if sitecol is None:
        sitecol = calculator.sitecol
    
    n_sites = len(sitecol)
    
    target_displacements = calculator.target_displacements
    n_displ = len(target_displacements)
    
    logger.info("Starting FDHA hazard calculation: %d sites, %d displacement "
                "levels", n_sites, n_displ)

    # Get max_distance from config (check multiple sections)
    max_dist = 50.0  # default
    for section in ['calculation', 'parameters', 'erf']:
        section_cfg = calculator.config.get(section, {})
        if 'max_distance_km' in section_cfg:
            max_dist = float(section_cfg['max_distance_km'])
            break
    
    # Create context maker with caching
    cmaker = FDHAContextMaker(
        sitecol=sitecol,
        fdha_params=calculator.get_fdha_params(),
        maximum_distance=max_dist,
    )
    
    # Check for Visini models
    use_visini, visini_calc = _setup_visini_calculator(calculator)
    
    # Initialize rate accumulators
    rate_principal = np.zeros((n_sites, n_displ), dtype=np.float64)
    rate_distributed = np.zeros((n_sites, n_displ), dtype=np.float64)

    # Get model adapters
    adapters = calculator.adapters

    logger.debug("Available adapters: %s", list(adapters.keys()))
    if not adapters:
        logger.warning("No adapters available - models may not have loaded correctly")
    
    # Reduction configs
    p_sr_red_cfg = calculator.p_sr_red_cfg
    s_sr_red_cfg = calculator.s_sr_red_cfg
    r_threshold_km = calculator.r_threshold_km

    # Rupture-location uncertainty (W_p).
    # See docs/design/rupture_location_uncertainty.md, section 2.
    r_sigma_km = calculator.r_sigma_km

    # Applicability advisory (C4): distributed FD models declare their
    # calibrated distance range; sites evaluated beyond it are collected
    # across the whole run and reported ONCE per model after the loop.
    # No behaviour change -- warning only.
    applicability_tracker = ApplicabilityTracker(
        r_threshold_km=r_threshold_km, r_sigma_km=r_sigma_km)
    secondary_fd_model = calculator.secondary_surf_displ_model

    # PMF time span for non-parametric (multiFaultSource) ruptures:
    # get_ctx converts their probs_occur into a Poisson-equivalent annual
    # rate over this investigation time. Parametric ruptures ignore it.
    # The authoritative value is per source (parsed from the NRML header);
    # the INI value is a fallback and must not silently contradict it.
    ini_investigation_time = calculator.config.get(
        'calculation', {}).get('investigation_time')

    # Process each fault source
    fault_sources = calculator.fault_sources.values()
    iterator = tqdm(fault_sources, desc="Computing hazard") \
        if show_progress else fault_sources

    for src in iterator:
        logger.info("Processing fault: %s", src.name)
        investigation_time = _resolve_investigation_time(
            src, ini_investigation_time)

        n_ruptures = 0
        n_surface_rupturing = 0

        for rup in src.iter_ruptures():
            n_ruptures += 1

            # Skip buried ruptures
            if not cmaker.is_surface_rupturing(rup):
                continue
            n_surface_rupturing += 1
            
            # Create context (returns None if all sites too far)
            ctx = cmaker.get_ctx(rup, investigation_time=investigation_time)
            if ctx is None:
                continue

            # Advisory extrapolation bookkeeping (once-per-run warning).
            applicability_tracker.observe(secondary_fd_model, ctx)

            # Calculate hazard contribution for this rupture
            principal_contrib, distributed_contrib = _compute_rupture_contribution(
                ctx=ctx,
                adapters=adapters,
                target_displacements=target_displacements,
                p_sr_red_cfg=p_sr_red_cfg,
                s_sr_red_cfg=s_sr_red_cfg,
                r_threshold_km=r_threshold_km,
                use_visini=use_visini,
                visini_calc=visini_calc,
                calculator=calculator,
                r_sigma_km=r_sigma_km,
            )
            
            # Accumulate by site ID using vectorized operations. The context
            # invariant guarantees one contribution row per ctx site and
            # sids within the sitecol; a violation means a broken
            # adapter/kernel/cmaker, not a condition to paper over.
            sids = ctx.sids
            if len(sids) != principal_contrib.shape[0]:
                raise AssertionError(
                    f"context/contribution shape mismatch: {len(sids)} sids "
                    f"vs {principal_contrib.shape[0]} contribution rows")
            if len(sids) == 0:
                continue
            if int(sids.max()) >= n_sites or int(sids.min()) < 0:
                raise AssertionError(
                    f"context sids outside the site collection: range "
                    f"[{sids.min()}, {sids.max()}] vs {n_sites} sites")

            np.add.at(rate_principal, sids, principal_contrib)
            np.add.at(rate_distributed, sids, distributed_contrib)

        logger.debug("  %s: %d ruptures, %d surface-rupturing",
                     src.name, n_ruptures, n_surface_rupturing)

    # Applicability advisory: one warning per model per run (C4).
    applicability_tracker.emit()

    # Log cache statistics
    stats = cmaker.get_cache_stats()
    logger.info("Distance cache: %d hits, %d misses, %.1f%% hit rate",
                stats['hits'], stats['misses'], 100.0 * stats['hit_rate'])

    return {
        'imls': np.asarray(target_displacements, dtype=np.float64),
        # 'poes' holds ANNUAL EXCEEDANCE RATES; the key name is historical.
        'poes': rate_principal + rate_distributed,
        'rate_principal': rate_principal,
        'rate_distributed': rate_distributed,
        'site_lons': np.asarray(cmaker._lons, dtype=np.float64),
        'site_lats': np.asarray(cmaker._lats, dtype=np.float64),
        'n_sites': n_sites,
        'n_displ': n_displ,
    }


def _resolve_investigation_time(src, ini_time) -> float:
    """Return the PMF time span (years) used to convert a non-parametric
    source's ``probs_occur`` into Poisson-equivalent annual rates.

    The span declared by the source itself (the ``investigation_time``
    attribute parsed from the NRML ``<sourceModel>``/``<geometryModel>``
    header, e.g. on a multiFaultSource) is authoritative: the PMF is
    *defined* over that span. The INI ``[calculation].investigation_time``
    may restate it, but a conflicting INI value would silently rescale every
    non-parametric rate (an XML span of 50 yr read with the 1-yr INI default
    inflates all rates 50x), so a mismatch is rejected loudly.

    Parametric sources carry no such attribute; they fall back to the INI
    value, then 1.0 - their ruptures have explicit annual rates and ignore
    this value anyway.
    """
    src_time = getattr(src, 'investigation_time', None)
    if src_time is None:
        return float(ini_time) if ini_time is not None else 1.0
    src_time = float(src_time)
    if ini_time is not None and abs(float(ini_time) - src_time) > \
            1e-9 * max(abs(src_time), 1.0):
        raise ValueError(
            f"Source '{getattr(src, 'source_id', src)}' declares "
            f"investigation_time={src_time} in its NRML header, but the job "
            f"INI sets [calculation].investigation_time={ini_time}. The PMF "
            "of a non-parametric source is defined over the NRML time span; "
            "remove the INI key or set it to the same value."
        )
    return src_time


def _setup_visini_calculator(
    calculator: 'BaseFaultRuptureCalculator'
) -> Tuple[bool, Optional[Any]]:
    """
    Check for Visini models and initialize calculator if needed.
    
    Args:
        calculator: FDHA calculator instance
        
    Returns:
        Tuple of (use_visini: bool, visini_calc: Optional[VisiniSecondaryCalculator])
    """
    sr_model = calculator.secondary_surf_rup_model
    fd_model = calculator.secondary_surf_displ_model

    # A model declares it needs the combined Visini pipeline via the
    # SECONDARY_PIPELINE class attribute (base default 'generic'); the kernel
    # never matches class names, so a Visini subclass or renamed variant keeps
    # the correct routing (see BaseSecondarySurfRup.SECONDARY_PIPELINE).
    use_visini = (
        getattr(sr_model, 'SECONDARY_PIPELINE', 'generic') == 'visini'
        or getattr(fd_model, 'SECONDARY_PIPELINE', 'generic') == 'visini'
    )

    if not use_visini:
        return False, None

    sr_name = sr_model.__class__.__name__ if sr_model else ''
    fd_name = fd_model.__class__.__name__ if fd_model else ''
    logger.info("Using Visini model: SR=%s, FD=%s", sr_name, fd_name)
    
    from openquake.fdha.calc.visini import VisiniSecondaryCalculator
    from openquake.fdha.calc.rank1p5_loader import attach_rank1p5_surfaces
    
    # Attach rank1p5 surfaces from XML file if available
    attach_rank1p5_surfaces(calculator, getattr(calculator, 'config_path', None))
    
    sec_rup_params = calculator.get_model_parameters('secondary_surf_rup')
    sec_displ_params = calculator.get_model_parameters('secondary_surf_displ')
    
    # Get rupture_traces and rank1p5_traces from model parameters and config
    rupture_traces = sec_rup_params.get('rupture_traces', [])
    
    # Try to get rank1p5_traces from config (TOML format)
    rank1p5_traces = calculator.config.get('rank1p5_ruptures', {}).get('trace', [])
    
    # If empty, try to get from loaded XML surfaces (INI format with rank1p5_traces_file)
    if not rank1p5_traces and hasattr(calculator, 'rank1p5_surface_by_name'):
        # Convert TraceOnlySurfaceAdapter objects to trace dict format
        for name, surface in calculator.rank1p5_surface_by_name.items():
            coords = [[p.longitude, p.latitude] for p in surface.coords_ll] if hasattr(surface, 'coords_ll') else []
            if not coords and hasattr(surface, '_line'):
                coords = [[p.longitude, p.latitude] for p in surface._line.points]
            rank1p5_traces.append({
                'name': name,
                'geometry': {'type': 'Line', 'coords': coords}
            })
    
    
    visini_calc = VisiniSecondaryCalculator(
        base_sec_rup_params={k: v for k, v in sec_rup_params.items() if k not in ('combination', 'case')},
        base_sec_displ_params={k: v for k, v in sec_displ_params.items() if k not in ('combination', 'case')},
        case_label=calculator.case_label,
        pixel_size=sec_rup_params.get('pixel_size', 100),
        along_strike_width=sec_rup_params.get('along_strike_width', None),
        near_far_threshold_km=calculator.near_far_threshold_km,
        rupture_traces=rupture_traces,
        rank1p5_traces=rank1p5_traces,
        segment_sampling=sec_rup_params.get('segment_sampling', 'truncated'),
        distribution_type=sec_rup_params.get('distribution_type', 'uniform'),
    )
    
    return True, visini_calc


def _compute_rupture_contribution(
    ctx: FDHAContext,
    adapters: Dict[str, Any],
    target_displacements: np.ndarray,
    p_sr_red_cfg: Dict[str, Any],
    s_sr_red_cfg: Dict[str, Any],
    r_threshold_km: float,
    use_visini: bool,
    visini_calc: Optional[Any],
    calculator: 'BaseFaultRuptureCalculator',
    r_sigma_km: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute principal and distributed hazard contributions for a rupture.

    Args:
        ctx: FDHA context for this rupture
        adapters: Model adapters dict
        target_displacements: Displacement levels array
        p_sr_red_cfg: Primary SR reduction config
        s_sr_red_cfg: Secondary SR reduction config
        r_threshold_km: W_p boxcar half-width h (sigma == 0 path only)
        use_visini: Whether to use Visini model
        visini_calc: Visini calculator instance (if use_visini)
        calculator: Parent calculator for model access
        r_sigma_km: Two-sided mapping-accuracy sigma for W_p. 0 selects the
            boxcar path; > 0 selects Petersen's pure Gaussian path (pinned,
            fixed +-2 sigma truncation; r_threshold_km plays no role there).

    The combination follows the W_p path: at sigma = 0 the historical
    COMPLEMENTARY boxcar split (inside h principal only, outside distributed
    only - Youngs 2003 / Takao 2013 either/or); at sigma > 0 principal and
    distributed are independent contributions of the same surface-rupturing
    event and are SUMMED (Petersen et al. 2011 eq. 1 + eq. 2; Fig. 10a
    "total hazard").

    Returns:
        Tuple of (principal_contrib, distributed_contrib) arrays, each shape (N_ctx, n_displ)
    """
    from openquake.fdha.calc.location_weight import location_weight
    N_ctx = len(ctx)
    n_displ = len(target_displacements)
    rate = ctx.occurrence_rate[0]
    
    # Initialize output arrays
    principal_contrib = np.zeros((N_ctx, n_displ), dtype=np.float64)
    distributed_contrib = np.zeros((N_ctx, n_displ), dtype=np.float64)
    
    # =========================================================================
    # PRIMARY SURFACE RUPTURE PROBABILITY
    # =========================================================================
    if 'primary_sr' in adapters:
        # Model errors propagate: a crashing model fails the job instead of
        # silently zeroing this rupture's hazard.
        P_sr = adapters['primary_sr'].compute_primary_sr(ctx, p_sr_red_cfg)
    else:
        # Default: assume surface rupture always occurs
        P_sr = np.ones(N_ctx, dtype=np.float64)

    # =========================================================================
    # PRIMARY FAULT DISPLACEMENT PROBABILITY
    # =========================================================================
    if 'primary_fd' in adapters:
        P_fd_primary = adapters['primary_fd'].compute_primary_fd(
            ctx, target_displacements, p_sr_red_cfg
        )
    else:
        P_fd_primary = np.zeros((N_ctx, n_displ), dtype=np.float64)
    
    # =========================================================================
    # AGGREGATE-DEFINITION PRIMARY MODEL: single-bucket path
    # =========================================================================
    # An aggregate-definition model (Sarmiento et al. 2025 Table 1: total
    # displacement across principal AND distributed ruptures in the
    # measurement aperture, e.g. Kuehn2024PrimaryFD or
    # Lavrentiadis2023PrimaryFD_aggregate -- the contract is STATIC, the class choice
    # IS the definition) already contains the distributed contribution, so
    # the split above does not apply:
    #
    #     lambda_total = rate * P_sr * P_fd_aggregate * W_p(r)
    #
    # -- ONE bucket, NO distributed term, W_p per the same two-path kernel
    # as below (docs/design/rupture_location_uncertainty.md, D8). Adding a
    # secondary-slot model on top would double count the off-fault hazard;
    # such chains are rejected up-front (logic-tree validator FDLT-013,
    # calculator guard in calculators._initialize_models), so by the time we
    # get here the secondary slot is empty and skipping it is a no-op.
    # OUTPUT BUCKETS: the aggregate contribution deliberately flows through
    # the existing principal bucket/columns (rate_principal, *_principal
    # outputs) and the distributed bucket stays exactly zero -- no new output
    # schema; consumers read the total as usual (principal + 0).
    if 'primary_fd' in adapters:
        from openquake.fdha.calc.model_adapter import (
            effective_displacement_definition)
        _pfd = adapters['primary_fd']
        if effective_displacement_definition(_pfd.model) == 'aggregate':
            W_p = location_weight(
                ctx.r,
                r_threshold_km=r_threshold_km,
                r_sigma_km=r_sigma_km,
            )
            principal_contrib = (
                rate * P_sr[:, np.newaxis] * P_fd_primary * W_p[:, np.newaxis]
            )
            return principal_contrib, distributed_contrib

    # =========================================================================
    # SECONDARY (DISTRIBUTED) CONTRIBUTION
    # =========================================================================
    if use_visini and visini_calc is not None:
        # Visini model computes combined SR × FD
        site_coords = ctx.site_coords
        if site_coords is None:
            # Fallback if coordinates not in context
            logger.warning("Site coordinates not available in context for Visini model")
            P_dist_combined = np.zeros((N_ctx, n_displ), dtype=np.float64)
        else:
            # r/x_L/L follow the Visini models' declared multi-fault
            # reference line ('segments': distance to the nearest
            # surface-reaching section, raw GC2 x/L). Single-strand
            # ruptures fall back to the canonical trace-based metrics.
            _sec_model = (calculator.secondary_surf_rup_model
                          or calculator.secondary_surf_displ_model)
            _method = getattr(_sec_model, 'MULTIFAULT_REFERENCE_LINE', 'lcp')
            r_sel, x_L_sel, L_sel = ctx.metrics_for(_method)
            # Style: an explicit model parameter wins; otherwise derive it
            # from the rupture rake, exactly like LegacyModelAdapter does.
            # The Visini coefficients are style-specific - silently
            # defaulting to 'normal' on a reverse fault shifts the FD median
            # by ~1.7x and swaps the SR occurrence tables.
            style = (
                calculator.get_model_parameters('secondary_surf_rup').get('style')
                or calculator.get_model_parameters('secondary_surf_displ').get('style')
                or str(ctx.style[0])
            )
            P_dist_combined = visini_calc.compute(
                mag=float(ctx.mag[0]),
                r=r_sel,
                rx=ctx.rx,
                L=L_sel,
                x_L=x_L_sel,
                dip=ctx.dip,
                target_displacements=target_displacements,
                sr_model=calculator.secondary_surf_rup_model,
                fd_model=calculator.secondary_surf_displ_model,
                s_sr_red_cfg=s_sr_red_cfg,
                site_coords=site_coords,
                style=style,
            )
    else:
        # Standard secondary: SR × FD
        if 'secondary_sr' in adapters:
            P_sr_sec = adapters['secondary_sr'].compute_secondary_sr(ctx, s_sr_red_cfg)
        else:
            P_sr_sec = np.zeros(N_ctx, dtype=np.float64)
        
        if 'secondary_fd' in adapters:
            P_fd_sec = adapters['secondary_fd'].compute_secondary_fd(
                ctx, target_displacements, s_sr_red_cfg,
            )
        else:
            P_fd_sec = np.zeros((N_ctx, n_displ), dtype=np.float64)
        
        # Combine: P(SR_sec) × P(FD_sec | SR_sec)
        P_dist_combined = P_sr_sec[:, np.newaxis] * P_fd_sec
    
    # =========================================================================
    # COMBINE CONTRIBUTIONS: per W_p path
    # =========================================================================
    #     lambda_principal   = rate * P_sr * P_fd_primary       * W_p(r)
    #     lambda_distributed = rate * P_sr * P_dist_combined(r) * G(r)
    #
    # W_p(r) is the probability that the site sits on the principal rupture
    # at across-strike distance r. Two separate paths (location_weight), each
    # with its own distributed weight G:
    #
    #   sigma = 0 -> W_p = boxcar |r| <= h (h = r_threshold_km) and
    #                G = 1 - W_p: the historical COMPLEMENTARY split - inside
    #                the principal zone only the principal component counts,
    #                outside it only the distributed component (Youngs 2003 /
    #                Takao 2013 per-fault either/or bookkeeping).
    #   sigma > 0 -> W_p = Petersen's pure Gaussian exp(-r^2/2 sigma^2),
    #                pinned, truncated at +-2 sigma (fixed; h plays no role)
    #                and G = 1: principal and distributed are independent and
    #                SUMMED (Petersen et al. 2011, eq. 1 + eq. 2; Fig. 10a
    #                "total hazard" = sum of its two contribution curves).
    #
    # abs() inside the helper keeps r symmetric about the trace, matching the
    # old np.abs(ctx.r) test bit-for-bit at sigma=0.
    W_p = location_weight(
        ctx.r,
        r_threshold_km=r_threshold_km,
        r_sigma_km=r_sigma_km,
    )
    if float(r_sigma_km) == 0.0:
        G = 1.0 - W_p           # complementary (legacy boxcar split, exact)
    else:
        G = np.ones_like(W_p)   # additive (Petersen eq. 1 + eq. 2)

    # Principal zone: uses primary SR and primary FD.
    # rate * P(SR_primary) * P(FD_primary | SR_primary) * W_p(r)
    principal_contrib = (
        rate * P_sr[:, np.newaxis] * P_fd_primary * W_p[:, np.newaxis]
    )

    # Distributed zone: uses secondary (distributed) models.
    # Visini's DR occurrence regressions are fit on the SURE database, which
    # contains only earthquakes with a mapped Rank-1 (principal) surface
    # rupture (Visini et al. 2025, Table 1) - i.e. P_dist_combined is already
    # conditional on the principal fault having reached the surface. It must
    # still be gated by P(SR_primary) here, same as the non-Visini branch,
    # to turn that conditional probability into a per-rupture rate
    # contribution (Visini et al. 2025 explicitly excludes both P_sr and the
    # earthquake rate from their worked example for this reason).
    distributed_contrib = (
        rate * P_sr[:, np.newaxis] * P_dist_combined * G[:, np.newaxis]
    )

    return principal_contrib, distributed_contrib