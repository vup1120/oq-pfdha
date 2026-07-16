# -*- coding: utf-8 -*-
"""
Unified FDHA hazard calculation module.

This module provides a single implementation for both hazard curve and
hazard map calculations, eliminating code duplication and ensuring
consistent results.
"""

import numpy as np
from tqdm import tqdm
import logging
from typing import Dict, Any, Optional, List, Tuple, TYPE_CHECKING

from openquake.fdha.calc.contexts import FDHAContext, FDHAContextMaker
from openquake.hazardlib.site import SiteCollection

if TYPE_CHECKING:
    from openquake.fdha.calc.calculators import BaseFaultRuptureCalculator

logger = logging.getLogger(__name__)

# Model names that indicate Visini et al. (2025) implementation
VISINI_SR_NAMES = {'Visini2025SecondarySR'}
VISINI_FD_NAMES = {'Visini2025SecondaryFD'}


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
        Dictionary containing:
        - 'imls': List of displacement levels (m)
        - 'poes': Exceedance rates array, shape (n_sites, n_displ)
        - 'rate_principal': Principal contribution, shape (n_sites, n_displ)
        - 'rate_distributed': Distributed contribution, shape (n_sites, n_displ)
        - 'each_fault': Per-fault contributions dict
        - 'site_lons': Site longitudes
        - 'site_lats': Site latitudes
    """
    # Use provided sitecol or calculator's sitecol
    if sitecol is None:
        sitecol = calculator.sitecol
    
    n_sites = len(sitecol)
    
    target_displacements = calculator.target_displacements
    n_displ = len(target_displacements)
    
    logger.info(f"Starting FDHA hazard calculation: {n_sites} sites, {n_displ} displacement levels")
    
    # Get max_distance from config (check multiple sections)
    max_dist = 50.0  # default
    for section in ['calculation', 'parameters', 'erf']:
        section_cfg = calculator.config.get(section, {})
        if 'max_distance_km' in section_cfg:
            max_dist = float(section_cfg['max_distance_km'])
            break
    
    # Create context maker with caching
    sitecol_for_cmaker = sitecol
    
    cmaker = FDHAContextMaker(
        sitecol=sitecol_for_cmaker,
        fdha_params=calculator.get_fdha_params(),
        maximum_distance=max_dist,
    )
    
    # Check for Visini models
    use_visini, visini_calc = _setup_visini_calculator(calculator)
    
    # Initialize rate accumulators
    rate_principal = np.zeros((n_sites, n_displ), dtype=np.float64)
    rate_distributed = np.zeros((n_sites, n_displ), dtype=np.float64)
    each_fault: Dict[str, List] = {}
    
    # Get model adapters
    adapters = calculator.adapters
    
    # Debug: log available adapters
    logger.debug(f"Available adapters: {list(adapters.keys())}")
    if not adapters:
        logger.warning("No adapters available - models may not have loaded correctly")
    for name, adapter in adapters.items():
        logger.debug(f"  Adapter {name}: model={adapter.model.__class__.__name__}")
    
    # Reduction configs
    p_sr_red_cfg = calculator.p_sr_red_cfg
    s_sr_red_cfg = calculator.s_sr_red_cfg
    r_threshold_km = calculator.r_threshold_km

    # Rupture-location uncertainty (W_p).
    # See docs/design/rupture_location_uncertainty.md, section 2.
    r_sigma_km = calculator.r_sigma_km

    # PMF time span for non-parametric (multiFaultSource) ruptures:
    # get_ctx converts their probs_occur into a Poisson-equivalent annual
    # rate over this investigation time. Parametric ruptures ignore it.
    # The authoritative value is per source (parsed from the NRML header);
    # the INI value is a fallback and must not silently contradict it.
    ini_investigation_time = calculator.config.get(
        'calculation', {}).get('investigation_time')

    # Process each fault source
    # Convert to list immediately to avoid iterator exhaustion issues
    fault_sources_list = list(calculator.fault_sources.values())

    iterator = tqdm(fault_sources_list, desc="Computing hazard") if show_progress else fault_sources_list

    for src in iterator:
        logger.info(f"Processing fault: {src.name}")
        investigation_time = _resolve_investigation_time(
            src, ini_investigation_time)
        
        fault_rate = np.zeros((n_sites, n_displ), dtype=np.float64)
        n_ruptures = 0
        n_surface_rupturing = 0
        
        # Convert iter_ruptures to list to avoid iterator exhaustion issues
        ruptures_list = list(src.iter_ruptures())
        
        for rup in ruptures_list:
            n_ruptures += 1
            
            # Skip buried ruptures
            if not cmaker.is_surface_rupturing(rup):
                continue
            n_surface_rupturing += 1
            
            # Create context (returns None if all sites too far)
            ctx = cmaker.get_ctx(rup, investigation_time=investigation_time)
            if ctx is None:
                continue
            
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
            
            # Accumulate by site ID using vectorized operations
            sids = ctx.sids
            # Ensure sids length matches principal_contrib first dimension
            if len(sids) != principal_contrib.shape[0]:
                # If mismatch, use first site ID for all contributions
                if principal_contrib.shape[0] == 1:
                    if len(sids) > 0:
                        sids = sids[:1]  # Take first site ID only
                    else:
                        # sids is empty but principal_contrib has shape (1, n_displ)
                        # This means all sites were filtered out, skip accumulation
                        continue
                else:
                    raise ValueError(f"sids length ({len(sids)}) doesn't match principal_contrib shape[0] ({principal_contrib.shape[0]})")
            
            # Skip if sids is empty (all sites filtered out)
            if len(sids) == 0:
                continue
                
            valid_mask = sids < n_sites
            
            if np.all(valid_mask):
                # For single site case, use direct array addition (matching old implementation)
                # This is more efficient and avoids potential broadcasting issues with np.add.at
                if n_sites == 1 and len(sids) == 1:
                    # Direct array addition for single site (matches old implementation)
                    rate_principal[0] += principal_contrib[0]
                    rate_distributed[0] += distributed_contrib[0]
                    fault_rate[0] += principal_contrib[0] + distributed_contrib[0]
                else:
                    # Multiple sites - use scatter-add
                    np.add.at(rate_principal, sids, principal_contrib)
                    np.add.at(rate_distributed, sids, distributed_contrib)
                    np.add.at(fault_rate, sids, principal_contrib + distributed_contrib)
            else:
                # Filter to valid sids only
                valid_sids = sids[valid_mask]
                np.add.at(rate_principal, valid_sids, principal_contrib[valid_mask])
                np.add.at(rate_distributed, valid_sids, distributed_contrib[valid_mask])
                np.add.at(fault_rate, valid_sids, (principal_contrib + distributed_contrib)[valid_mask])
        
        logger.debug(f"  {src.name}: {n_ruptures} ruptures, {n_surface_rupturing} surface-rupturing")
        
        each_fault[f"fault_{src.name}"] = fault_rate.tolist()
    
    # Log cache statistics
    stats = cmaker.get_cache_stats()
    logger.info(
        f"Distance cache: {stats['hits']} hits, {stats['misses']} misses, "
        f"{stats['hit_rate']:.1%} hit rate"
    )
    
    # Total exceedance rate
    total_rate = rate_principal + rate_distributed
    
    # Extract site coordinates
    site_lons = cmaker._lons.tolist()
    site_lats = cmaker._lats.tolist()
    
    return {
        # Primary keys (new format)
        'imls': target_displacements.tolist(),
        'poes': total_rate.tolist(),
        'rate_principal': rate_principal.tolist(),
        'rate_distributed': rate_distributed.tolist(),
        'each_fault': each_fault,
        'site_lons': site_lons,
        'site_lats': site_lats,
        'n_sites': n_sites,
        'n_displ': n_displ,
        # Backward compatibility aliases (old hazard_map_calculator.py format)
        'displacements': target_displacements.tolist(),
        'annual_rate_total': total_rate.tolist(),
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
    value, then 1.0 — their ruptures have explicit annual rates and ignore
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
    
    sr_name = sr_model.__class__.__name__ if sr_model else ''
    fd_name = fd_model.__class__.__name__ if fd_model else ''
    
    use_visini = sr_name in VISINI_SR_NAMES or fd_name in VISINI_FD_NAMES
    
    if not use_visini:
        return False, None
    
    logger.info(f"Using Visini model: SR={sr_name}, FD={fd_name}")
    
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

    Principal and distributed are independent contributions of the same
    surface-rupturing event (both carry rate * P_sr) and are SUMMED by the
    caller (Petersen et al. 2011 eq. 1 + eq. 2; Fig. 10a "total hazard").

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
        P_sr = adapters['primary_sr'].compute_primary_sr(ctx, p_sr_red_cfg)
        if P_sr is None:
            return principal_contrib, distributed_contrib
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
        if P_fd_primary is None:
            P_fd_primary = np.zeros((N_ctx, n_displ), dtype=np.float64)
    else:
        P_fd_primary = np.zeros((N_ctx, n_displ), dtype=np.float64)
    
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
            # The Visini coefficients are style-specific — silently
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
    # COMBINE CONTRIBUTIONS: PRINCIPAL * W_p + DISTRIBUTED
    # =========================================================================
    # Both contributions belong to the same surface-rupturing event and are
    # independent, so they are summed (Petersen et al. 2011, eq. 1 + eq. 2;
    # Fig. 10a "total hazard" = sum of its two contribution curves):
    #     lambda_principal   = rate * P_sr * P_fd_primary * W_p(r)
    #     lambda_distributed = rate * P_sr * P_dist_combined(r)
    #
    # W_p(r) is the probability that the site sits on the principal rupture
    # at across-strike distance r. Two separate paths (location_weight):
    # sigma=0 -> the boxcar |r| <= h (h = r_threshold_km); sigma>0 ->
    # Petersen's pure Gaussian exp(-r^2/2 sigma^2), pinned, truncated at
    # +-2 sigma (fixed), with h playing no role. abs() inside the helper
    # keeps r symmetric about the trace, matching the old np.abs(ctx.r)
    # test bit-for-bit at sigma=0.
    W_p = location_weight(
        ctx.r,
        r_threshold_km=r_threshold_km,
        r_sigma_km=r_sigma_km,
    )

    # Principal zone: uses primary SR and primary FD.
    # rate * P(SR_primary) * P(FD_primary | SR_primary) * W_p(r)
    principal_contrib = (
        rate * P_sr[:, np.newaxis] * P_fd_primary * W_p[:, np.newaxis]
    )

    # Distributed zone: uses secondary (distributed) models.
    # Visini's DR occurrence regressions are fit on the SURE database, which
    # contains only earthquakes with a mapped Rank-1 (principal) surface
    # rupture (Visini et al. 2025, Table 1) — i.e. P_dist_combined is already
    # conditional on the principal fault having reached the surface. It must
    # still be gated by P(SR_primary) here, same as the non-Visini branch,
    # to turn that conditional probability into a per-rupture rate
    # contribution (Visini et al. 2025 explicitly excludes both P_sr and the
    # earthquake rate from their worked example for this reason).
    distributed_contrib = (
        rate * P_sr[:, np.newaxis] * P_dist_combined
    )

    return principal_contrib, distributed_contrib


def rates_to_poes(
    rates: np.ndarray,
    investigation_time: float = 1.0
) -> np.ndarray:
    """
    Convert annual exceedance rates to probabilities of exceedance.
    
    Uses Poisson model: POE = 1 - exp(-rate × time)
    
    Args:
        rates: Annual exceedance rates
        investigation_time: Time period in years
        
    Returns:
        Probabilities of exceedance
    """
    return 1.0 - np.exp(-rates * investigation_time)


def poes_to_rates(
    poes: np.ndarray,
    investigation_time: float = 1.0
) -> np.ndarray:
    """
    Convert probabilities of exceedance to annual rates.
    
    Inverse of Poisson model: rate = -ln(1 - POE) / time
    
    Args:
        poes: Probabilities of exceedance
        investigation_time: Time period in years
        
    Returns:
        Annual exceedance rates
    """
    # Clip POE to avoid log(0)
    poes_clipped = np.clip(poes, 0.0, 1.0 - 1e-10)
    return -np.log(1.0 - poes_clipped) / investigation_time


def interpolate_hazard_map(
    rates: np.ndarray,
    imls: np.ndarray,
    return_period: float
) -> np.ndarray:
    """
    Interpolate hazard to specific return period.
    
    Uses vectorized operations for efficient computation across all sites.
    
    Args:
        rates: Exceedance rates, shape (n_sites, n_imls)
        imls: Intensity measure levels (displacement values)
        return_period: Target return period in years
        
    Returns:
        Interpolated displacement values for each site
    """
    target_rate = 1.0 / return_period
    n_sites = rates.shape[0]
    imls = np.asarray(imls)
    
    # Pre-compute log values (vectorized)
    log_imls = np.log(imls)
    log_imls_reversed = log_imls[::-1]
    log_target = np.log(target_rate)
    
    # Get first and last rates for all sites (vectorized)
    first_rates = rates[:, 0]
    last_rates = rates[:, -1]
    
    # Initialize result
    result = np.zeros(n_sites, dtype=np.float64)
    
    # Case 1: Rate higher than all IMLs - use min IML
    mask_high = target_rate >= first_rates
    result[mask_high] = imls[0]
    
    # Case 2: Rate lower than all IMLs - use max IML
    mask_low = target_rate <= last_rates
    result[mask_low] = imls[-1]
    
    # Case 3: Interpolate in log-log space
    mask_interp = ~mask_high & ~mask_low
    
    if np.any(mask_interp):
        # Process interpolation sites
        interp_indices = np.where(mask_interp)[0]
        
        # Vectorized log-space interpolation
        for idx in interp_indices:
            site_rates = rates[idx, :]
            log_rates_reversed = np.log(site_rates[::-1] + 1e-30)
            result[idx] = np.exp(np.interp(log_target, log_rates_reversed, log_imls_reversed))
    
    return result