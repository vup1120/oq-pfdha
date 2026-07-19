"""
Hazard-map calculation: build a site grid around the fault sources, run
the hazard kernel and interpolate displacement at target return periods.
"""
import numpy as np
import logging
from typing import Optional, Literal

_HazardMapRateComponent = Literal["total", "principal", "distributed"]
from openquake.hazardlib import valid
from openquake.hazardlib.geo import Point
from openquake.hazardlib.site import Site, SiteCollection
from openquake.hazardlib.geo.surface.simple_fault import SimpleFaultSurface
from openquake.fdha.calc.utils.interpolation import get_map_from_curves
from openquake.fdha.calc.utils.parsing import parse_source_model_faults
from openquake.fdha.calc.config_loader import load_config, get_max_distance_km
from .hazard import calculate_fdha_hazard
from .calculators import BaseFaultRuptureCalculator

logger = logging.getLogger(__name__)


def compute_hazard_map(
    config_path: str,
    source_model_path: str,
    hdf5path: str = None,
    rupture_mesh_spacing: float = 2.0,
    complex_fault_mesh_spacing: float = None,
    width_of_mfd_bin: float = 0.1,
    return_period: float = None,
    rate_component: _HazardMapRateComponent = "total",
):
    cfg = load_config(config_path)
    geom = cfg['geometry']
    spacing = float(geom.get('region_grid_spacing', 0.01))
    max_dist = get_max_distance_km(cfg, default=10.0)
    vs30 = cfg.get('site_location', {}).get('vs30', None)

    para = cfg['parameters']
    calc = cfg.get('calculation', {})
    if return_period is None:
        return_period = float(calc.get('return_period', para.get('return_period', 100000)))

    # Create grid of points. The corner list uses the engine's job.ini
    # grammar for `region`, so parse it with the engine's validator too:
    # lon/lat range checks, duplicate-corner rejection and the engine's
    # 5-digit (~1 m) coordinate rounding, instead of a bare float() split.
    corner_coords = np.array(valid.coordinates(geom['region']))[:, :2]
    lon_vals, lat_vals = corner_coords[:, 0], corner_coords[:, 1]
    lons = np.arange(lon_vals.min(), lon_vals.max() + spacing, spacing)
    lats = np.arange(lat_vals.min(), lat_vals.max() + spacing, spacing)
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    # Generate all grid point coordinates
    grid_points = np.column_stack((lon_grid.ravel(), lat_grid.ravel()))
    # Create SiteCollection for all grid points
    if vs30 is not None:
        sites = [Site(Point(lon, lat), vs30=vs30) for lon, lat in grid_points]
    else:
        sites = [Site(Point(lon, lat)) for lon, lat in grid_points]
    sitecol = SiteCollection(sites)

    # Distance-based filtering of active grid sites
    converter_params = dict(
        rupture_mesh_spacing=rupture_mesh_spacing,
        width_of_mfd_bin=width_of_mfd_bin,
    )
    if complex_fault_mesh_spacing is not None:
        converter_params['complex_fault_mesh_spacing'] = complex_fault_mesh_spacing
    fault_sources = parse_source_model_faults(
        source_model_path,
        hdf5path=hdf5path,
        **converter_params
    )

    logger.info("Loaded %d fault sources", len(fault_sources))

    def build_surface(source) -> Optional[object]:
        """Return an OpenQuake surface instance suitable for distance queries."""
        if hasattr(source, 'surface'):
            return getattr(source, 'surface')
        if all(
            hasattr(source, attr)
            for attr in (
                'fault_trace',
                'upper_seismogenic_depth',
                'lower_seismogenic_depth',
                'dip',
            )
        ):
            try:
                return SimpleFaultSurface.from_fault_data(
                    source.fault_trace,
                    source.upper_seismogenic_depth,
                    source.lower_seismogenic_depth,
                    source.dip,
                    rupture_mesh_spacing,
                )
            except Exception as exc:
                logger.warning(
                    "Failed to build SimpleFaultSurface for source '%s': %s",
                    getattr(source, 'source_id', getattr(source, 'name', 'unknown')),
                    exc,
                )
        if hasattr(source, 'iter_ruptures'):
            try:
                if getattr(source, 'rupture_idxs', None) is not None:
                    # multiFaultSource: ruptures cover different section
                    # subsets; the site pre-filter and trace overlay need the
                    # widest footprint, not whichever rupture comes first.
                    return max((r.surface for r in source.iter_ruptures()),
                               key=lambda s: len(getattr(s, 'surfaces', ())))
                return next(source.iter_ruptures()).surface
            except StopIteration:
                logger.warning(
                    "Source '%s' has no ruptures to derive a surface",
                    getattr(source, 'source_id', getattr(source, 'name', 'unknown')),
                )
        return None

    surface_cache = {}
    for source_id, source in fault_sources.items():
        surface = build_surface(source)
        if surface is not None:
            surface_cache[source_id] = surface
        else:
            logger.warning(
                "No usable surface for source '%s' - distances will fallback to all grid sites",
                getattr(source, 'source_id', getattr(source, 'name', 'unknown')),
            )

    from openquake.fdha.calc.utils.rupture_distance import (
        VectorizedRuptureDistanceCalculator,
        _sections_info, SURFACE_DEPTH_TOLERANCE_KM,
        trace_polyline_for_source, resample_polyline,
    )

    # Sample principal-zone (on-trace) sites at the map grid resolution, so the
    # principal band stays as dense as the distributed grid regardless of the
    # (possibly coarse) ERF rupture_mesh_spacing. Convert the degree grid step
    # to km at the region's mean latitude.
    _mean_lat = float(np.mean(lats)) if len(lats) else 0.0
    trace_step_km = float(spacing) * 111.32 * max(np.cos(np.radians(_mean_lat)), 0.1)

    dist_arrays = []
    for source_id, surface in surface_cache.items():
        try:
            # 'segments' = distance to the nearest section trace: the natural
            # rupture-proximity measure for the active-site pre-filter, and
            # independent of any reference-line smoothing (an ECS/LCP line can
            # bulge away from the sections and skew the cutoff).
            calc = VectorizedRuptureDistanceCalculator(
                sitecol, surface, reference_line_method='segments')
            distances = calc.calculate_site_to_trace_distances()
            dist_arrays.append(distances)
            logger.debug("Got distances for source %s", source_id)
        except Exception as exc:
            logger.warning(
                "Failed to compute distances for source '%s': %s",
                source_id,
                exc,
            )

    # Compute minimum distance with robust handling
    if len(dist_arrays) == 0:
        logger.warning("No valid surfaces found for distance calculation. All sites will be considered active.")
        active_mask = np.ones(len(sites), dtype=bool)
    else:
        dist_matrix = np.vstack(dist_arrays)  # shape (n_sources, n_sites)
        min_dist = dist_matrix.min(axis=0)
        active_mask = min_dist <= max_dist
        logger.info("Distance filtering: %d/%d sites within %s km",
                    int(active_mask.sum()), len(sites), max_dist)
    # Extract active grid SiteCollection
    active_grid_sites = [sites[i] for i in np.where(active_mask)[0]]
    # Extract trace points at zero depth from rupture mesh (not XML fault_trace)
    # Using the actual rupture mesh ensures trace points align with the
    # rupture surfaces used in hazard calculation, avoiding misclassification
    # of on-fault sites as distributed due to coordinate misalignment.
    trace_coords = []
    for src_id, src in fault_sources.items():
        # Use the same polyline that distance calculators use, ensuring
        # consistency between trace overlay and r/x_L computation.
        # Also avoids extracting thousands of zero-depth mesh nodes for
        # fine rupture_mesh_spacing (e.g. 0.02 km on an 80 km fault).
        if src_id in surface_cache:
            surf = surface_cache[src_id]
            sections = _sections_info(surf)
            if sections is not None:
                # multi-section rupture: one top-edge trace per section;
                # buried sections produce no surface displacement, so they
                # contribute no on-trace display sites.
                for sec_lons, sec_lats, dep in sections:
                    if dep <= SURFACE_DEPTH_TOLERANCE_KM:
                        trace_coords.extend(zip(sec_lons, sec_lats))
            else:
                # Single-strand fault: use the exact original trace (mesh-
                # independent) and resample it to the grid step, so principal
                # is neither coarsened by a large rupture_mesh_spacing nor
                # over-sampled at the trace's native vertex spacing.
                polyline = trace_polyline_for_source(src, surf)
                if polyline is not None and len(polyline):
                    polyline = resample_polyline(polyline, trace_step_km)
                    trace_coords.extend((lon, lat) for lon, lat in polyline)
        elif hasattr(src, 'fault_trace'):
            polyline = resample_polyline(
                np.asarray([(pt.longitude, pt.latitude) for pt in src.fault_trace],
                           dtype=float),
                trace_step_km,
            )
            trace_coords.extend((lon, lat) for lon, lat in polyline)
    # Create Sites for trace points
    trace_sites = [
        Site(Point(lon, lat), vs30=vs30) if vs30 is not None else Site(Point(lon, lat))
        for lon, lat in trace_coords
    ]

    # Combine active grid and trace sites
    combined_sites = active_grid_sites + trace_sites

    if len(combined_sites) == 0:
        raise ValueError(
            f"No sites found within {max_dist} km of any fault. "
            f"Check your configuration: region={geom['region']}, "
            f"max_distance_km={max_dist}, and ensure your source model has valid fault geometries."
        )

    combined_sitecol = SiteCollection(combined_sites)
    logger.info("Computing hazard map for %d sites (%d grid + %d trace)",
                len(combined_sites), len(active_grid_sites), len(trace_sites))
    # Compute hazard map rates. Reuse the sources parsed above for the
    # distance pre-filter instead of re-parsing the source model XML.
    calculator = BaseFaultRuptureCalculator(
        config_path,
        [],
        hdf5path=hdf5path,
        fault_sources=fault_sources,
        **converter_params
    )
    result = calculate_fdha_hazard(calculator, combined_sitecol)
    rate_keys = {
        "total": "rates",
        "principal": "rate_principal",
        "distributed": "rate_distributed",
    }
    if rate_component not in rate_keys:
        raise ValueError(
            f"rate_component must be one of {list(rate_keys)}, got {rate_component!r}"
        )
    rates = np.asarray(result[rate_keys[rate_component]], dtype=float)
    logger.info(
        "Hazard map displacement inversion uses %s exceedance rates (%s)",
        rate_component,
        rate_keys[rate_component],
    )

    # Invert hazard curves to displacements for target return period
    imls = np.array(cfg['parameters']['target_displacement'])
    target_rate = 1.0 / return_period
    displacements = get_map_from_curves(imls, rates, target_rate)

    # Reconstruct full grid displacement map
    full_disp = np.zeros(lon_grid.size)
    full_disp[active_mask] = displacements[: active_mask.sum()]
    hazard_map = full_disp.reshape(lon_grid.shape)
    # Separate fault trace displacements
    if trace_coords:
        fault_lons, fault_lats = zip(*trace_coords)
        trace_disp = displacements[active_mask.sum():]
        return hazard_map, lons, lats, list(fault_lons), list(fault_lats), trace_disp.tolist()
    else:
        # No trace coordinates
        return hazard_map, lons, lats, [], [], []
