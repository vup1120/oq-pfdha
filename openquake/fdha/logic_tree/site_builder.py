"""Build the map-mode SiteCollection (grid + fault-trace + distance filter).

This replicates the site-building logic inside ``hazard_map.compute_hazard_map``
so the logic-tree driver can obtain per-site annual rates on the *same* grid
used by the single-branch map, without calling ``compute_hazard_map`` itself
(the driver needs raw rates, not displacements).

The sole scientific kernel call stays as ``calculate_fdha_hazard`` in
``hazard.py`` (read-only per the v3 integration rules).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Optional

import numpy as np

from openquake.fdha.calc.utils.parsing import parse_source_model_faults
from openquake.fdha.calc.utils.rupture_distance import (
    VectorizedRuptureDistanceCalculator,
    resample_polyline,
    trace_polyline_for_source,
)
from openquake.hazardlib.geo import Point
from openquake.hazardlib.geo.surface.simple_fault import SimpleFaultSurface
from openquake.hazardlib.site import Site, SiteCollection

logger = logging.getLogger(__name__)


@dataclass
class HazardMapSites:
    combined_sitecol: Any  # SiteCollection
    active_mask: np.ndarray  # (n_grid,) bool
    lon_grid: np.ndarray  # (n_lat, n_lon)
    lat_grid: np.ndarray  # (n_lat, n_lon)
    lons: np.ndarray  # 1D grid lons
    lats: np.ndarray  # 1D grid lats
    trace_coords: list[tuple[float, float]]
    n_active_grid: int
    n_trace: int


def _build_surface(source: Any, rupture_mesh_spacing: float) -> Optional[Any]:
    if hasattr(source, "surface"):
        return source.surface
    needed = ("fault_trace", "upper_seismogenic_depth", "lower_seismogenic_depth", "dip")
    if all(hasattr(source, a) for a in needed):
        try:
            return SimpleFaultSurface.from_fault_data(
                source.fault_trace,
                source.upper_seismogenic_depth,
                source.lower_seismogenic_depth,
                source.dip,
                rupture_mesh_spacing,
            )
        except Exception as exc:  # pragma: no cover - defensive
            logger.warning(
                "Failed to build SimpleFaultSurface for source '%s': %s",
                getattr(source, "source_id", getattr(source, "name", "?")),
                exc,
            )
    if hasattr(source, "iter_ruptures"):
        try:
            return next(source.iter_ruptures()).surface
        except StopIteration:
            return None
    return None


def build_hazard_map_sites(
    *,
    region: str,
    spacing: float,
    max_distance_km: float,
    vs30: Optional[float],
    source_model_paths: Any,
    hdf5path: Optional[str] = None,
    rupture_mesh_spacing: float = 2.0,
    complex_fault_mesh_spacing: Optional[float] = None,
    width_of_mfd_bin: float = 0.1,
    fault_sources: Optional[dict] = None,
) -> HazardMapSites:
    corner_coords = np.array([
        list(map(float, p.strip().split())) for p in region.split(",")
    ])
    lon_vals, lat_vals = corner_coords[:, 0], corner_coords[:, 1]
    lons = np.arange(lon_vals.min(), lon_vals.max() + spacing, spacing)
    lats = np.arange(lat_vals.min(), lat_vals.max() + spacing, spacing)
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    grid_points = np.column_stack((lon_grid.ravel(), lat_grid.ravel()))

    if vs30 is not None:
        grid_sites = [Site(Point(float(lon), float(lat)), vs30=float(vs30))
                      for lon, lat in grid_points]
    else:
        grid_sites = [Site(Point(float(lon), float(lat))) for lon, lat in grid_points]

    converter_params: dict[str, Any] = dict(
        rupture_mesh_spacing=rupture_mesh_spacing,
        width_of_mfd_bin=width_of_mfd_bin,
    )
    if complex_fault_mesh_spacing is not None:
        converter_params["complex_fault_mesh_spacing"] = complex_fault_mesh_spacing

    if fault_sources is None:
        fault_sources = parse_source_model_faults(
            source_model_paths, hdf5path=hdf5path or "", **converter_params
        )

    surface_cache: dict[str, Any] = {}
    for src_id, src in fault_sources.items():
        surf = _build_surface(src, rupture_mesh_spacing)
        if surf is not None:
            surface_cache[src_id] = surf

    dist_arrays = []
    sitecol_grid = SiteCollection(grid_sites)
    for src_id, surface in surface_cache.items():
        try:
            dcalc = VectorizedRuptureDistanceCalculator(sitecol_grid, surface)
            dist_arrays.append(dcalc.calculate_site_to_trace_distances())
        except Exception as exc:  # pragma: no cover - defensive
            logger.warning("distance fail for %s: %s", src_id, exc)

    if not dist_arrays:
        active_mask = np.ones(len(grid_sites), dtype=bool)
    else:
        dist_matrix = np.vstack(dist_arrays)
        active_mask = dist_matrix.min(axis=0) <= max_distance_km

    active_grid_sites = [grid_sites[i] for i in np.where(active_mask)[0]]

    # Principal-zone (on-trace) sites are sampled along strike at the map's own
    # grid resolution (``region_grid_spacing``), so the principal band matches
    # the distributed grid it overlays and stays independent of the (possibly
    # coarse) ERF ``rupture_mesh_spacing``. ``spacing`` is in degrees; convert
    # to km at the region's mean latitude (longitudes shrink by cos φ, the
    # finer of the two grid axes, keeping trace sites at least grid-dense).
    #
    # Resample (not merely densify): NRML traces are often digitised at
    # sub-kilometre vertex spacing, so keeping every native vertex would place
    # thousands of principal sites on a grid that cannot resolve them — the 48
    # onshore Taiwan faults gave 3 246 trace sites against a 1 040-site 0.1°
    # grid, a ~36x cost for no extra information.
    mean_lat = float(np.mean(lats)) if len(lats) else 0.0
    trace_step_km = float(spacing) * 111.32 * max(np.cos(np.radians(mean_lat)), 0.1)

    trace_coords: list[tuple[float, float]] = []
    for src_id, src in fault_sources.items():
        polyline = trace_polyline_for_source(src, surface_cache.get(src_id))
        if polyline is None or len(polyline) == 0:
            continue
        polyline = resample_polyline(polyline, trace_step_km)
        trace_coords.extend((float(lon), float(lat)) for lon, lat in polyline)

    trace_sites = [
        Site(Point(lon, lat), vs30=vs30) if vs30 is not None else Site(Point(lon, lat))
        for lon, lat in trace_coords
    ]
    combined = active_grid_sites + trace_sites
    if not combined:
        raise ValueError(
            f"No sites within {max_distance_km} km of any fault trace in region "
            f"{region!r}; check grid / max_distance_km / source geometry."
        )

    return HazardMapSites(
        combined_sitecol=SiteCollection(combined),
        active_mask=active_mask,
        lon_grid=lon_grid,
        lat_grid=lat_grid,
        lons=lons,
        lats=lats,
        trace_coords=trace_coords,
        n_active_grid=len(active_grid_sites),
        n_trace=len(trace_coords),
    )
