"""Plotting and JSON-export helpers for hazard curves and hazard maps."""
import os
import numpy as np
import matplotlib.pyplot as plt
import logging
import json
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from typing import Union, List, Dict, Any, Optional
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatter, LogLocator, FuncFormatter

logger = logging.getLogger(__name__)

def plot_hazard_map(
    hazard_map: np.ndarray,
    lons: np.ndarray,
    lats: np.ndarray,
    fault_lons: Union[List[float], np.ndarray],
    fault_lats: Union[List[float], np.ndarray],
    trace_disp: Union[List[float], np.ndarray],
    plot_file: Optional[str] = None,
    title: Optional[str] = None,
) -> None:
    """
    Plot a fault-displacement hazard map.

    :param hazard_map: 2D array of shape (n_rows, n_cols) with displacement (m) values
    :param lons: 1D array of longitudes for the grid columns
    :param lats: 1D array of latitudes for the grid rows
    :param fault_lons: list or 1D array of longitudes for fault-trace sites
    :param fault_lats: list or 1D array of latitudes for fault-trace sites
    :param trace_disp: list or 1D array of displacement (m) values at fault-trace sites
    :param plot_file: optional path to save the plot (if None, show interactively)
    :param title: optional figure title (default: generic hazard map title)
    """
    # Ensure numpy arrays for masking
    fault_lons = np.array(fault_lons)
    fault_lats = np.array(fault_lats)
    trace_disp = np.array(trace_disp)

    # Flatten grid for scatter plotting
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    disp_flat = hazard_map.flatten()

    # Mask out zero-displacement cells
    grid_mask = disp_flat > 0
    grid_lons = lon_grid.flatten()[grid_mask]
    grid_lats = lat_grid.flatten()[grid_mask]
    grid_disp = disp_flat[grid_mask]

    # Mask trace sites
    trace_mask = trace_disp > 0
    trace_lons = fault_lons[trace_mask]
    trace_lats = fault_lats[trace_mask]
    trace_vals = trace_disp[trace_mask]

    # Nothing to plot?
    if grid_disp.size == 0 and trace_vals.size == 0:
        logger.warning("No non-zero displacements to plot.")
        return

    # Color normalization
    all_vals = np.concatenate([grid_disp, trace_vals]) if trace_vals.size > 0 else grid_disp
    norm = LogNorm(vmin=max(1e-5, np.min(all_vals)), vmax=np.max(all_vals))
    
    # Plot setup with larger figure size
    fig = plt.figure(figsize=(14, 10))
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    # Enhanced geographic background with high-resolution NaturalEarth features
    try:
        # Use high-resolution (10m) NaturalEarth features for better coastline quality
        ocean_feature = cfeature.NaturalEarthFeature(
            'physical', 'ocean', '10m',
            edgecolor='none', facecolor='#cce5ff', zorder=0
        )
        land_feature = cfeature.NaturalEarthFeature(
            'physical', 'land', '10m',
            edgecolor='none', facecolor='#f5f5f5', zorder=0
        )
        ax.add_feature(ocean_feature)
        ax.add_feature(land_feature)
    except Exception as e:
        # Fallback to medium resolution if 10m not available
        try:
            ocean_feature = cfeature.NaturalEarthFeature(
                'physical', 'ocean', '50m',
                edgecolor='none', facecolor='#cce5ff', zorder=0
            )
            land_feature = cfeature.NaturalEarthFeature(
                'physical', 'land', '50m',
                edgecolor='none', facecolor='#f5f5f5', zorder=0
            )
            ax.add_feature(ocean_feature)
            ax.add_feature(land_feature)
        except Exception as e2:
            # Final fallback to simple features
            logger.debug(f"NaturalEarth features not available, using simple features: {e2}")
            ax.add_feature(cfeature.OCEAN, facecolor='#cce5ff', zorder=0)
            ax.add_feature(cfeature.LAND, facecolor='#f5f5f5', zorder=0)
    
    # High-resolution coastline and borders
    try:
        # Try high-resolution coastline first
        coastline_10m = cfeature.NaturalEarthFeature(
            'physical', 'coastline', '10m',
            edgecolor='#2c3e50', facecolor='none', linewidth=1.0, zorder=2
        )
        ax.add_feature(coastline_10m)
    except Exception:
        # Fallback to standard coastline
        ax.add_feature(cfeature.COASTLINE, linewidth=1.0, edgecolor='#2c3e50', zorder=2)
    
    try:
        # Try high-resolution borders
        borders_10m = cfeature.NaturalEarthFeature(
            'cultural', 'admin_0_boundary_lines_land', '10m',
            edgecolor='#7f8c8d', facecolor='none', linestyle=':', linewidth=0.6, zorder=2
        )
        ax.add_feature(borders_10m)
    except Exception:
        # Fallback to standard borders
        ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.6, edgecolor='#7f8c8d', zorder=2)

    # Use professional colormap (viridis for better perception)
    cmap = plt.cm.get_cmap('viridis')
    
    # Combine grid sites and fault trace sites with same plot style
    all_lons = grid_lons
    all_lats = grid_lats
    all_disp = grid_disp
    
    if trace_vals.size > 0:
        # Combine fault trace sites with grid sites
        all_lons = np.concatenate([grid_lons, trace_lons])
        all_lats = np.concatenate([grid_lats, trace_lats])
        all_disp = np.concatenate([grid_disp, trace_vals])
    
    # Scatter all sites with unified styling
    grid_scatter = ax.scatter(
        all_lons, all_lats, c=all_disp,
        transform=ccrs.PlateCarree(),
        norm=norm, cmap=cmap, s=15, edgecolors='none',
        alpha=0.8, zorder=3
    )

    # Enhanced colorbar with more tick labels
    cbar = plt.colorbar(grid_scatter, ax=ax, extend='both', pad=0.02, shrink=0.8)
    cbar.set_label('Displacement (m)', fontsize=13, fontweight='bold', labelpad=15)
    
    # Calculate appropriate number of ticks based on data range
    vmin = norm.vmin
    vmax = norm.vmax
    log_range = np.log10(vmax) - np.log10(vmin)
    
    # Set more tick labels - aim for ~10-15 major ticks
    num_ticks = max(10, min(20, int(log_range * 3) + 1))
    cbar.ax.yaxis.set_major_locator(LogLocator(base=10, numticks=num_ticks))
    cbar.ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
    
    # Custom formatter: use simplified scientific notation (10^-1 instead of 1×10^-1)
    def scientific_formatter(x, pos):
        """Format all values in simplified scientific notation"""
        if x == 0:
            return '0'
        # Use scientific notation - only show 10^exp, skip mantissa if it's 1
        exp = int(np.floor(np.log10(abs(x))))
        mantissa = x / (10 ** exp)
        
        # If mantissa is close to 1, just show 10^exp
        if abs(mantissa - 1.0) < 0.01:
            if exp == 0:
                return '1'
            else:
                return f'10$^{{{exp}}}$'
        # Otherwise show mantissa × 10^exp
        elif exp == 0:
            return f'{mantissa:.1f}'
        else:
            return f'{mantissa:.1f}×10$^{{{exp}}}$'
    
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(scientific_formatter))
    
    # Improve tick label appearance: bold, larger font
    cbar.ax.tick_params(labelsize=13, length=4, width=1)
    cbar.ax.tick_params(which='minor', length=2, width=0.5)
    # Make labels bold
    for label in cbar.ax.yaxis.get_ticklabels():
        label.set_fontweight('bold')
        label.set_fontsize(13)
    
    # Enhanced axis labels with degree symbols
    ax.set_xlabel('Longitude (°)', fontsize=12, fontweight='bold', labelpad=10)
    ax.set_ylabel('Latitude (°)', fontsize=12, fontweight='bold', labelpad=10)
    
    # Enhanced title
    ax.set_title(
        title if title is not None else 'Fault Displacement Hazard Map',
        fontsize=15,
        fontweight='bold',
        pad=20,
    )
    
    # Enhanced gridlines
    gl = ax.gridlines(draw_labels=True, linestyle='--', alpha=0.6, linewidth=0.8, 
                      color='gray', zorder=1)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 10, 'weight': 'normal'}
    gl.ylabel_style = {'size': 10, 'weight': 'normal'}
    
    # No legend needed since all sites use same style

    # Set extent
    all_lons = np.concatenate([lons, fault_lons])
    all_lats = np.concatenate([lats, fault_lats])
    ax.set_extent([all_lons.min() - 0.5, all_lons.max() + 0.5,
                   all_lats.min() - 0.5, all_lats.max() + 0.5],
                  crs=ccrs.PlateCarree())

    if plot_file:
        # Save with high resolution (at least 600 DPI as requested)
        plt.savefig(plot_file, dpi=600, bbox_inches='tight', facecolor='white')
        logger.info(f"Hazard map plot saved to {plot_file} at 600 DPI")
    else:
        plt.show()

    plt.close()

def _save_single_curve(
    imls: np.ndarray,
    poe: np.ndarray,
    plot_file: Optional[str],
    title_suffix: str = "",
) -> None:
    """Render a single hazard curve and save (or show) it."""
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.loglog(
        imls, poe,
        color='C0', linestyle='-', linewidth=2, marker='o', markersize=4,
        label="Hazard Curve", alpha=0.8,
    )
    ax.set_xlabel("Displacement (m)", fontsize=12, fontweight='bold')
    ax.set_ylabel("Annual Exceedance Rate", fontsize=12, fontweight='bold')
    title = "Fault Displacement Hazard Curve"
    if title_suffix:
        title = f"{title} - {title_suffix}"
    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
    ax.grid(True, which='major', linestyle='-', alpha=0.3, linewidth=0.8)
    ax.grid(True, which='minor', linestyle='--', alpha=0.2, linewidth=0.5)
    plt.tight_layout()

    if plot_file:
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        logger.info(f"Hazard curve plot saved to {plot_file}")
    else:
        plt.show()
    plt.close()


def plot_fault_displacement_hazard(
    results: Dict[str, Any],
    plot_file: Optional[str] = None
) -> None:
    """
    Plot fault displacement hazard curves.

    For single-site results, produces one figure (identical to legacy behaviour).
    For multi-site results (N > 1) with ``plot_file`` set, produces N separate
    figures with ``_site{i}`` inserted before the extension
    (e.g. ``hazard_curve.png`` -> ``hazard_curve_site0.png``,
    ``hazard_curve_site1.png``).  Each figure carries a label with the site's
    longitude and latitude when available.

    :param results: Dictionary containing 'imls' and 'poes' arrays
    :param plot_file: optional path to save the plot (if None, show interactively)
    """
    imls = np.array(results["imls"])
    poes = np.array(results["poes"])

    # Handle both 1D and 2D poes arrays
    if poes.ndim == 1:
        poes = poes.reshape(1, -1)

    site_lons = results.get("site_lons") or []
    site_lats = results.get("site_lats") or []

    n_sites = len(poes)

    if n_sites == 1:
        # Single-site path - unchanged behaviour
        _save_single_curve(imls, poes[0], plot_file, title_suffix="")
        return

    # Multi-site path - one figure per site
    if plot_file:
        base, ext = os.path.splitext(plot_file)
        for idx, poe in enumerate(poes):
            per_site_path = f"{base}_site{idx}{ext}"
            if idx < len(site_lons) and idx < len(site_lats):
                suffix = f"site {idx} (lon={site_lons[idx]:.4f}, lat={site_lats[idx]:.4f})"
            else:
                suffix = f"site {idx}"
            _save_single_curve(imls, poe, per_site_path, title_suffix=suffix)
    else:
        # Interactive display: still show each site in its own figure
        for idx, poe in enumerate(poes):
            if idx < len(site_lons) and idx < len(site_lats):
                suffix = f"site {idx} (lon={site_lons[idx]:.4f}, lat={site_lats[idx]:.4f})"
            else:
                suffix = f"site {idx}"
            _save_single_curve(imls, poe, None, title_suffix=suffix)

def save_results_to_json(
    results: Dict[str, Any],
    output_file: str
) -> None:
    """
    Serialise hazard-curve results to JSON.

    For multi-site results, adds a ``sites`` list indexing each curve by its
    site id (``sid``) with ``(lon, lat)`` metadata. Row ``i`` of ``poes`` /
    ``rate_principal`` / ``rate_distributed`` corresponds to ``sites[i]``.
    Existing keys (``poes``, ``site_lons``, ``site_lats``, ...) are left
    unchanged, so single-site consumers are unaffected.
    """
    # Additive per-site index: does not mutate the caller's dict
    payload = dict(results)
    if 'site_lons' in payload and 'site_lats' in payload and 'poes' in payload:
        poes_arr = np.asarray(payload['poes'])
        site_lons = payload['site_lons']
        site_lats = payload['site_lats']
        if poes_arr.ndim == 2 and len(site_lons) == poes_arr.shape[0]:
            payload['sites'] = [
                {'sid': i, 'lon': site_lons[i], 'lat': site_lats[i]}
                for i in range(poes_arr.shape[0])
            ]

    with open(output_file, "w") as f:
        json.dump(payload, f, indent=2)
    logger.info(f"Hazard curve results saved to {output_file}")

def save_map_to_json(
    hazard_map: np.ndarray,
    lons: np.ndarray,
    lats: np.ndarray,
    fault_lons: Union[List[float], np.ndarray],
    fault_lats: Union[List[float], np.ndarray],
    trace_disp: Union[List[float], np.ndarray],
    output_file: str
) -> None:
    """Save hazard-map data to JSON."""
    data = {
        "hazard_map": hazard_map.tolist(),
        "lons": lons.tolist(),
        "lats": lats.tolist(),
        "fault_lons": np.array(fault_lons).tolist(),
        "fault_lats": np.array(fault_lats).tolist(),
        "trace_disp": np.array(trace_disp).tolist()
    }
    with open(output_file, "w") as f:
        json.dump(data, f, indent=2)
    logger.info(f"Hazard map data saved to {output_file}")