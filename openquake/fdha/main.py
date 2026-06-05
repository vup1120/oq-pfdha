# -*- coding: utf-8 -*-
"""
FDHA Command Line Interface.

OpenQuake-style invocation:

  fdha job.ini
  fdha job.ini --plot
  fdha job.ini --output results.json

"""

import os
import sys
import json
import argparse
import logging
from pathlib import Path
from typing import Optional, Dict, Any, TYPE_CHECKING

# Apply the h3 v3/v4 compatibility shim before any OpenQuake (hazardlib)
# module that uses h3 is imported by the calculation path.
from openquake.fdha import h3_compat  # noqa: E402,F401

if TYPE_CHECKING:  # avoid importing heavy modules at CLI startup
    from openquake.fdha.logic_tree.driver import LogicTreeResult

logger = logging.getLogger(__name__)


def setup_logging(verbose: bool = False):
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def parse_config_file(config_path: str) -> Dict[str, Any]:
    """
    Parse configuration file (INI or TOML) using unified loader.
    
    Args:
        config_path: Path to configuration file
        
    Returns:
        Configuration dictionary
    """
    from openquake.fdha.calc.config_loader import load_config
    return load_config(config_path)


def run_calculation(
    config_path: str,
    output_path: Optional[str] = None,
    plot: bool = False,
    plot_file: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Run FDHA calculation.
    
    Args:
        config_path: Path to configuration file
        output_path: Optional output file path
        plot: Whether to generate plots
        plot_file: Optional plot file path
        
    Returns:
        Results dictionary
    """
    # Parse config
    config = parse_config_file(config_path)

    # Logic-tree dispatch (new path)
    from openquake.fdha.calc.config_loader import calculation_requests_fdha_logic_tree

    calc_cfg = config.get("calculation", {})
    has_fdha_lt = (
        calculation_requests_fdha_logic_tree(calc_cfg) if isinstance(calc_cfg, dict) else False
    )

    if has_fdha_lt:
        return _run_logic_tree(
            config_path=config_path,
            output_path=output_path,
            plot=plot,
            plot_file=plot_file,
        )

    from openquake.fdha.calc.config_loader import ConfigurationError

    raise ConfigurationError(
        "Configuration must provide canonical [calculation].fdha_logic_tree_file "
        "and [calculation].source_model_logic_tree_file."
    )


def save_results(results: Dict[str, Any], output_path: str):
    """Save results to file."""
    import numpy as np
    
    ext = os.path.splitext(output_path)[1].lower()
    
    if ext == '.json':
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
    elif ext == '.csv':
        imls = results.get('imls', [])
        poes = np.array(results.get('poes', []))
        
        with open(output_path, 'w') as f:
            f.write('site_id,' + ','.join(f'd={d}' for d in imls) + '\n')
            for i, row in enumerate(poes):
                f.write(f'{i},' + ','.join(f'{v:.6e}' for v in row) + '\n')
    else:
        # Default to JSON
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)


def _run_logic_tree(
    config_path: str,
    output_path: Optional[str],
    plot: bool,
    plot_file: Optional[str],
) -> Dict[str, Any]:
    """Run the logic-tree driver and collect convenience paths for the CLI."""
    from openquake.fdha.logic_tree.driver import FdhaLogicTree

    lt = FdhaLogicTree.from_ini(config_path)
    result = lt.run()
    outdir = Path(result.outdir)

    # Infer shape metadata so that the CLI summary in main() prints meaningful
    # numbers for logic-tree jobs (legacy dict keys are not produced here).
    n_sites = len(result.site_lons) if result.site_lons else (
        len(result.mean_rates) if result.mean_rates else 1
    )
    n_displ = len(result.d0) if result.d0 else 0

    results: Dict[str, Any] = {
        "logic_tree": True,
        "outdir": str(outdir),
        "mode": result.mode,
        "n_sites": n_sites,
        "n_displ": n_displ,
        "manifest_json": str(outdir / "manifest.json"),
        "validator_report": str(outdir / "validator_report.txt"),
    }

    if result.mode == "hazard_curve":
        results["aggregate_hazard_csv"] = str(outdir / "aggregate_hazard.csv")
    elif result.mode == "hazard_map":
        agg = outdir / "aggregate"
        results["rates_mean_h5"] = str(agg / "rates_mean.h5")
        results["rates_fractiles_h5"] = str(agg / "rates_fractiles.h5")
        results["displacement_map_mean_csv"] = str(agg / "displacement_map_mean.csv")
        if result.target_return_period is not None:
            results["target_return_period"] = float(result.target_return_period)

    if plot or plot_file:
        try:
            _plot_logic_tree_result(result, plot_file=plot_file)
        except Exception as e:  # never let plotting break the CLI
            logger.warning(f"Plot generation failed (logic-tree mode): {e}")

    if output_path:
        save_results(results, output_path)
        logger.info(f"Results saved to: {output_path}")
    return results


def _save_or_show(fig, plot_file: Optional[str]) -> None:
    """Either save the figure to ``plot_file`` or display it interactively."""
    import matplotlib.pyplot as plt

    fig.tight_layout()
    if plot_file:
        fig.savefig(plot_file, dpi=160, bbox_inches="tight")
        plt.close(fig)
        logger.info(f"Plot saved to: {plot_file}")
    else:
        plt.show()


def _plot_logic_tree_result(
    result: "LogicTreeResult",
    plot_file: Optional[str] = None,
) -> None:
    """
    Plot logic-tree outputs.

    - ``hazard_curve``: mean annual exceedance rate vs displacement for site 0,
      with an optional 16–84% fractile band when the driver produced fractiles.
    - ``hazard_map``:   scatter plot of per-site mean displacement on a
      longitude/latitude plane.

    Matplotlib-only by design, so the CLI also works on environments without
    optional GIS dependencies (e.g., cartopy).
    """
    try:
        import matplotlib.pyplot as plt  # noqa: F401  (import check only)
    except ImportError:
        logger.warning("matplotlib not available, skipping plot")
        return

    if getattr(result, "mode", "hazard_curve") == "hazard_map":
        _plot_lt_hazard_map(result, plot_file)
    else:
        _plot_lt_hazard_curve(result, plot_file)


def _plot_lt_hazard_map(result: "LogicTreeResult", plot_file: Optional[str]) -> None:
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors

    lons = np.asarray(result.site_lons or [], dtype=float)
    lats = np.asarray(result.site_lats or [], dtype=float)
    displ = np.asarray(result.displ_mean or [], dtype=float)
    if lons.size == 0 or lats.size == 0 or displ.size == 0:
        logger.warning("No map data found in logic-tree result; skipping plot")
        return

    fig, ax = plt.subplots(1, 1, figsize=(10, 8), dpi=110)
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(True, alpha=0.3, linestyle=":")

    pos = displ > 0
    if np.any(pos):
        vmin = max(1e-5, float(np.nanmin(displ[pos])))
        vmax = float(np.nanmax(displ[pos]))
        if vmax <= vmin:
            vmax = vmin * 10.0
        norm = mcolors.LogNorm(vmin=vmin, vmax=vmax)
        sc = ax.scatter(
            lons[pos], lats[pos], c=displ[pos], s=18,
            cmap="magma_r", norm=norm, edgecolors="none",
        )
        cbar = plt.colorbar(sc, ax=ax, shrink=0.85, pad=0.02)
        rp = result.target_return_period
        if rp:
            cbar.set_label(f"Mean displacement at RP={int(rp)} yr (m)")
        else:
            cbar.set_label("Mean displacement (m)")
    if np.any(~pos):
        ax.scatter(lons[~pos], lats[~pos], c="0.85", s=4, label="no hazard")

    ax.set_xlabel("Longitude (deg)")
    ax.set_ylabel("Latitude (deg)")
    ax.set_title("Logic-tree hazard map (mean displacement)")
    _save_or_show(fig, plot_file)


def _plot_lt_hazard_curve(result: "LogicTreeResult", plot_file: Optional[str]) -> None:
    import numpy as np
    import matplotlib.pyplot as plt

    d0 = np.asarray(result.d0 or [], dtype=float)
    mean_rates = np.asarray(result.mean_rates or [], dtype=float)
    if d0.size == 0 or mean_rates.size == 0:
        logger.warning("No curve data found in logic-tree result; skipping plot")
        return

    if mean_rates.ndim == 2:
        sites = mean_rates
    else:
        sites = mean_rates[np.newaxis, :]
    primary = sites[0]

    fig, ax = plt.subplots(1, 1, figsize=(7.5, 6), dpi=110)

    # Background traces for additional sites, if any.
    for i in range(1, sites.shape[0]):
        y = sites[i]
        m = y > 0
        if np.any(m):
            ax.loglog(d0[m], y[m], lw=0.6, color="0.7", alpha=0.6)

    # Optional 16–84% fractile band for site 0, when available.
    frac = getattr(result, "fractiles", None) or {}
    q_lo = np.asarray(frac.get(0.16)[0]) if 0.16 in frac else None
    q_hi = np.asarray(frac.get(0.84)[0]) if 0.84 in frac else None
    if q_lo is not None and q_hi is not None and q_lo.shape == q_hi.shape == d0.shape:
        band = (q_lo > 0) & (q_hi > 0)
        if np.any(band):
            ax.fill_between(
                d0[band], q_lo[band], q_hi[band],
                alpha=0.18, color="black", label="16–84% band (site 0)",
            )

    m = primary > 0
    if np.any(m):
        label = "LT mean" if sites.shape[0] == 1 else "LT mean (site 0)"
        ax.loglog(d0[m], primary[m], lw=2.2, color="black", label=label)

    ax.set_xlabel("Displacement (m)")
    ax.set_ylabel("Annual exceedance rate (1/yr)")
    ax.set_title("Logic-tree hazard curve")
    ax.grid(True, which="both", alpha=0.25, linestyle=":")
    ax.legend(loc="best")
    _save_or_show(fig, plot_file)


def create_parser() -> argparse.ArgumentParser:
    """Create argument parser."""
    parser = argparse.ArgumentParser(
        prog='fdha',
        description='Probabilistic Fault Displacement Hazard Analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  fdha job.ini
  fdha job.ini --plot
  fdha job.ini --plot curve.png
  fdha job.ini --output results.json
"""
    )
    
    # Positional argument for config file
    parser.add_argument(
        'job_ini',
        help='Configuration file (v5 canonical job.ini)'
    )
    
    # Output options
    parser.add_argument(
        '--output', '-o',
        help='Save results to file (JSON or CSV)'
    )
    
    parser.add_argument(
        '--plot', '-p',
        nargs='?',
        const='show',
        default=None,
        metavar='FILE',
        help='Generate plot (display if no file specified, save if file given)'
    )

    # Verbosity
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Verbose output'
    )
    
    return parser


def main():
    """Main entry point."""
    parser = create_parser()
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose)
    
    # Check config file exists
    if not os.path.exists(args.job_ini):
        print(f"Error: Configuration file not found: {args.job_ini}")
        sys.exit(1)
    
    # Handle --plot argument
    plot_show = args.plot == 'show'
    plot_file = args.plot if args.plot and args.plot != 'show' else None
    
    try:
        results = run_calculation(
            config_path=args.job_ini,
            output_path=args.output,
            plot=plot_show,
            plot_file=plot_file,
        )
        
        # Print summary
        n_sites = results.get('n_sites', len(results.get('poes', [[]])))
        n_displ = results.get('n_displ', len(results.get('imls', [])))
        print(f"\nCalculation complete:")
        print(f"  Sites: {n_sites}")
        print(f"  Displacement levels: {n_displ}")

        if results.get('logic_tree'):
            print(f"  Mode: {results.get('mode', 'hazard_curve')}")
            print(f"  Output dir: {results.get('outdir')}")

        if args.output:
            print(f"  Results: {args.output}")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
