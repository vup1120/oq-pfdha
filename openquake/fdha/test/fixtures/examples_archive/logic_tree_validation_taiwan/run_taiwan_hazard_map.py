"""Run the unified Taiwan FDHA logic-tree hazard map and plot the mean map.

One INI (``fdha_map_taiwan.ini``) drives the whole thing; the driver merges
the three LT XMLs (reverse / strike-slip / normal), enumerates end-branches
per (source, selection), subsets the source model per source_id, runs the
FDHA kernel once per end-branch, and then aggregates rates as

    total_mean(site, D0) = sum_over_sources( LT_mean_per_source(site, D0) )

The mean hazard map is the per-site inversion of that total_mean curve at
the configured return period.

Usage
-----
    cd /home/ychen/GIT/pfdha/examples/logic_tree_validation_taiwan
    python run_taiwan_hazard_map.py
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from openquake.fdha.logic_tree.driver import FdhaLogicTree

HERE = Path(__file__).resolve().parent
INI = HERE / "fdha_map_taiwan.ini"
OUT = HERE / "out_taiwan_lt"


def _trace_polyline_for_plot(source_model_xml: Path):
    """Extract (name, lons, lats) for each source's surface trace for plotting."""
    from openquake.fdha.calc.utils.parsing import parse_source_model_faults
    from openquake.fdha.calc.utils.rupture_distance import (
        _extract_fault_trace_from_mesh,
    )
    from openquake.hazardlib.geo.surface.simple_fault import SimpleFaultSurface

    faults = parse_source_model_faults(
        str(source_model_xml), hdf5path="",
        rupture_mesh_spacing=2.0, width_of_mfd_bin=0.1,
        complex_fault_mesh_spacing=2.0,
    )
    out = []
    for sid, src in faults.items():
        surf = None
        if hasattr(src, "surface"):
            surf = src.surface
        elif all(hasattr(src, a) for a in (
            "fault_trace", "upper_seismogenic_depth",
            "lower_seismogenic_depth", "dip",
        )):
            try:
                surf = SimpleFaultSurface.from_fault_data(
                    src.fault_trace, src.upper_seismogenic_depth,
                    src.lower_seismogenic_depth, src.dip, 2.0,
                )
            except Exception:
                surf = None
        if surf is not None:
            poly = _extract_fault_trace_from_mesh(surf)
            lons = [p[0] for p in poly]
            lats = [p[1] for p in poly]
        elif hasattr(src, "fault_trace"):
            lons = [p.longitude for p in src.fault_trace]
            lats = [p.latitude for p in src.fault_trace]
        else:
            continue
        out.append((sid, getattr(src, "name", sid), lons, lats))
    return out


def plot_hazard_map(result, outpath: Path, source_xml: Path, title: str):
    lons = np.asarray(result.site_lons, dtype=float)
    lats = np.asarray(result.site_lats, dtype=float)
    displ = np.asarray(result.displ_mean, dtype=float)

    fig, ax = plt.subplots(1, 1, figsize=(9, 11), dpi=110)

    # Full bounding box context.
    ax.set_xlim(119.5, 122.5)
    ax.set_ylim(21.5, 25.7)
    ax.set_aspect(1.0 / np.cos(np.deg2rad(23.5)))
    ax.grid(True, alpha=0.3, linestyle=":")

    # Fault traces first for layering.
    traces = _trace_polyline_for_plot(source_xml)
    for sid, name, tlons, tlats in traces:
        ax.plot(tlons, tlats, "-", color="0.25", lw=1.6, zorder=2)
        if tlons:
            ax.annotate(
                f"{name}",
                xy=(tlons[len(tlons) // 2], tlats[len(tlats) // 2]),
                xytext=(6, 6), textcoords="offset points",
                fontsize=8, color="0.2",
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                          edgecolor="0.5", alpha=0.85),
                zorder=5,
            )

    # Sites coloured by mean displacement.
    # Log scale is useful for PFDH where values span many orders of magnitude.
    pos = displ > 0
    sc = ax.scatter(
        lons[pos], lats[pos], c=displ[pos],
        cmap="magma_r", s=18,
        norm=__import__("matplotlib").colors.LogNorm(
            vmin=max(1e-3, float(np.nanmin(displ[pos]))),
            vmax=float(np.nanmax(displ[pos])),
        ),
        edgecolors="none", zorder=3,
    )
    # Show zero-displ sites as small grey dots so the grid footprint is visible.
    ax.scatter(lons[~pos], lats[~pos], c="0.85", s=4, zorder=1, label="no hazard")

    cbar = plt.colorbar(sc, ax=ax, shrink=0.72, pad=0.02)
    cbar.set_label(f"Mean displacement at RP={int(result.target_return_period)} yr  [m]")

    ax.set_xlabel("Longitude (deg)")
    ax.set_ylabel("Latitude (deg)")
    ax.set_title(title, fontsize=12)

    fig.tight_layout()
    fig.savefig(outpath, dpi=160)
    plt.close(fig)


def main() -> None:
    # Quiet the kernel: we monitor only top-level LT progress.
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")
    # Silence the per-rupture tqdm pulse from calculate_fdha_hazard.
    import os as _os
    _os.environ.setdefault("TQDM_DISABLE", "1")
    # Let the driver emit high-level progress at INFO level.
    logging.getLogger("openquake.fdha.logic_tree.driver").setLevel(logging.INFO)
    OUT.mkdir(parents=True, exist_ok=True)

    print(f"Running LT: {INI}", flush=True)
    t0 = time.time()
    lt = FdhaLogicTree.from_ini(INI)
    result = lt.run(outdir=OUT)
    dt = time.time() - t0
    print(f"Done in {dt:.1f} s. mode={result.mode}, "
          f"n_sites={len(result.site_lons)}, n_d0={len(result.d0)}", flush=True)

    manifest = json.loads((OUT / "manifest.json").read_text())
    n_br = len(manifest["branches"])
    per_src = manifest.get("per_source_ids", [])
    print(f"Manifest: {n_br} end-branches; per_source_ids={per_src}")

    displ = np.asarray(result.displ_mean, dtype=float)
    print(
        f"displ_mean  min={displ[displ>0].min():.4f} m  "
        f"max={displ.max():.4f} m  "
        f"n_nonzero={int((displ>0).sum())}/{displ.size}"
    )

    source_xml = HERE / manifest.get("grid", {}).get("source_model", "taiwan_representative_3faults.xml")
    if not source_xml.exists():
        source_xml = HERE / "taiwan_representative_3faults.xml"

    png = HERE / "taiwan_hazard_map_mean.png"
    plot_hazard_map(
        result, png, source_xml,
        title=(
            f"Taiwan PFDH logic-tree mean displacement map\n"
            f"RP = {int(result.target_return_period)} yr   "
            f"(3 representative faults, {n_br} end-branches, full LT)"
        ),
    )
    print(f"Wrote {png}")
    print(f"Wrote {OUT / 'aggregate' / 'rates_mean.h5'}")
    print(f"Wrote {OUT / 'aggregate' / 'rates_fractiles.h5'}")
    print(f"Wrote {OUT / 'aggregate' / 'displacement_map_mean.csv'}")
    print(f"Wrote {OUT / 'manifest.json'}")


if __name__ == "__main__":
    main()
