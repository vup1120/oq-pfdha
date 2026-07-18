"""Plot hazard maps from the Logic Tree Validation demo (map mode).

Mirror of ``../plot_hazard_curves.py`` but for hazard maps. Overlays, in a
single figure:

  1. single_bilinear_map       - weight-1.0 tree on Petersen2011PrimaryFD (version=bilinear)
  2. single_elliptical_map     - weight-1.0 tree on Petersen2011PrimaryFD (version=elliptical)
  3. logic-tree blend          - displacement_map_mean.csv from blend_50_50_map
  4. analytical post-processing - invert  0.5 * rates_A + 0.5 * rates_B
                                  at the same return period

Panels (3) and (4) must coincide to numerical noise because aggregation
happens in rate space and inversion is deterministic; ``verify_map.py``
already proves that to 1e-10 m.

One fault, one return period, one colorbar shared across all panels.
The strike-slip source trace from ``../source_model_ss.xml`` is overlaid
on every panel so the reader can see which sites are on-trace vs off-trace.
"""
from __future__ import annotations

import csv
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

from openquake.fdha.calc.utils.interpolation import get_map_from_curves

HERE = Path(__file__).resolve().parent
RUNS = {
    "single_bilinear_map": HERE / "single_bilinear_map",
    "single_elliptical_map": HERE / "single_elliptical_map",
    "blend_50_50_map": HERE / "blend_50_50_map",
}
SOURCE_XML = HERE.parent / "source_model_ss.xml"
OUT_PNG = HERE / "hazard_maps_comparison.png"


# ------------------------------------------------------------- readers
def _read_rates_mean(path: Path):
    with h5py.File(path, "r") as f:
        rates = np.asarray(f["rates_mean"][()], dtype=float)
        d0 = np.asarray(f["d0"][()], dtype=float)
        lons = np.asarray(f["site_lons"][()], dtype=float)
        lats = np.asarray(f["site_lats"][()], dtype=float)
    return rates, d0, lons, lats


def _read_displ_csv(path: Path):
    rp = None
    rows: list[tuple[int, float, float, int, float]] = []
    with path.open() as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith("#") and "return_period" in line:
            rp = float(line.split("=")[1].strip())
            break
    reader = csv.reader((ln for ln in lines if not ln.startswith("#")))
    next(reader)  # header
    for row in reader:
        rows.append(
            (int(row[0]), float(row[1]), float(row[2]), int(row[3]), float(row[4]))
        )
    return rp, rows


def _read_fault_trace(xml_path: Path) -> tuple[list[float], list[float]]:
    """Pull the raw <gml:posList> polyline; good enough for plotting."""
    import re
    text = xml_path.read_text()
    m = re.search(r"<gml:posList>(.*?)</gml:posList>", text, re.DOTALL)
    if not m:
        return [], []
    toks = m.group(1).split()
    lons = [float(toks[i]) for i in range(0, len(toks), 2)]
    lats = [float(toks[i + 1]) for i in range(0, len(toks), 2)]
    return lons, lats


# ------------------------------------------------------------- main
def main() -> None:
    # 1) Individual model runs (the "component" step).
    rA, d0, lonsA, latsA = _read_rates_mean(
        RUNS["single_bilinear_map"] / "out" / "aggregate" / "rates_mean.h5"
    )
    rB, d0B, lonsB, latsB = _read_rates_mean(
        RUNS["single_elliptical_map"] / "out" / "aggregate" / "rates_mean.h5"
    )
    rC_lt, d0C, lonsC, latsC = _read_rates_mean(
        RUNS["blend_50_50_map"] / "out" / "aggregate" / "rates_mean.h5"
    )

    # Pre-condition: grid consistency (verify_map.py also checks this).
    if not (np.allclose(d0, d0B) and np.allclose(d0, d0C)
            and np.allclose(lonsA, lonsB) and np.allclose(lonsA, lonsC)
            and np.allclose(latsA, latsB) and np.allclose(latsA, latsC)):
        raise SystemExit(
            "Grid / D0 mismatch between runs; rerun verify_map.py first."
        )

    # 2) Post-processing weighted average (rate space).
    rC_manual = 0.5 * rA + 0.5 * rB
    max_rate_diff = float(np.max(np.abs(rC_lt - rC_manual)))

    # 3) LT's own displacement map (what the framework outputs).
    rp_lt, rows_lt = _read_displ_csv(
        RUNS["blend_50_50_map"] / "out" / "aggregate" / "displacement_map_mean.csv"
    )
    displ_lt = np.array([r[4] for r in rows_lt])
    # A & B displacement maps for the top row.
    _, rows_A = _read_displ_csv(
        RUNS["single_bilinear_map"] / "out" / "aggregate" / "displacement_map_mean.csv"
    )
    displ_A = np.array([r[4] for r in rows_A])
    _, rows_B = _read_displ_csv(
        RUNS["single_elliptical_map"] / "out" / "aggregate" / "displacement_map_mean.csv"
    )
    displ_B = np.array([r[4] for r in rows_B])

    # 4) Analytical post-processing in displacement space: invert the
    #    manual blend rate curve at the same return period.
    assert rp_lt is not None
    pex = 1.0 / rp_lt
    displ_manual = get_map_from_curves(d0, rC_manual, pex)
    max_displ_diff = float(np.max(np.abs(displ_lt - displ_manual)))

    # ------------------------------------------------------------- plot
    lons, lats = lonsA, latsA
    trace_lons, trace_lats = _read_fault_trace(SOURCE_XML)

    # Shared log-scale colorbar over every panel.
    stacked = np.concatenate([displ_A, displ_B, displ_lt, displ_manual])
    pos = stacked > 0
    vmin = max(1e-4, float(np.min(stacked[pos])))
    vmax = float(np.max(stacked[pos]))
    norm = plt.matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)

    fig, axes = plt.subplots(
        2, 2, figsize=(11.5, 10.2), dpi=130,
        constrained_layout=True,
    )
    fig.suptitle(
        "Logic Tree Validation demo - hazard maps\n"
        "Petersen2011PrimaryFD bilinear vs elliptical, 50/50 blend\n"
        f"(strike-slip source, one-fault demo; RP = {int(rp_lt)} yr)",
        fontsize=12,
    )

    panels = [
        (axes[0, 0], displ_A,
         "(a) single_bilinear_map  (weight = 1.0)",
         "tab:blue"),
        (axes[0, 1], displ_B,
         "(b) single_elliptical_map (weight = 1.0)",
         "tab:orange"),
        (axes[1, 0], displ_lt,
         "(c) LT framework  blend_50_50.mean",
         "tab:red"),
        (axes[1, 1], displ_manual,
         r"(d) analytical  invert$\,[0.5\cdot\lambda_A + 0.5\cdot\lambda_B\,]$",
         "tab:green"),
    ]

    for ax, displ, title, _color in panels:
        # Fault trace.
        if trace_lons:
            ax.plot(trace_lons, trace_lats, "-", color="0.15", lw=1.8,
                    zorder=3, label="fault trace")
        pos_mask = displ > 0
        # Non-zero sites on the log colormap.
        sc = ax.scatter(
            lons[pos_mask], lats[pos_mask], c=displ[pos_mask],
            cmap="magma_r", s=42, norm=norm,
            edgecolors="0.3", linewidths=0.3, zorder=4,
        )
        # Zero-displ sites (off-fault) as tiny grey dots.
        ax.scatter(lons[~pos_mask], lats[~pos_mask],
                   c="0.85", s=10, zorder=2)
        ax.set_title(title, fontsize=10)
        ax.set_aspect(1.0 / np.cos(np.deg2rad(float(lats.mean()))))
        ax.grid(True, ls=":", alpha=0.5)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

    cbar = fig.colorbar(
        sc, ax=axes, shrink=0.9, pad=0.02,
        location="right", aspect=35,
    )
    cbar.set_label(f"Mean displacement at RP = {int(rp_lt)} yr  [m]")

    # On-figure annotation with the two numerical checks (printed inside
    # panel (d) so constrained_layout does not have to reserve space).
    txt = (
        r"$\max\,|\,\mathrm{LT\ blend} - 0.5\cdot\lambda_A - 0.5\cdot\lambda_B\,|$"
        f"\n(rate space) = {max_rate_diff:.2e}  (tol $10^{{-12}}$)"
        "\n"
        r"$\max\,|\,\mathrm{LT\ displ} - \mathrm{analytical\ displ}\,|$"
        f"\n(displacement space) = {max_displ_diff:.2e}  (tol $10^{{-10}}$)"
    )
    axes[1, 1].text(
        0.02, 0.98, txt, transform=axes[1, 1].transAxes,
        ha="left", va="top", fontsize=8,
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.6", alpha=0.95),
        zorder=10,
    )

    fig.savefig(OUT_PNG, dpi=160, bbox_inches="tight")
    print(f"Wrote {OUT_PNG.relative_to(HERE.parent)}")
    print(f"max |LT blend - analytical| (rate space)        = {max_rate_diff:.3e}  1/yr")
    print(f"max |LT blend - analytical| (displacement space) = {max_displ_diff:.3e}  m")


if __name__ == "__main__":
    main()
