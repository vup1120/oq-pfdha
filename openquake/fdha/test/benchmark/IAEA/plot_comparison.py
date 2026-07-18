#!/usr/bin/env python
"""Plot the IAEA exercise benchmark: oq-pfdha curves vs the published curves.

Reads the ``computed/<case>_<job>.csv`` files written by ``run_all.py`` and
overlays them (dashed) on the digitized paper curves (solid), reproducing the
layout of Figs 4 and 6 of the Valentini et al. IAEA exercise paper.

Usage (repo root)::

    PYTHONPATH=. python openquake/fdha/test/benchmark/IAEA/plot_comparison.py

Writes ``Figures/iaea_fig4_principal_comparison.png`` and
``Figures/iaea_fig6_distributed_comparison.png`` next to this script.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from manifest import MANIFEST  # noqa: E402

COLORS = {"P11": "#d62728", "T13": "#000000", "L23": "#c724c7",
          "C24": "#2ca02c", "K24": "#1f77b4", "M11": "#17becf",
          "Y03": "#d3a020", "V24": "#77ac30"}

PANELS_FIG4 = [
    ("fig4a_kumamoto_principal.csv", "kumamoto", "principal",
     "(a) Kumamoto - principal (single segment)"),
    ("fig4b_kumamoto_principal_floating.csv", "kumamoto", "floating",
     "(b) Kumamoto - principal (floating)"),
    ("fig4c_leteil_principal.csv", "le_teil", "principal",
     "(c) Le Teil - principal"),
    ("fig4d_norcia_principal.csv", "norcia", "principal",
     "(d) Norcia - principal"),
]
PANELS_FIG6 = [
    ("fig6a_kumamoto_distributed.csv", "kumamoto", "distributed",
     "(a) Kumamoto - distributed (r=5.2 km)"),
    ("fig6b_leteil_distributed.csv", "le_teil", "distributed",
     "(b) Le Teil - distributed (r=0.6 km)"),
    ("fig6c_norcia_distributed.csv", "norcia", "distributed",
     "(c) Norcia - distributed (r=7.6 km)"),
]


def read_csv_columns(path: Path):
    with path.open() as f:
        rows = list(csv.DictReader(f))
    cols = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
    return cols


def plot_panel(ax, figure_csv: str, case: str, job_prefix: str, title: str,
               ylim):
    ref = read_csv_columns(HERE / "reference" / figure_csv)
    disp_cm = ref.pop("disp_cm")
    for model, vals in ref.items():
        ax.loglog(disp_cm, vals, color=COLORS.get(model, "grey"), lw=2.2,
                  alpha=0.9, label=f"{model} (paper)")

    for entry in MANIFEST:
        if entry.case != case or not entry.job.startswith(job_prefix):
            continue
        # Skip the base-case-tree diagnostic variants: they share a column
        # (e.g. M11, T13) and figure with the single-scenario team chain,
        # but do NOT correspond to the published curve (see README), so
        # plotting them here would draw a second, identically labelled
        # dashed curve. Their own comparison is documented separately.
        if entry.job.endswith("_basecase"):
            continue
        if entry.figure_csv != figure_csv and not (
                figure_csv.startswith("fig4c") and
                entry.figure_csv.startswith("fig4c")):
            continue
        comp_path = HERE / "computed" / f"{entry.case}_{entry.job}.csv"
        if not comp_path.exists():
            continue
        comp = read_csv_columns(comp_path)
        pos = comp["rate_post_factor"] > 0
        ax.loglog(comp["disp_m"][pos] * 100.0, comp["rate_post_factor"][pos],
                  "--", color=COLORS.get(entry.column, "grey"), lw=1.8,
                  label=f"{entry.column} (oq-pfdha)")

    ax.set_xlim(1, 1000)
    ax.set_ylim(*ylim)
    ax.set_xlabel("Displacement (cm)")
    ax.set_ylabel("AFOE (yr$^{-1}$)")
    ax.set_title(title, fontsize=10)
    ax.grid(True, which="both", ls=":", alpha=0.35)
    ax.legend(fontsize=6.5, ncol=2, loc="lower left")


def main() -> None:
    figures_dir = HERE / "Figures"
    figures_dir.mkdir(exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(11, 10))
    for ax, (csv_name, case, prefix, title) in zip(axes.flat, PANELS_FIG4):
        plot_panel(ax, csv_name, case, prefix, title, ylim=(1e-8, 1e-3))
        if csv_name.startswith("fig4c"):
            # overlay the M11 curve, which lives in its own reference CSV
            m11 = read_csv_columns(HERE / "reference" /
                                   "fig4c_leteil_principal_M11.csv")
            ax.loglog(m11["disp_cm"], m11["M11"], color=COLORS["M11"], lw=2.2,
                      alpha=0.9, label="M11 (paper)")
            ax.legend(fontsize=6.5, ncol=2, loc="lower left")
    fig.suptitle("IAEA PFDHA exercise - principal fault displacement "
                 "(paper Fig. 4 vs oq-pfdha)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = figures_dir / "iaea_fig4_principal_comparison.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print("saved", out)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, (csv_name, case, prefix, title) in zip(axes.flat, PANELS_FIG6):
        # paper axis is 1e-10..1e-5; extend down one decade so the Le Teil
        # V24 curve computed with the final (TECDOC) Visini model stays visible
        plot_panel(ax, csv_name, case, prefix, title, ylim=(1e-11, 1e-5))
    fig.suptitle("IAEA PFDHA exercise - distributed fault displacement "
                 "(paper Fig. 6 vs oq-pfdha)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = figures_dir / "iaea_fig6_distributed_comparison.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print("saved", out)


if __name__ == "__main__":
    main()
