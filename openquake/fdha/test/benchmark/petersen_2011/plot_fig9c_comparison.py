#!/usr/bin/env python
"""Overlay the computed across-strike hazard profile on Petersen et al.
(2011) Fig. 9c.

The curves are produced by the very functions the pytest module asserts on,
so the figure cannot drift from the test. Two constructions are shown per
mapping-accuracy class: the tool's stated rupture-location weight, and the
documented figure-equivalent weight 0.90*exp(-r^2 / 2(1.65 sigma)^2) that
reproduces the published curve (see README, "Reading the offset").

Usage (repo root)::

    PYTHONPATH=. python openquake/fdha/test/benchmark/petersen_2011/plot_fig9c_comparison.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from test_petersen2011_fig9c import (            # noqa: E402
    SIGMA_KM, _R_GRID_M, stated_method_curve, figure_offset_curve)

COLOURS = {"Accurate": "tab:blue", "Approximate": "tab:orange",
           "Concealed": "tab:green", "Inferred": "tab:red"}


def main():
    dig = np.load(HERE / "fig9c_digitized.npz")
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.0), sharey=True)

    for ax, (title, fn) in zip(axes, (
            ("(a) tool's stated rupture-location weight", stated_method_curve),
            ("(b) documented figure-equivalent weight", figure_offset_curve))):
        for name, sigma in SIGMA_KM.items():
            r_dig, d_dig = dig[name]
            ax.plot(r_dig, d_dig, ".", ms=2.2, alpha=0.45,
                    color=COLOURS[name])
            curve = fn(_R_GRID_M, sigma)
            ax.plot(_R_GRID_M, curve, "-", lw=1.7, color=COLOURS[name],
                    label=f"{name} ($\\sigma$ = {sigma*1000:.0f} m)")
            ax.plot(-_R_GRID_M, curve, "-", lw=1.7, color=COLOURS[name])
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Across-strike distance from the mapped trace, r (m)")
        ax.set_xlim(-300, 300)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("Displacement at 10 % in 50 yr (cm)")
    axes[0].legend(fontsize=8, loc="upper right")
    fig.suptitle("Petersen et al. (2011) Fig. 9c — computed profile "
                 "(lines) vs digitized published curves (points)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    outdir = HERE / "Figures"
    outdir.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(outdir / f"petersen2011_fig9c_comparison.{ext}", dpi=170)
    print("wrote", outdir / "petersen2011_fig9c_comparison.png")


if __name__ == "__main__":
    sys.exit(main())
