#!/usr/bin/env python
"""Overlay Moss2013PrimarySR on both digitisations of the published curves.

Writes ``Figures/moss2013_fig3_comparison.png``: the library model (lines)
against the curves digitised from Moss et al. (2013) Fig. 3 (filled
markers) and from the GIRS-2022-05 Fig. 3.2 redraw (open markers), for the
reverse and strike-slip panels.
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from openquake.fdha.primary_surf_rup import Moss2013PrimarySR

HERE = os.path.dirname(os.path.abspath(__file__))
VS30 = {"stiff": 760.0, "soft": 200.0}
COLOURS = {"stiff": "#7a25a6", "soft": "#6e9529"}


def load(name):
    with open(os.path.join(HERE, "reference", name)) as f:
        rows = list(csv.DictReader(f))
    out = {}
    for r in rows:
        out.setdefault((r["panel"], r["curve"]), []).append(
            (float(r["mw"]), float(r["p"])))
    return {k: np.array(v).T for k, v in out.items()}


def main():
    original = load("moss2013_fig3.csv")
    girs = load("girs_fig3_2.csv")
    model = Moss2013PrimarySR()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    for ax, panel, xlim in zip(axes, ("reverse", "strike-slip"),
                               ((5, 9), (4, 9))):
        for curve in ("stiff", "soft"):
            c = COLOURS[curve]
            mw, p = original[(panel, curve)]
            ax.plot(mw, p, "o", ms=2.5, color=c, alpha=0.6,
                    label=f"Moss et al. (2013) Fig. 3, {curve}")
            mw, p = girs[(panel, curve)]
            ax.plot(mw, p, "o", ms=3, mfc="none", mec=c, mew=0.6, alpha=0.5,
                    label=f"GIRS-2022-05 Fig. 3.2 redraw, {curve}")
            m = np.linspace(4.2, 8.7, 200)
            ax.plot(m, model.get_prob(m, style=panel, vs30=VS30[curve]),
                    "-", lw=1.2, color=c, label=f"Moss2013PrimarySR, {curve}")
        ax.set_title(f"{panel.capitalize()} faulting")
        ax.set_xlim(*xlim)
        ax.set_ylim(0, 1)
        ax.set_xlabel("Mw")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=6.5, loc="upper left")
    axes[0].set_ylabel("P(sr ≠ 0 | Mw)")
    fig.tight_layout()
    os.makedirs(os.path.join(HERE, "Figures"), exist_ok=True)
    out = os.path.join(HERE, "Figures", "moss2013_fig3_comparison.png")
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
