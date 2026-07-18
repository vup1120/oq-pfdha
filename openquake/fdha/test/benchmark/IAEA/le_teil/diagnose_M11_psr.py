#!/usr/bin/env python
"""Decompose the Le Teil M11 plateau gap: it is entirely the P_sr term.

The published M11 curve (Fig 4c) plateaus ~13x below our M11 chain. This
script shows that the *plateau* difference is 100% a surface-rupture-
probability (P_sr) effect, not a displacement-distribution effect:

- both curves are flat at small displacement (P(D > 0.01 m) ~ 1), so the
  plateau equals rate x P_sr;
- rate = 4.6e-5 is common to all teams (K24, with P_sr = 1, plateaus at
  the rate);
- scaling our curve down to the paper's plateau makes the two overlay
  exactly up to ~10 cm (ratio 1.00-1.06), i.e. identical shape there -
  the gap is a pure vertical scale = a P_sr ratio.

Attribution of the ~13x plateau gap:

    our P_sr (Moss 2013 soft, vs30 200) = 15.6%
    our P_sr (Moss 2013 stiff, vs30 760) = 9.9%   <- matches the paper's
                                                     stated "10% at Mw5.5"
    paper's plotted effective P_sr        = 1.22%

    soft-vs-stiff soil class explains only 1.58x;
    the residual 8.1x (stiff 9.9% vs plotted 1.2%) is NOT explained by the
    P_sr model - our Moss 2013 reproduces the value the paper reports; the
    paper's *plotted* curve uses an effective P_sr 8x below its own text.

(Beyond ~30 cm the published curve also has a fatter tail - a separate
magnitude-scaling issue, see diagnose_M11.py: its normalized shape implies
a median displacement at Mw ~6.6, not 5.5.)

Writes ``../Figures/iaea_leteil_M11_psr_decomposition.png``.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from openquake.fdha.primary_surf_rup import Moss2013PrimarySR
from openquake.fdha.primary_surf_displ import MossRoss2011PrimaryFD

HERE = Path(__file__).resolve().parent
RATE = 4.6e-5
XL = 0.46
MAG = 5.5


def m11_curve(vs30, d):
    """rate x P_sr(vs30) x P(D>d | Moss&Ross AD, M5.5, x/L)."""
    psr = float(Moss2013PrimarySR().get_prob(MAG, "reverse", vs30))
    fd = MossRoss2011PrimaryFD(n_sigma=5.0)
    pexc = np.atleast_1d(np.asarray(fd.get_prob(d, XL, MAG, "AD")).squeeze())
    return psr, RATE * psr * pexc


def main():
    ref = list(csv.DictReader(
        (HERE.parent / "reference" / "fig4c_leteil_principal_M11.csv").open()))
    dp = np.array([float(r["disp_cm"]) for r in ref]) / 100.0
    rp = np.array([float(r["M11"]) for r in ref])
    eff_psr = rp[0] / RATE

    d = np.logspace(-4, 1.0, 200)
    psr_soft, c_soft = m11_curve(200.0, d)
    psr_stiff, c_stiff = m11_curve(760.0, d)

    fig, ax = plt.subplots(figsize=(8.5, 6))
    ax.plot(dp * 100, rp, "-", color="k", lw=2.4,
            label=f"paper M11 (plateau ⇒ P_sr = {eff_psr*100:.1f}%)")
    ax.plot(d * 100, c_soft, "--", color="tab:red", lw=1.8,
            label=f"ours, Moss2013 soft vs30=200 (P_sr = {psr_soft*100:.1f}%)")
    ax.plot(d * 100, c_stiff, "--", color="tab:blue", lw=1.8,
            label=f"ours, Moss2013 stiff vs30=760 (P_sr = {psr_stiff*100:.1f}%)")
    # our stiff curve rescaled to the paper plateau: overlays through ~10 cm
    ax.plot(d * 100, c_stiff * eff_psr / psr_stiff, ":", color="0.5", lw=1.6,
            label="ours(stiff) × (1.2%/9.9%) - pure P_sr rescale")

    # P_sr reference lines as horizontal plateaus (rate x P_sr)
    for psr, c in [(0.10, "tab:green")]:
        ax.axhline(RATE * psr, color=c, ls="-.", lw=1,
                   label=f"rate × paper-text 10% = {RATE*psr:.1e}")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-1, 1e3)
    ax.set_ylim(1e-9, 1e-4)
    ax.set_xlabel("Displacement (cm)")
    ax.set_ylabel("Annual rate of exceedance")
    ax.set_title("Le Teil M11 - the plateau gap is the P_sr term\n"
                 "(identical shape ≤10 cm; fatter paper tail >30 cm is a "
                 "separate magnitude issue)")
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    out = HERE.parent / "Figures" / "iaea_leteil_M11_psr_decomposition.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"paper effective P_sr = {eff_psr*100:.2f}%  "
          f"(stiff {psr_stiff*100:.1f}%, soft {psr_soft*100:.1f}%)")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
