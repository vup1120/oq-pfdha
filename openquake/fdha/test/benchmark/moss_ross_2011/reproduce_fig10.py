#!/usr/bin/env python
"""Reproduce the Moss & Ross (2011) Los Osos example (their Figs 7 and 10).

Moss & Ross (2011, BSSA 101(4), 1542-1553) illustrate their reverse-fault
PFDHA on the Los Osos fault zone: slip rate 0.5 mm/yr, fault area
44 km x 14 km, shear modulus 3.75e11 dyne/cm2, b = 0.8, Mmax = Mw 7,
x/L = 0.25, AD method with the gamma D/AD distribution (their Eq. 7), and
events below Mw 5 neglected. Their hazard integral (Eqs. 2-3) evaluates
every event at the fixed site position x/L = 0.25 -- no rupture floating --
so this benchmark reproduces it as the direct integral

    nu(d) = alpha * int f(m) * P_sr(m) * P(D>d | m, x/L) dm ,

with f(m) a truncated-exponential (Gutenberg-Richter) density on
Mw 5.0-7.0 and alpha moment-balanced through log M0 = 1.5 Mw + 16.05.
The reference curve (their Fig. 7) uses their reverse P_sr (Eq. 5 =
``MossRoss2011PrimarySR``); the Fig. 10 comparison swaps in the
Wells & Coppersmith (1993) all-slip-types model (``WC1993PrimarySR``).

Published anchors used for the comparison (paper text):

  reference (reverse):  2%-in-50-yr (1/2475 yr)  -> 55 cm
                        1%-in-50-yr (1/4975 yr)  -> 105 cm
  Fig. 10 (all types):  2%-in-50-yr              -> 95 cm
                        1%-in-50-yr              -> 142 cm
  plateau shift: all-slip-types plateau ~ +45% above the reverse plateau

Writes ``Figures/mossross2011_fig10_reproduction.png`` and
``fig10_agreement.json`` next to this script.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("NUMBA_DISABLE_JIT", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from openquake.fdha.primary_surf_displ import MossRoss2011PrimaryFD
from openquake.fdha.primary_surf_rup import MossRoss2011PrimarySR, WC1993PrimarySR

HERE = Path(__file__).resolve().parent

# ---- Los Osos source characterization (paper, Example section)
MU = 3.75e11                # shear modulus, dyne/cm^2
AREA_CM2 = 44e5 * 14e5      # 44 km x 14 km in cm^2
SLIP_CM_YR = 0.05           # 0.5 mm/yr
B_VALUE, M_MIN, M_MAX = 0.8, 5.0, 7.0
HK_CONST = 16.05            # log M0 = 1.5 Mw + 16.05 (dyne-cm)
X_L = 0.25

NU_2PC_50YR = 1.0 / 2475.0
NU_1PC_50YR = 1.0 / 4975.0

ANCHORS = {  # (hazard level, published displacement in m)
    "reverse": [(NU_2PC_50YR, 0.55), (NU_1PC_50YR, 1.05)],
    "all_slip_types": [(NU_2PC_50YR, 0.95), (NU_1PC_50YR, 1.42)],
}


def hazard_curves():
    m = np.arange(M_MIN + 0.005, M_MAX, 0.01)
    beta = B_VALUE * np.log(10.0)
    pdf = beta * np.exp(-beta * (m - M_MIN))
    pdf /= np.trapezoid(pdf, m)

    moment_rate = MU * AREA_CM2 * SLIP_CM_YR
    alpha = moment_rate / np.trapezoid(pdf * 10 ** (1.5 * m + HK_CONST), m)

    d = np.logspace(-2, 1.02, 90)  # 0.01 - 10.5 m
    fd = MossRoss2011PrimaryFD(n_sigma=5.0)  # paper's MC is untruncated
    p_exceed = np.stack(
        [np.asarray(fd.get_prob(d=d, X_L_ratio=X_L, mag=float(mm),
                                norm_disp_type="AD")).reshape(len(d), -1)[:, 0]
         for mm in m],
        axis=1,
    )
    psr_rev = np.array([float(MossRoss2011PrimarySR().get_prob(mag=float(mm)))
                        for mm in m])
    psr_all = np.array([WC1993PrimarySR().get_prob(mag=float(mm)) for mm in m])

    nu_rev = np.trapezoid(alpha * pdf * psr_rev * p_exceed, m, axis=1)
    nu_all = np.trapezoid(alpha * pdf * psr_all * p_exceed, m, axis=1)
    return d, nu_rev, nu_all, alpha


def displacement_at(d, nu, target):
    """Displacement (m) at hazard level ``target`` by log-log interpolation."""
    return float(np.exp(np.interp(np.log(target), np.log(nu[::-1]),
                                  np.log(d[::-1]))))


def main() -> None:
    d, nu_rev, nu_all, alpha = hazard_curves()
    print(f"alpha(M>=5) = {alpha:.4e} /yr")

    agreement = {"plateau_ratio": round(float(nu_all[0] / nu_rev[0]), 3),
                 "plateau_ratio_paper": "~1.45 ('increased by nearly 45%')"}
    for name, nu in (("reverse", nu_rev), ("all_slip_types", nu_all)):
        for level, d_paper in ANCHORS[name]:
            ours = displacement_at(d, nu, level)
            key = f"{name}_d_at_{'2' if level == NU_2PC_50YR else '1'}pc50yr"
            agreement[key] = {
                "computed_cm": round(ours * 100, 1),
                "paper_cm": d_paper * 100,
                "ratio": round(ours / d_paper, 3),
            }
            print(f"{name:15s} nu={level:.3e}: computed {ours*100:6.1f} cm, "
                  f"paper {d_paper*100:.0f} cm, ratio {ours/d_paper:.3f}")
    print(f"plateau ratio all/reverse = {agreement['plateau_ratio']}")

    figures = HERE / "Figures"
    figures.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.5, 6))
    ax.loglog(d, nu_rev, "-", color="#d62728", lw=2.2,
              label="reverse P$_{sr}$ (Eq. 5) — reference curve (paper Fig. 7)")
    ax.loglog(d, nu_all, "-", color="#1f77b4", lw=2.2,
              label="all-slip-types P$_{sr}$ (WC93) — comparison (paper Fig. 10)")
    for name, color in (("reverse", "#d62728"), ("all_slip_types", "#1f77b4")):
        pts = ANCHORS[name]
        ax.loglog([p[1] for p in pts], [p[0] for p in pts], "X",
                  color=color, ms=11, mec="k", mew=0.8,
                  label=f"paper anchors ({name.replace('_', ' ')})")
    for level, lbl in ((NU_2PC_50YR, "2% in 50 yr"), (NU_1PC_50YR, "1% in 50 yr")):
        ax.axhline(level, color="grey", ls=":", lw=1)
        ax.text(0.011, level * 1.1, lbl, fontsize=8, color="grey")
    ax.set_xlabel("Displacement (m)")
    ax.set_ylabel("Annual probability of exceedance")
    ax.set_xlim(0.01, 10)
    ax.set_ylim(1e-7, 1e-2)
    ax.set_title("Moss & Ross (2011) Los Osos example — oq-pfdha reproduction\n"
                 "of the Fig. 10 surface-rupture-distribution comparison")
    ax.grid(True, which="both", ls=":", alpha=0.35)
    ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    out = figures / "mossross2011_fig10_reproduction.png"
    fig.savefig(out, dpi=150)
    print("saved", out)

    # All-slip-types-only variant (no reverse curve/anchors)
    fig, ax = plt.subplots(figsize=(8.5, 6))
    ax.loglog(d, nu_all, "-", color="#1f77b4", lw=2.2,
              label="all-slip-types P$_{sr}$ (WC93) — oq-pfdha")
    pts = ANCHORS["all_slip_types"]
    ax.loglog([p[1] for p in pts], [p[0] for p in pts], "X",
              color="#1f77b4", ms=11, mec="k", mew=0.8,
              label="paper anchors (all slip types)")
    for level, lbl in ((NU_2PC_50YR, "2% in 50 yr"), (NU_1PC_50YR, "1% in 50 yr")):
        ax.axhline(level, color="grey", ls=":", lw=1)
        ax.text(0.011, level * 1.1, lbl, fontsize=8, color="grey")
    ax.set_xlabel("Displacement (m)")
    ax.set_ylabel("Annual probability of exceedance")
    ax.set_xlim(0.01, 10)
    ax.set_ylim(1e-7, 1e-2)
    ax.set_title("Moss & Ross (2011) Los Osos example — all-slip-types "
                 "P$_{sr}$ curve\n(paper Fig. 10 comparison curve) vs oq-pfdha")
    ax.grid(True, which="both", ls=":", alpha=0.35)
    ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    out = figures / "mossross2011_fig10_reproduction_allslip.png"
    fig.savefig(out, dpi=150)
    print("saved", out)

    (HERE / "fig10_agreement.json").write_text(
        json.dumps(agreement, indent=2) + "\n")
    print("saved", HERE / "fig10_agreement.json")


if __name__ == "__main__":
    main()
