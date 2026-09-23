#!/usr/bin/env python
"""Reproduce the Takao et al. (2013) distributed-faulting example (their Fig. 11 case (b)).

Takao et al. (2013, J. Japan Assoc. Earthq. Eng. 13(1), 17-36) illustrate
their PFDHA formulation on two model cases. Case (b) is the distributed
(secondary) faulting case: the site is NOT on an active fault, and four
surrounding faults contribute:

    fault 1: r =  5 km, recurrence  3,000 yr
    fault 2: r =  5 km, recurrence 30,000 yr
    fault 3: r = 10 km, recurrence  3,000 yr
    fault 4: r = 10 km, recurrence 30,000 yr

all with a single characteristic magnitude Mw 6.8. Their Eq. (13) evaluates

    nu_i(d) = (1/R_i) * P1p(m) * P2d(m, r_i) * P3d(D>d | m, r_i)

(the 1/nf rupture-placement sum collapses because the site-to-rupture
distance equals the site-to-fault distance for every placement). The three
probability terms map onto the library models

    P1p -> Takao2013PrimarySR    (their Eq. 4;  P1p(6.8) = 0.784 in the paper)
    P2d -> Takao2013SecondarySR  (their Eq. 14)
    P3d -> Takao2013SecondaryFD  (their Eqs. 15-17, AD normalization as in
                                  their case study)

Published anchors used for the comparison (paper text, Section 4(2)):

    sum over the four faults at d = 0.01 m  ->  7.0e-7 /yr  ("once in
    1.43 million years")
    P1p(Mw 6.8) = 0.784

plus two structural features of their Fig. 11 case (b): the fault-2 and
fault-3 curves cross ("hazard curves intersect and reverse"), and same-
distance curves differ exactly by the 10x recurrence ratio.

Writes ``Figures/takao2013_fig11b_reproduction.png`` and
``fig11b_agreement.json`` next to this script.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from openquake.pfd.primary_surf_rup import Takao2013PrimarySR
from openquake.pfd.secondary_surf_rup import Takao2013SecondarySR
from openquake.pfd.secondary_surf_displ import Takao2013SecondaryFD

HERE = Path(__file__).resolve().parent

MAG = 6.8
FAULTS = {  # name: (distance km, recurrence yr)
    "fault 1": (5.0, 3000.0),
    "fault 2": (5.0, 30000.0),
    "fault 3": (10.0, 3000.0),
    "fault 4": (10.0, 30000.0),
}

PAPER_P1P = 0.784
PAPER_SUM_001M = 7.0e-7  # /yr at d = 0.01 m


def hazard_curves(norm_disp_type="AD"):
    """Fig. 11 case (b) hazard curves; returns (d, {name: nu}, sum).

    n_sigma = 5 approximates the paper's untruncated Eq. 12 integral (the
    curves are converged: n_sigma 5 and 10 are indistinguishable); the
    default 3 visibly starves the tails below AFOE ~ 1e-9.
    """
    p1p_model = Takao2013PrimarySR()
    p2d_model = Takao2013SecondarySR()
    p3d_model = Takao2013SecondaryFD(n_sigma=5.0)

    d = np.logspace(-2, 1, 61)  # 0.01 - 10 m, as plotted in the paper
    p1p = float(p1p_model.get_prob(MAG))

    curves = {}
    for name, (r_km, recurrence) in FAULTS.items():
        p2d = float(np.atleast_1d(p2d_model.get_prob(r_km, MAG))[0])
        p3d = p3d_model.get_prob(d, MAG, r_km, norm_disp_type)[0]
        curves[name] = (1.0 / recurrence) * p1p * p2d * p3d

    total = np.sum(list(curves.values()), axis=0)
    return d, p1p, curves, total


def agreement():
    import sys
    sys.path.insert(0, str(HERE))
    from figcompare import load_points, match_stats

    d, p1p, curves, total = hazard_curves()
    sum_001 = float(np.interp(0.01, d, total))
    diff23 = curves["fault 2"] - curves["fault 3"]

    # full-curve comparison against the digitized Fig 11(b) points (all
    # black ink; the five curves are resolved by nearest-curve assignment)
    named = {k.replace(" ", ""): (d, v) for k, v in curves.items()}
    named["sum"] = (d, total)
    fig11b_match = match_stats(load_points("fig11b"), named, y_floor=1.05e-10)

    return {
        "P1p_M6.8": {"computed": p1p, "paper": PAPER_P1P,
                     "ratio": p1p / PAPER_P1P},
        "sum_at_0.01m": {"computed": sum_001, "paper": PAPER_SUM_001M,
                         "ratio": sum_001 / PAPER_SUM_001M},
        "fault2_fault3_curves_cross": bool(np.nanmin(diff23) < 0 < np.nanmax(diff23)),
        "fault1_over_fault2": float(np.median(curves["fault 1"] / curves["fault 2"])),
        "fig11b_match": fig11b_match,
    }


def make_figure():
    d, _, curves, total = hazard_curves()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))

    # ---- left: Eqs 15/16 90% non-exceedance levels vs the implemented gamma
    p3d = Takao2013SecondaryFD()
    r = np.linspace(0.0, 20.0, 200)
    from scipy.stats import gamma as gamma_dist
    for ndt, c, color in [("MD", 0.55, "tab:blue"), ("AD", 1.9, "tab:red")]:
        level = c * np.exp(-0.17 * r)
        b = level / p3d._p90_factor
        model_p90 = gamma_dist.ppf(0.90, p3d._GAMMA_SHAPE, scale=b)
        ax1.plot(r, level, color=color, lw=2,
                 label=f"Eq. {15 if ndt == 'MD' else 16}: DD/P{ndt}")
        ax1.plot(r, model_p90, "--", color="k", lw=1)
    ax1.set_yscale("log")
    ax1.set_xlabel("Distance r (km)")
    ax1.set_ylabel("90% non-exceedance DD/PMD, DD/PAD")
    ax1.set_title("Eqs. 15-16 vs implemented gamma p90 (dashed)")
    ax1.legend()
    ax1.grid(alpha=0.3, which="both")

    # ---- right: Fig 11 case (b)
    styles = {"fault 1": ":", "fault 2": "--", "fault 3": "-.", "fault 4": (0, (3, 1, 1, 1))}
    for name, nu in curves.items():
        ax2.plot(d, nu, linestyle=styles[name], color="k", label=name)
    ax2.plot(d, total, "-", color="k", lw=2, label="sum.")
    import sys
    sys.path.insert(0, str(HERE))
    from figcompare import load_points
    pts = load_points("fig11b")
    ax2.plot(pts[:, 0], pts[:, 1], ".", ms=3, color="tab:red", zorder=5,
             label="paper (digitized)")
    ax2.plot([0.01], [PAPER_SUM_001M], "o", color="tab:blue", zorder=6,
             label="paper text: 7.0e-7 at 0.01 m")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(1e-2, 1e1)
    ax2.set_ylim(1e-10, 1e-3)
    ax2.set_xlabel("Displacement (m)")
    ax2.set_ylabel("Annual rate of exceedance")
    ax2.set_title("Takao et al. (2013) Fig. 11 case (b)")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3, which="both")

    fig.tight_layout()
    figdir = HERE / "Figures"
    figdir.mkdir(exist_ok=True)
    out = figdir / "takao2013_fig11b_reproduction.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main():
    result = agreement()
    out_png = make_figure()
    (HERE / "fig11b_agreement.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    print(f"figure: {out_png}")


if __name__ == "__main__":
    import os

    os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
    main()
