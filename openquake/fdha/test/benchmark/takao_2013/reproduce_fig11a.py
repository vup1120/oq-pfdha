#!/usr/bin/env python
"""Reproduce the Takao et al. (2013) principal-faulting example (Fig. 11 case (a)).

Case (a): the site sits at the midpoint of a 22-km active fault
(11 km from the ends). Six hazard curves: recurrence 3,000 / 30,000 yr x
Mw 6.6 / 6.8 / 7.0. Their Eq. (1):

    nu(d) = nu0 * P1p(m) * sum_n (1/nf) sum_models 0.25
            sum_l (1/lp) I(cover) * P3p(D>d | m, x_l/L_l)

with the machinery the paper documents in Section 4(2):

- rupture length from Takemura (1998), log10 L = 0.5 Mw - 2.07
  (gives the paper's 21.38 km at Mw 6.8), discretized to 1 km and capped
  at the 22-km fault; nf = 22 - L + 1 rupture placements;
- principal-segment length ratios Ls/Lm from the four lines of their
  Fig. 3 (equal weight 0.25), anchored at (5.9, 50%), (6.2, 35%),
  (6.5, 20%), (6.8, 5%) and reaching 100% at Mw +0.7, flat outside;
- lp = L - L_l + 1 segment placements per rupture placement, 1-km grid,
  half-open [s, s+L_l) site-coverage convention — this reproduces all
  eight placement counts of the paper's Mw 6.8 worked example;
- P3p from the ``Takao2013PrimaryFD`` building blocks (their Eqs. 6/8 via
  ``get_prob_D_AD`` with the *actual segment length* as srl, and Eq. 10 AD
  scaling via ``get_prob_avg_displacement``), AD normalization as stated
  in the paper.

Published anchors (paper text, Section 4(2), all for the worked example):

    P1p(Mw 6.8) = 0.784        P2p(Mw 6.8) = 0.751
    placement counts (0-21 km rupture): 1/21, 10/11, 3/3, 1/1
    placement counts (1-22 km rupture): 1/21, 11/11, 3/3, 1/1
    nu(0.01 m | Mw 6.6, R = 30,000 yr) = 1.2e-5 /yr

plus the full six digitized curves of Fig. 11(a) in
``reference/fig11a_points.csv`` (gray ink = R 3,000 yr, black = 30,000 yr).

Writes ``Figures/takao2013_fig11a_reproduction.png`` and
``fig11a_agreement.json``.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from openquake.fdha.primary_surf_rup import Takao2013PrimarySR
from openquake.fdha.primary_surf_displ import Takao2013PrimaryFD

from figcompare import load_points, match_stats  # noqa: E402 (script-local)

HERE = Path(__file__).resolve().parent

FAULT_KM = 22
SITE_KM = 11.0
MAGS = (6.6, 6.8, 7.0)
RECURRENCES = (3000.0, 30000.0)

# Fig. 3 segment-length-ratio lines: (Mw_low, ratio_low %); each rises
# linearly to 100% at Mw_low + 0.7 and is flat outside. Weight 0.25 each.
FIG3_LINES = [(5.9, 50.0), (6.2, 35.0), (6.5, 20.0), (6.8, 5.0)]

PAPER_P1P = 0.784
PAPER_P2P = 0.751
PAPER_NU_001_M66_R30000 = 1.2e-5


def rupture_length_km(mag):
    """Takemura (1998) magnitude-length, 1-km discretized, fault-capped."""
    return min(int(10 ** (0.5 * mag - 2.07)), FAULT_KM)


def segment_ratio_pct(mag, m_low, r_low):
    if mag <= m_low:
        return r_low
    return min(100.0, r_low + (100.0 - r_low) * (mag - m_low) / 0.7)


def placements(mag):
    """Yield (weight, segment_length_km, x_over_L or None) covering the site.

    Weight = 1/nf * 0.25 * 1/lp per covered placement; x_over_L is the
    site's normalized distance from the closest segment end.
    """
    L = rupture_length_km(mag)
    nf = FAULT_KM - L + 1
    for n0 in range(nf):  # rupture spans [n0, n0 + L)
        for m_low, r_low in FIG3_LINES:
            Ll = max(1, int(round(segment_ratio_pct(mag, m_low, r_low)
                                  / 100.0 * L)))
            lp = L - Ll + 1
            for s in range(n0, n0 + lp):  # segment spans [s, s + Ll)
                if s <= SITE_KM < s + Ll:
                    x_over_L = min(SITE_KM - s, s + Ll - SITE_KM) / Ll
                    yield (1.0 / nf) * 0.25 * (1.0 / lp), Ll, x_over_L


def paper_placement_counts():
    """Placement counts of the Mw 6.8 worked example, per rupture position."""
    L = rupture_length_km(6.8)
    counts = {}
    for n0 in range(FAULT_KM - L + 1):
        per_model = []
        for m_low, r_low in FIG3_LINES:
            Ll = max(1, int(round(segment_ratio_pct(6.8, m_low, r_low)
                                  / 100.0 * L)))
            lp = L - Ll + 1
            cov = sum(1 for s in range(n0, n0 + lp) if s <= SITE_KM < s + Ll)
            per_model.append((Ll, cov, lp))
        counts[n0] = per_model
    return counts


def p2p(mag):
    return sum(w for w, _, _ in placements(mag))


def p3p_weighted(d, mag, n_sigma=5.0):
    """sum over placements of weight * P3p(D>d | m, x_l/L_l), AD method."""
    model = Takao2013PrimaryFD(n_sigma=n_sigma)
    log_mean = -4.80 + 0.69 * mag
    sigma = 0.36
    lower = 10 ** (log_mean - n_sigma * sigma)
    upper = 10 ** (log_mean + n_sigma * sigma)
    grid = np.logspace(np.log10(lower), np.log10(upper), 1000)
    pdf = np.array([model.get_prob_avg_displacement(ad, mag) for ad in grid])

    # group identical (segment length, x/L) placements
    groups = {}
    for w, Ll, xl in placements(mag):
        groups[(Ll, round(xl, 9))] = groups.get((Ll, round(xl, 9)), 0.0) + w

    d = np.atleast_1d(np.asarray(d, dtype=float))
    total = np.zeros_like(d)
    for (Ll, xl), w in groups.items():
        # P(D > d | AD) marginalized over AD; srl = the actual segment length
        exceed = np.zeros_like(d)
        for ad, p in zip(grid, pdf):
            exceed += p * model.get_prob_D_AD(d / ad, xl, Ll)
        total += w * exceed
    return total


def hazard_curves():
    p1p_model = Takao2013PrimarySR()
    d = np.logspace(-2, 1, 61)
    curves = {}
    for mag in MAGS:
        p1p = float(p1p_model.get_prob(mag))
        p3_weighted = p3p_weighted(d, mag)
        for rec in RECURRENCES:
            curves[(rec, mag)] = (1.0 / rec) * p1p * p3_weighted
    return d, curves


def agreement():
    p1p = float(Takao2013PrimarySR().get_prob(6.8))
    counts = paper_placement_counts()
    # paper's published counts: rupture 0-21: 1/21, 10/11, 3/3, 1/1 for
    # segment lengths 1, 11, 19, 21 km; rupture 1-22: 1/21, 11/11, 3/3, 1/1
    paper_counts = {
        0: [(1, 1, 21), (11, 10, 11), (19, 3, 3), (21, 1, 1)],
        1: [(1, 1, 21), (11, 11, 11), (19, 3, 3), (21, 1, 1)],
    }
    counts_sorted = {n0: sorted(v) for n0, v in counts.items()}
    paper_sorted = {n0: sorted(v) for n0, v in paper_counts.items()}

    d, curves = hazard_curves()
    nu_001_66_30000 = float(np.interp(0.01, d, curves[(30000.0, 6.6)]))

    result = {
        "P1p_M6.8": {"computed": p1p, "paper": PAPER_P1P,
                     "ratio": p1p / PAPER_P1P},
        "P2p_M6.8": {"computed": p2p(6.8), "paper": PAPER_P2P,
                     "ratio": p2p(6.8) / PAPER_P2P},
        "placement_counts_match_paper": counts_sorted == paper_sorted,
        "placement_counts": {str(k): v for k, v in counts_sorted.items()},
        "rupture_length_M6.8_km": rupture_length_km(6.8),
        "nu_001m_M6.6_R30000": {"computed": nu_001_66_30000,
                                "paper": PAPER_NU_001_M66_R30000,
                                "ratio": nu_001_66_30000
                                / PAPER_NU_001_M66_R30000},
    }

    # full-curve comparison against the digitized Fig 11(a) points;
    # gray ink = R 3,000 yr curves, black ink = R 30,000 yr
    for color, rec in [("gray", 3000.0), ("black", 30000.0)]:
        pts = load_points("fig11a", color=color)
        named = {f"Mw{m}_R{int(rec)}": (d, curves[(rec, m)]) for m in MAGS}
        result[f"fig11a_{color}_match"] = match_stats(
            pts, named, y_floor=1.05e-10)
    return result, d, curves


def make_figure(d, curves):
    fig, ax = plt.subplots(figsize=(7, 5.6))
    styles = {6.6: "--", 6.8: "-", 7.0: "-."}
    for (rec, mag), nu in curves.items():
        color = "0.55" if rec == 3000.0 else "k"
        ax.plot(d, nu, styles[mag], color=color, lw=1.4,
                label=f"R={int(rec)}, Mw{mag}")
    for color, mcolor in [("gray", "tab:orange"), ("black", "tab:red")]:
        pts = load_points("fig11a", color=color)
        ax.plot(pts[:, 0], pts[:, 1], ".", ms=3, color=mcolor,
                label=f"paper ({color} ink)", zorder=5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-2, 1e1)
    ax.set_ylim(1e-10, 1e-3)
    ax.set_xlabel("Displacement (m)")
    ax.set_ylabel("Annual rate of exceedance")
    ax.set_title("Takao et al. (2013) Fig. 11 case (a) — computed vs digitized")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    figdir = HERE / "Figures"
    figdir.mkdir(exist_ok=True)
    out = figdir / "takao2013_fig11a_reproduction.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main():
    result, d, curves = agreement()
    out_png = make_figure(d, curves)
    (HERE / "fig11a_agreement.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    print(f"figure: {out_png}")


if __name__ == "__main__":
    import os
    import sys

    os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
    sys.path.insert(0, str(HERE))
    main()
