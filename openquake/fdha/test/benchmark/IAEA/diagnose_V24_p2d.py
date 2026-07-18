#!/usr/bin/env python
"""Decompose the V24 distributed plateau gap: it is the P2d term, and it is
NOT reconcilable by any single scalar - unlike the M11 P_sr gap.

The V24 distributed chain is `FixedPrimarySR(1.0)` (gate) x
`Visini2025SecondarySR` (P2d, the probability of a distributed rupture in
a 100 m cell at distance r) x `Visini2025SecondaryFD` (P3d). The plateau
(d -> 0, where P3d -> 1 and the gate = 1) is therefore

    plateau = rate x P2d(r, mechanism)

so the rate-normalized plateau is the effective P2d. Comparing to the
published V24 curves:

                       Le Teil (r=0.6, reverse)   Norcia (r=7.6, normal)
    our plateau              4.17e-8                    5.02e-8
    paper plateau            2.57e-9                    2.98e-7
    ratio ours/paper           16.2x  (HIGH)             0.17x  (LOW)
    our effective P2d        9.07e-4                    1.24e-4
    paper effective P2d      5.59e-5                    7.39e-4

The two gaps point in OPPOSITE directions: to match Le Teil our curve
would need x0.062, to match Norcia x5.94 - a 96x spread. So no single
gate / P_sr convention can reconcile both (contrast M11, where one P_sr
scalar explained the whole plateau gap). The discrepancy lives in the
distance- and mechanism-dependent P2d *function*: our Visini 2025 P2d
decays faster with distance than the paper's V24 (too high near the
fault, too low far from it).

This confirms the documented conclusion (README "Known deviations"): the
paper's Fig-6 V24 curves were produced with the earlier (2024, under-
review) revision of the Visini et al. model, whereas `Visini2025Secondary*`
implements the final published coefficients and reproduces the
TECDOC-2092 versions of the V24 curves instead (validated directly in
`benchmark/visini_et_al_2025`).

Writes ``Figures/iaea_V24_p2d_decomposition.png``.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent

RATE = {"le_teil": 4.6e-5, "norcia": 4.03281e-4}
CASES = [
    ("le_teil", "Le Teil - distributed (r = 0.6 km, reverse, FW)",
     "fig6b_leteil_distributed.csv"),
    ("norcia", "Norcia - distributed (r = 7.6 km, normal, HW)",
     "fig6c_norcia_distributed.csv"),
]


def read_curve(path, col=None):
    rows = list(csv.DictReader(open(path)))
    x = np.array([float(r.get("disp_m", None) or float(r["disp_cm"]) / 100.0)
                  for r in rows])
    y = np.array([float(r[col if col else "rate"]) for r in rows])
    return x, y


def main():
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    summary = {}
    for ax, (case, title, ref_csv) in zip(axes, CASES):
        do, co = read_curve(HERE / "computed" / f"{case}_distributed_V24.csv")
        dp, cp = read_curve(HERE / "reference" / ref_csv, col="V24")
        rate = RATE[case]
        our_p2d = co[0] / rate
        pap_p2d = cp[0] / rate
        summary[case] = (co[0], cp[0], our_p2d, pap_p2d, cp[0] / co[0])

        ax.loglog(dp * 100, cp, "-", color="tab:green", lw=2.4,
                  label=f"paper V24  (eff. P2d {pap_p2d:.2e})")
        ax.loglog(do * 100, co, "--", color="tab:green", lw=1.9,
                  label=f"ours V24   (eff. P2d {our_p2d:.2e})")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1, 1e3)
        ax.set_ylim(1e-11, 1e-6)
        ax.set_xlabel("Displacement (cm)")
        ax.set_ylabel("AFOE (yr$^{-1}$)")
        ratio = co[0] / cp[0]
        arrow = "HIGH" if ratio > 1 else "LOW"
        ax.set_title(f"{title}\nours {ratio:.1f}x {arrow} at plateau")
        ax.legend(fontsize=8, loc="lower left")
        ax.grid(alpha=0.3, which="both")

    g_lt = summary["le_teil"][4]
    g_no = summary["norcia"][4]
    fig.suptitle(
        "V24 distributed plateau = rate x P2d - the gap is the P2d function, "
        "not a scalar gate\n"
        f"(gate to match Le Teil x{g_lt:.3f}, to match Norcia x{g_no:.2f}: "
        "opposite directions ⇒ model-version difference, not P_sr)",
        fontsize=11)
    fig.tight_layout()
    out = HERE / "Figures" / "iaea_V24_p2d_decomposition.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    for case, (co, cp, op, pp, g) in summary.items():
        print(f"{case:8}: ours {co:.3e}  paper {cp:.3e}  "
              f"ratio {co/cp:.2f}x  gate-to-match x{g:.3f}")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
