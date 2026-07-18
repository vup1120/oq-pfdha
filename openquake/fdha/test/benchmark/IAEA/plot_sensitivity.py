#!/usr/bin/env python
"""Overlay the IAEA sensitivity-case demonstration curves per study.

Reads ``computed/<case>_<job>.csv`` (written by ``run_sensitivity.py``)
and draws one panel per case study, writing
``Figures/iaea_sensitivity_cases.png``.
"""

from __future__ import annotations

import csv
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
from sensitivity import SENSITIVITY_JOBS  # noqa: E402


def load(tag):
    rows = list(csv.DictReader((HERE / "computed" / f"{tag}.csv").open()))
    return (np.array([float(r["disp_m"]) for r in rows]),
            np.array([float(r["rate"]) for r in rows]))


def main():
    cases = ["kumamoto", "le_teil", "norcia"]
    titles = {"kumamoto": "Kumamoto (strike-slip)",
              "le_teil": "Le Teil (reverse)",
              "norcia": "Norcia (normal)"}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    for ax, case in zip(axes, cases):
        for sj in SENSITIVITY_JOBS:
            if sj.case != case:
                continue
            d, r = load(f"{sj.case}_{sj.job}")
            m = r > 0
            ax.plot(d[m], r[m], lw=1.6,
                    label=sj.job.replace("distributed_", "").replace(
                        "principal_", ""))
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1e-3, 1e1)
        ax.set_ylim(1e-12, 1e-3)
        ax.set_xlabel("Displacement (m)")
        ax.set_ylabel("Annual rate of exceedance")
        ax.set_title(f"{titles[case]} - sensitivity cases")
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3, which="both")
    fig.suptitle("IAEA PFDHA exercise - sensitivity-case demonstration curves "
                 "(no published reference; snapshot-validated)", y=1.02)
    fig.tight_layout()
    figdir = HERE / "Figures"
    figdir.mkdir(exist_ok=True)
    out = figdir / "iaea_sensitivity_cases.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
