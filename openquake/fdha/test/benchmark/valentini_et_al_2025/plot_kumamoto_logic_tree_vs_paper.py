# -*- coding: utf-8 -*-
"""
Generate overlay and relative-error plots for Kumamoto Case 2 (Chiou 2025).

Compares the paper reference vector (REF_EXCEED in test_chiou2025_case2.py) against
the logic-tree job.ini mean annual exceedance rates (aggregate hazard).

Usage (from repo root, PYTHONPATH=.):

  python openquake/fdha/test/benchmark/valentini_et_al_2025/plot_kumamoto_logic_tree_vs_paper.py

Figures are written under openquake/fdha/test/figures/:
  - kumamoto_case2_Chiou2025_lt_comparison.png
  - kumamoto_case2_Chiou2025_lt_relerr.png
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from openquake.fdha.logic_tree.driver import FdhaLogicTree

# Same reference as test_chiou2025_case2.py (Valentini et al. 2025 / Kumamoto case 2)
REF_EXCEED = np.array(
    [
        1.31775e-04,
        1.31572e-04,
        1.30400e-04,
        1.28706e-04,
        1.26877e-04,
        1.20965e-04,
        1.12737e-04,
        1.02679e-04,
        9.33119e-05,
        7.71261e-05,
        4.54213e-05,
        2.48653e-05,
        1.32235e-05,
        7.71793e-06,
        4.62204e-07,
        8.05555e-08,
        1.63559e-08,
        4.71018e-09,
    ],
    dtype=float,
)


def main() -> None:
    here = Path(__file__).resolve().parent
    ini = here / "data" / "configuration" / "job_kumamoto_case2_Chiou2025.ini"
    if not ini.is_file():
        raise SystemExit(f"Missing logic-tree job.ini: {ini}")

    figures_dir = here.parent.parent / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="valentini_lt_plot_") as tmp:
        res = FdhaLogicTree.from_ini(ini).run(outdir=Path(tmp))
        d0 = np.asarray(res.d0, dtype=float)
        poes = np.asarray(res.mean_rates, dtype=float)
        if poes.ndim == 2:
            poes = poes[0]

    assert d0.size == REF_EXCEED.size == poes.size

    overlay_png = figures_dir / "kumamoto_case2_Chiou2025_lt_comparison.png"
    plt.figure(figsize=(9, 5.5))
    plt.loglog(d0, REF_EXCEED, color="#D39200", lw=2.5, label="Reference (paper)")
    plt.loglog(d0, poes, "--", color="#1f77b4", lw=2.5, label="Computed (logic tree mean)")
    plt.xlabel("Displacement, d (m)")
    plt.ylabel("Annual exceedance rate (1/yr)")
    plt.grid(True, which="both", ls=":", alpha=0.35)
    plt.legend(loc="lower left")
    plt.tight_layout()
    plt.savefig(overlay_png, dpi=150)
    plt.close()

    rel_err = 100.0 * (poes - REF_EXCEED) / REF_EXCEED
    max_idx = int(np.nanargmax(np.abs(rel_err)))
    relerr_png = figures_dir / "kumamoto_case2_Chiou2025_lt_relerr.png"
    plt.figure(figsize=(9, 5.5))
    plt.semilogx(d0, rel_err, "-o", color="#D39200", ms=4)
    plt.xlabel("Displacement, d (m)")
    plt.ylabel("Relative error (%)")
    plt.title(
        f"Logic tree vs paper - max |rel err| = {abs(rel_err[max_idx]):.3f}% at d={d0[max_idx]:g} m"
    )
    plt.grid(True, which="both", ls=":", alpha=0.35)
    plt.tight_layout()
    plt.savefig(relerr_png, dpi=150)
    plt.close()

    print(f"Saved:\n  {overlay_png}\n  {relerr_png}")


if __name__ == "__main__":
    main()
