# -*- coding: utf-8 -*-
"""
Overlay + relative-error plots for Norcia Youngs2003 AD 85 benchmark.

Uses the canonical logic-tree INI (``job_norcia_youngs2003_AD_85.ini``) and
compares the aggregated mean rates to the tabulated paper reference used in
``plot_results.py`` / ``compare_results.py``.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from openquake.fdha.logic_tree.driver import FdhaLogicTree

# Reference from plot_results.py (Youngs et al. 2003, AD 85 %ile context)
REF_RATES = np.array(
    [
        3.29e-04,
        3.25e-04,
        3.16e-04,
        3.08e-04,
        3.01e-04,
        2.85e-04,
        2.67e-04,
        2.49e-04,
        2.34e-04,
        2.08e-04,
        1.56e-04,
        1.14e-04,
        8.16e-05,
        6.15e-05,
        1.31e-05,
        4.72e-06,
        1.79e-06,
        8.09e-07,
    ],
    dtype=float,
)


def main() -> None:
    here = Path(__file__).resolve().parent
    ini = here / "job_norcia_youngs2003_AD_85.ini"
    if not ini.is_file():
        raise SystemExit(f"Missing job.ini: {ini}")

    figures_dir = here.parent.parent / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="norcia_lt_plot_") as tmp:
        res = FdhaLogicTree.from_ini(ini).run(outdir=Path(tmp))
        d0 = np.asarray(res.d0, dtype=float)
        poes = np.asarray(res.mean_rates, dtype=float)
        if poes.ndim == 2:
            poes = poes[0]

    assert d0.size == REF_RATES.size == poes.size

    overlay_png = figures_dir / "norcia_youngs2003_AD_85_lt_comparison.png"
    plt.figure(figsize=(9, 5.5))
    plt.loglog(d0, REF_RATES, color="#D39200", lw=2.5, label="Reference (paper)")
    plt.loglog(d0, poes, "--", color="#1f77b4", lw=2.5, label="Computed (logic tree mean)")
    plt.xlabel("Displacement, d (m)")
    plt.ylabel("Annual exceedance rate (1/yr)")
    plt.grid(True, which="both", ls=":", alpha=0.35)
    plt.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(overlay_png, dpi=150)
    plt.close()

    rel_err = 100.0 * (poes - REF_RATES) / REF_RATES
    max_idx = int(np.nanargmax(np.abs(rel_err)))
    relerr_png = figures_dir / "norcia_youngs2003_AD_85_lt_relerr.png"
    plt.figure(figsize=(9, 5.5))
    plt.semilogx(d0, rel_err, "-o", color="#D39200", ms=4)
    plt.xlabel("Displacement, d (m)")
    plt.ylabel("Relative error (%)")
    plt.title(
        f"Norcia vs paper - max |rel err| = {abs(rel_err[max_idx]):.3f}% at d={d0[max_idx]:g} m"
    )
    plt.grid(True, which="both", ls=":", alpha=0.35)
    plt.tight_layout()
    plt.savefig(relerr_png, dpi=150)
    plt.close()

    print(f"Saved:\n  {overlay_png}\n  {relerr_png}")


if __name__ == "__main__":
    main()
