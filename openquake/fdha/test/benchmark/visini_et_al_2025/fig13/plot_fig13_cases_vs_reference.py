# -*- coding: utf-8 -*-
"""
Visini et al. (2025) Fig. 13 — logic-tree hazard curves vs digitised reference CSVs.

Runs each ``job_case*.ini`` via :class:`FdhaLogicTree` (same path as CLI) and plots
mean annual exceedance rates against ``reference_data/visini2025_case*.csv``.

Figures under ``openquake/fdha/test/figures/``:
  - visini2025_fig13_case{1,2,3}_lt_comparison.png
  - visini2025_fig13_case{1,2,3}_lt_relerr.png
  - visini2025_fig13_all_cases_ref_vs_impl.png  (combined log-log panel)
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from openquake.fdha.logic_tree.driver import FdhaLogicTree


def _load_reference_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    arr = np.loadtxt(path, delimiter=",")
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return arr[:, 0].astype(float), arr[:, 1].astype(float)


def _mean_curve_from_ini(ini: Path) -> tuple[np.ndarray, np.ndarray]:
    with tempfile.TemporaryDirectory(prefix="visini_fig13_lt_") as tmp:
        res = FdhaLogicTree.from_ini(ini).run(outdir=Path(tmp))
        d0 = np.asarray(res.d0, dtype=float)
        p = np.asarray(res.mean_rates, dtype=float)
        if p.ndim == 2:
            p = p[0]
    return d0, p


def _interp_loglog(x_need: np.ndarray, x_have: np.ndarray, y_have: np.ndarray) -> np.ndarray:
    """Interpolate y_have(x_have) to x_need in log-log space (strictly positive values)."""
    lx = np.log(x_have)
    ly = np.log(y_have)
    return np.exp(np.interp(np.log(x_need), lx, ly))


def main() -> None:
    base = Path(__file__).resolve().parent
    ref_dir = base / "reference_data"
    figures_dir = base.parent.parent.parent / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    cases = [
        ("case1", base / "job_case1.ini", ref_dir / "visini2025_case1.csv"),
        ("case2", base / "job_case2.ini", ref_dir / "visini2025_case2.csv"),
        ("case3", base / "job_case3.ini", ref_dir / "visini2025_case3.csv"),
    ]

    combined_fig, combined_ax = plt.subplots(figsize=(9, 6))
    colors = {"case1": "#DC143C", "case2": "#4169E1", "case3": "#32CD32"}
    ref_colors = {"case1": "#8B0000", "case2": "#00008B", "case3": "#006400"}

    for case_key, ini_path, ref_path in cases:
        if not ini_path.is_file():
            print(f"[skip] missing INI {ini_path}")
            continue
        d_impl, p_impl = _mean_curve_from_ini(ini_path)
        combined_ax.loglog(
            d_impl,
            p_impl,
            "-",
            color=colors[case_key],
            lw=2,
            label=f"{case_key} computed",
        )

        if ref_path.is_file():
            d_ref, p_ref = _load_reference_csv(ref_path)
            combined_ax.loglog(
                d_ref,
                p_ref,
                "--",
                color=ref_colors[case_key],
                lw=2.5,
                label=f"{case_key} reference",
            )

            p_impl_on_ref = _interp_loglog(d_ref, d_impl, p_impl)
            rel_err = 100.0 * (p_impl_on_ref - p_ref) / (p_ref + 1e-30)
            max_idx = int(np.nanargmax(np.abs(rel_err)))

            plt.figure(figsize=(9, 5.5))
            plt.loglog(d_impl, p_impl, "-", color="#1f77b4", lw=2, label="Computed (LT mean)")
            plt.loglog(d_ref, p_ref, "o", color="#D39200", ms=6, label="Reference (CSV)")
            plt.xlabel("Displacement, d (m)")
            plt.ylabel("Annual exceedance rate (1/yr)")
            plt.title(f"Visini Fig.13 {case_key}")
            plt.grid(True, which="both", ls=":", alpha=0.35)
            plt.legend(loc="best")
            plt.tight_layout()
            cmp_path = figures_dir / f"visini2025_fig13_{case_key}_lt_comparison.png"
            plt.savefig(cmp_path, dpi=150)
            plt.close()

            plt.figure(figsize=(9, 5.5))
            plt.semilogx(d_ref, rel_err, "-o", color="#D39200", ms=5)
            plt.xlabel("Displacement, d (m)")
            plt.ylabel("Relative error (%)")
            plt.title(
                f"{case_key}: max |rel err| = {abs(rel_err[max_idx]):.2f}% "
                f"at d={d_ref[max_idx]:g} m (interp computed to ref d)"
            )
            plt.grid(True, which="both", ls=":", alpha=0.35)
            plt.tight_layout()
            rel_path_fig = figures_dir / f"visini2025_fig13_{case_key}_lt_relerr.png"
            plt.savefig(rel_path_fig, dpi=150)
            plt.close()

            print(f"Saved {cmp_path.name}, {rel_path_fig.name}")
        else:
            print(f"[warn] no reference CSV {ref_path}")

    combined_ax.set_xlabel("Displacement (m)")
    combined_ax.set_ylabel("Annual exceedance rate (1/yr)")
    combined_ax.set_title("Visini et al. (2025) Fig. 13 — all cases")
    combined_ax.grid(True, alpha=0.3, which="both")
    combined_ax.legend(fontsize=9, loc="best")
    combined_fig.tight_layout()
    combo_out = figures_dir / "visini2025_fig13_all_cases_ref_vs_impl.png"
    combined_fig.savefig(combo_out, dpi=150, bbox_inches="tight")
    plt.close(combined_fig)
    print(f"Saved {combo_out}")


if __name__ == "__main__":
    main()
