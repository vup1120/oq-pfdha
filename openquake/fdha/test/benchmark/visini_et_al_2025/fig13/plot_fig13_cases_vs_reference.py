# -*- coding: utf-8 -*-
"""
Visini et al. (2025) Fig. 13 — logic-tree hazard curves vs digitised reference CSVs.

Runs each ``job_case*.ini`` via :class:`FdhaLogicTree` (same path as CLI) and plots
mean annual exceedance rates against ``reference_data/visini2025_case*.csv``.

Figures under ``benchmark/visini_et_al_2025/Figures/``:
  - visini2025_fig13_case{1,2,3}_lt_comparison.png
  - visini2025_fig13_case{1,2,3}_lt_relerr.png
  - visini2025_fig13_all_cases_ref_vs_impl.png  (combined log-log panel,
    the headline figure reproducing the paper's Fig. 13)

Also writes ``fig13_agreement.json`` next to this script with per-case
computed/reference ratio statistics over the same comparison window as
``test_fig13_reproduction.py`` (inside the digitized range and below the
+3-sigma truncation cliff at d <= 2 m).

Note on the 1-3 m gap in the job's displacement grid: ``job_case*.ini``'s
``displacement_measure_levels`` jumps straight from 1.0 m to 3.0 m. The FD
model (``Visini2025SecondaryFD``) is a *truncated* log-normal (``n_sigma =
3``, matching the FDHLab MATLAB reference convention) — beyond its 3-sigma
bound the exceedance probability is exactly zero by construction, and for
these scenarios that bound falls around d ~ 2.1-2.6 m depending on the case.
Plotting only the coarse grid draws a straight (and on a log-y axis,
near-vertical) line from the last nonzero point at 1.0 m to the exact zero
at 3.0 m, which looks like a modelling error but is a sampling artifact: the
model matches the digitized reference right up to where the reference itself
stops (see REFERENCE.md / README.md). ``_mean_curve_from_ini`` below inserts
extra points across that gap for the **plot only** — it does not touch the
golden ``job_case*.ini`` files or the pytest comparison grid.
"""

from __future__ import annotations

import json
import re
import shutil
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from openquake.fdha.logic_tree.driver import FdhaLogicTree

# Extra displacement levels (m) inserted only for plotting, to resolve the
# gaps in job_case*.ini's displacement_measure_levels (1.0->3.0 and
# 3.0->10.0 m) across each case's Visini2025SecondaryFD truncation cliff.
# The cliff location scales with each case's TPFm/combination-dependent
# median displacement -- case1 (A+B+C) truncates near 8-9 m, case2/case3
# (A+B / A) near 2.5 m -- so this covers both regions log-spaced.
_PLOT_FILL_LEVELS = sorted(set(
    [round(x, 3) for x in np.geomspace(1.05, 2.9, 20)]
    + [round(x, 3) for x in np.geomspace(3.1, 9.8, 20)]
))


def _load_reference_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    arr = np.loadtxt(path, delimiter=",")
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return arr[:, 0].astype(float), arr[:, 1].astype(float)


def _densify_ini_text(text: str) -> str:
    """Merge ``_PLOT_FILL_LEVELS`` into the ini's ``displacement_measure_levels``."""
    match = re.search(r'displacement_measure_levels\s*=\s*(\{.*\})', text)
    if not match:
        return text
    dml = json.loads(match.group(1))
    levels = sorted(set(dml["FD"]) | set(_PLOT_FILL_LEVELS))
    dml["FD"] = levels
    return text[:match.start(1)] + json.dumps(dml) + text[match.end(1):]


def _mean_curve_from_ini(ini: Path, *, densify: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """Run ``ini`` and return its (displacement, mean-rate) curve.

    ``densify=True`` merges ``_PLOT_FILL_LEVELS`` into the displacement grid
    first (see module docstring) — used only for the plotted line, never for
    the reported agreement statistics, so those stay tied to the exact grid
    ``test_fig13_reproduction.py`` asserts against.
    """
    with tempfile.TemporaryDirectory(prefix="visini_fig13_lt_") as tmp:
        tmp_dir = Path(tmp)
        if densify:
            # Copy the whole case folder so relative includes (fdha/source-
            # model logic-tree files, rank1p5 traces, ...) resolve unchanged.
            for f in ini.parent.iterdir():
                if f.is_file():
                    shutil.copy(f, tmp_dir / f.name)
            patched_ini = tmp_dir / ini.name
            patched_ini.write_text(_densify_ini_text(ini.read_text()))
            ini = patched_ini

        res = FdhaLogicTree.from_ini(ini).run(outdir=tmp_dir / "out")
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
    figures_dir = base.parent / "Figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    agreement: dict[str, dict] = {}

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
        # Coarse grid: exactly what test_fig13_reproduction.py runs and
        # asserts against — used only for the reported agreement stats.
        d_coarse, p_coarse = _mean_curve_from_ini(ini_path, densify=False)
        # Densified grid: fills the coarse grid's 1-3 m gap so the plotted
        # line shows the true roll-off instead of a straight line jumping to
        # the model's exact truncation zero (see module docstring). Used for
        # the plot and the relative-error panel only, never for `agreement`.
        d_impl, p_impl = _mean_curve_from_ini(ini_path, densify=True)
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

            # Ratio statistics over the same window as the pytest benchmark
            # (calculation levels inside the digitized range, d <= 2 m, below
            # the +3-sigma truncation cliff of Combinations A/B), computed
            # from the coarse (un-densified) grid so these numbers match what
            # test_fig13_reproduction.py actually asserts.
            mask = (d_coarse >= d_ref.min()) & (d_coarse <= min(2.0, d_ref.max())) \
                & (p_coarse > 0)
            ref_on_coarse = _interp_loglog(d_coarse[mask], d_ref,
                                           np.maximum(p_ref, 1e-12))
            ratio = p_coarse[mask] / ref_on_coarse
            agreement[case_key] = {
                "n_points": int(mask.sum()),
                "ratio_min": round(float(ratio.min()), 4),
                "ratio_median": round(float(np.median(ratio)), 4),
                "ratio_max": round(float(ratio.max()), 4),
            }
            print(f"{case_key}: computed/reference ratio over d <= 2 m — "
                  f"min {ratio.min():.3f}, median {np.median(ratio):.3f}, "
                  f"max {ratio.max():.3f} ({mask.sum()} points)")

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
    # Scenario rate = 1/yr and P(SR_primary) = 1, so the mean rate equals the
    # conditional probability of exceedance published in Fig. 13.
    combined_ax.set_ylabel("Conditional probability of exceedance")
    combined_ax.set_title(
        "Visini et al. (2025) Fig. 13 — decision-tree cases 1-3:\n"
        "published curves (dashed) vs oq-pfdha (solid)")
    combined_ax.grid(True, alpha=0.3, which="both")
    combined_ax.legend(fontsize=9, loc="best")
    combined_fig.tight_layout()
    combo_out = figures_dir / "visini2025_fig13_all_cases_ref_vs_impl.png"
    combined_fig.savefig(combo_out, dpi=150, bbox_inches="tight")
    plt.close(combined_fig)
    print(f"Saved {combo_out}")

    agreement_path = base / "fig13_agreement.json"
    agreement_path.write_text(json.dumps(agreement, indent=2) + "\n")
    print(f"Saved {agreement_path}")


if __name__ == "__main__":
    main()
