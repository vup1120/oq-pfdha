"""Plot hazard curves from the Logic Tree Validation demo.

Overlays, in a single figure:
  1. single_bilinear      - weight-1.0 tree on Petersen2011PrimaryFD_bilinear
  2. single_elliptical    - weight-1.0 tree on Petersen2011PrimaryFD_elliptical
  3. analytical average   - 0.5 * bilinear + 0.5 * elliptical (computed here)
  4. logic-tree blend     - mean column from blend_50_50/out/aggregate_hazard.csv

Curves (3) and (4) must visually coincide; ``verify.py`` already proves they
agree to ~1e-16. The 5-95% fractile band from the logic-tree run is also
drawn so the reader can see that for a 2-branch 50/50 tree it simply tracks
(min, max) of the individual branches under the linear-interpolation
convention used by ``aggregation.weighted_fractiles`` (design memo §9).

Output: hazard_curves_comparison.png next to this script.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def read_col(path: Path, col: str) -> np.ndarray:
    with path.open() as f:
        return np.asarray([float(row[col]) for row in csv.DictReader(f)], dtype=float)


def _positive(x: np.ndarray) -> np.ndarray:
    """Replace non-positive values with NaN so log-scale plots skip them."""
    y = x.astype(float).copy()
    y[y <= 0] = np.nan
    return y


def main() -> None:
    # --- read all four datasets ---
    d0 = read_col(HERE / "single_bilinear" / "out" / "hazard_curves" / "branch_0000.csv", "D0")

    bilinear = read_col(
        HERE / "single_bilinear" / "out" / "hazard_curves" / "branch_0000.csv",
        "annual_rate",
    )
    elliptical = read_col(
        HERE / "single_elliptical" / "out" / "hazard_curves" / "branch_0000.csv",
        "annual_rate",
    )

    blend_agg_path = HERE / "blend_50_50" / "out" / "aggregate_hazard.csv"
    blend_d0 = read_col(blend_agg_path, "D0")
    blend_mean = read_col(blend_agg_path, "mean")
    blend_p05 = read_col(blend_agg_path, "p05")
    blend_p95 = read_col(blend_agg_path, "p95")
    blend_p16 = read_col(blend_agg_path, "p16")
    blend_p84 = read_col(blend_agg_path, "p84")

    if not np.allclose(d0, blend_d0):
        raise SystemExit("D0 grids differ across runs (rerun verify.py first).")

    analytical_avg = 0.5 * bilinear + 0.5 * elliptical

    max_err_vs_lt = np.max(np.abs(blend_mean - analytical_avg))

    # --- plot ---
    fig, ax = plt.subplots(figsize=(8.5, 6.0))

    ax.fill_between(
        d0,
        _positive(blend_p05),
        _positive(blend_p95),
        color="tab:red",
        alpha=0.10,
        label="LT 5-95% fractile band (blend_50_50)",
    )
    ax.fill_between(
        d0,
        _positive(blend_p16),
        _positive(blend_p84),
        color="tab:red",
        alpha=0.15,
        label="LT 16-84% fractile band (blend_50_50)",
    )

    ax.plot(
        d0, _positive(bilinear),
        color="tab:blue", lw=1.8, marker="o", ms=5,
        label="single_bilinear  (weight = 1.0)",
    )
    ax.plot(
        d0, _positive(elliptical),
        color="tab:orange", lw=1.8, marker="s", ms=5,
        label="single_elliptical (weight = 1.0)",
    )
    ax.plot(
        d0, _positive(analytical_avg),
        color="tab:green", lw=2.4, ls="--",
        label=r"analytical  $0.5\cdot\lambda_A + 0.5\cdot\lambda_B$",
    )
    ax.plot(
        d0, _positive(blend_mean),
        color="tab:red", lw=1.2, ls=":", marker="x", ms=7, mew=1.6,
        label="LT framework  blend_50_50.mean",
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Displacement threshold $D_0$ (m)")
    ax.set_ylabel(r"Annual rate of exceedance  $\lambda(D_0)$ (1/yr)")
    ax.set_title(
        "Logic Tree Validation demo — hazard curves\n"
        "Petersen2011PrimaryFD bilinear vs elliptical, 50/50 blend\n"
        r"(strike-slip source, on-trace site; Primary-SR = Pizza2023)"
    )
    ax.grid(True, which="both", ls=":", alpha=0.6)
    ax.legend(loc="lower left", fontsize=9, framealpha=0.95)

    ax.annotate(
        (
            r"$\max\,|\,\mathrm{LT\ mean} - 0.5\cdot\lambda_A - 0.5\cdot\lambda_B\,|$"
            f"\n= {max_err_vs_lt:.2e}  (tol $10^{{-12}}$)"
        ),
        xy=(0.98, 0.98),
        xycoords="axes fraction",
        ha="right",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.6", alpha=0.9),
    )

    out_path = HERE / "hazard_curves_comparison.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Wrote {out_path.relative_to(HERE.parent)}")
    print(f"max |LT.mean - analytical_avg| = {max_err_vs_lt:.3e}")


if __name__ == "__main__":
    main()
