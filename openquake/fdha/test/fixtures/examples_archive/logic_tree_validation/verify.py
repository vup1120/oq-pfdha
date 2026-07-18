"""Logic Tree Validation demo - verification script.

Runs the three sibling ``job.ini`` configs (single_bilinear, single_elliptical,
blend_50_50) through ``FdhaLogicTree`` if their outputs are missing, then
checks the aggregation arithmetic demanded by the design memo.

Passes if:
  - D0 grid is consistent across all three runs.
  - mean column of blend_50_50 equals 0.5 * bilinear.mean + 0.5 * elliptical.mean
    elementwise (tolerance 1e-12).
  - 5/16/50 fractiles equal min(A, B) and 84/95 fractiles match the linear
    interpolation used by ``openquake.fdha.logic_tree.aggregation``
    (this is the convention called out in the design memo §9 when
    aggregation.py uses linear interpolation, not step-function).
  - blend_50_50 manifest.json lists both branches with weight 0.5.
  - If ``expected/girs_reference.csv`` exists, single_bilinear.mean agrees
    with it within GIRS_TOL (default 1e-6). Otherwise skip with a [SKIP].

Exit code 0 on success, 1 on any failure.
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
RUNS = {
    "single_bilinear": HERE / "single_bilinear",
    "single_elliptical": HERE / "single_elliptical",
    "blend_50_50": HERE / "blend_50_50",
}
OUT_SUBDIR = "out"  # FdhaLogicTree default out dir relative to the INI


class CheckFailed(SystemExit):
    def __init__(self, msg: str) -> None:
        super().__init__(f"[FAIL] {msg}")


def _banner(title: str) -> None:
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def _ensure_runs() -> None:
    """Run any INI whose out/aggregate_hazard.csv is missing."""
    from openquake.fdha.logic_tree.driver import FdhaLogicTree

    for name, run_dir in RUNS.items():
        out = run_dir / OUT_SUBDIR / "aggregate_hazard.csv"
        if out.is_file():
            print(f"  [skip run] {name}: found existing {out.relative_to(HERE)}")
            continue
        print(f"  [run]      {name}: computing ...")
        FdhaLogicTree.from_ini(str(run_dir / "job.ini")).run(
            outdir=run_dir / OUT_SUBDIR
        )


def _read_aggregate(csv_path: Path) -> dict[str, np.ndarray]:
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        cols = reader.fieldnames or []
        rows = list(reader)
    data: dict[str, list[float]] = {c: [] for c in cols}
    for row in rows:
        for c in cols:
            data[c].append(float(row[c]))
    return {c: np.asarray(v, dtype=float) for c, v in data.items()}


def _load_manifest(run_dir: Path) -> dict:
    return json.loads((run_dir / OUT_SUBDIR / "manifest.json").read_text())


def main() -> int:
    _banner("Logic Tree Validation demo - verify.py")

    _ensure_runs()

    A = _read_aggregate(RUNS["single_bilinear"] / OUT_SUBDIR / "aggregate_hazard.csv")
    B = _read_aggregate(RUNS["single_elliptical"] / OUT_SUBDIR / "aggregate_hazard.csv")
    X = _read_aggregate(RUNS["blend_50_50"] / OUT_SUBDIR / "aggregate_hazard.csv")

    required_cols = {"D0", "mean", "p05", "p16", "p50", "p84", "p95"}
    for label, data in (("single_bilinear", A), ("single_elliptical", B), ("blend_50_50", X)):
        missing = required_cols - set(data.keys())
        if missing:
            raise CheckFailed(f"{label}/aggregate_hazard.csv missing columns: {missing}")

    # --- D0 grid consistency ---
    _banner("Check 1 | D0 grid consistency across runs")
    if not (np.allclose(A["D0"], X["D0"]) and np.allclose(B["D0"], X["D0"])):
        raise CheckFailed("D0 grids differ across runs")
    print(f"  [PASS] D0 grid identical. n_D0 = {len(A['D0'])}")

    # --- Check 2: aggregation arithmetic ---
    _banner("Check 2 | Aggregation arithmetic: mean(X) = 0.5 * mean(A) + 0.5 * mean(B)")
    expected_mean = 0.5 * A["mean"] + 0.5 * B["mean"]
    diff = np.abs(X["mean"] - expected_mean)
    if diff.max() >= 1e-12:
        raise CheckFailed(
            f"Aggregation arithmetic failed. max|X.mean - 0.5*A - 0.5*B| "
            f"= {diff.max():.3e}"
        )
    print(f"  [PASS] max|diff| = {diff.max():.3e}  (tol 1e-12)")

    # --- Check 3: fractile behaviour for the 2-branch 50/50 tree ---
    #
    # aggregation.py uses linear interpolation of the weighted empirical CDF
    # (see openquake/fdha/logic_tree/aggregation.py::weighted_fractiles).
    # For two equal-weight branches the CDF breakpoints are at cum = [0.5, 1.0]
    # so the expected fractile values are:
    #
    #   q <= 0.5   -> min(A, B)
    #   q >  0.5   -> min + (q - 0.5)/(1.0 - 0.5) * (max - min)
    #
    # which gives:
    #   p05 = p16 = p50 = min
    #   p84 = min + 0.68 * (max - min)
    #   p95 = min + 0.90 * (max - min)
    #
    # (The design memo §9 explicitly permits updating the CHECK to match the
    # implementation when linear interpolation is used.)
    _banner("Check 3 | Fractile reduction for 2-branch 50/50 (linear interpolation)")
    lo = np.minimum(A["mean"], B["mean"])
    hi = np.maximum(A["mean"], B["mean"])
    expected = {
        "p05": lo,
        "p16": lo,
        "p50": lo,
        "p84": lo + 0.68 * (hi - lo),
        "p95": lo + 0.90 * (hi - lo),
    }
    for q, expected_vals in expected.items():
        err = np.abs(X[q] - expected_vals)
        if err.max() >= 1e-12:
            raise CheckFailed(
                f"{q} fractile mismatch. max|err| = {err.max():.3e}"
            )
        print(f"  [PASS] {q}: max|err| = {err.max():.3e}")

    # --- Check 4: blend_50_50 manifest lists both branches with weight 0.5 ---
    _banner("Check 4 | blend_50_50 manifest: 2 branches, each weight 0.5")
    manifest = _load_manifest(RUNS["blend_50_50"])
    branches = manifest.get("branches", [])
    if len(branches) != 2:
        raise CheckFailed(
            f"blend_50_50 manifest has {len(branches)} branches (expected 2)"
        )
    weights = sorted(b["weight"] for b in branches)
    if not np.allclose(weights, [0.5, 0.5], atol=1e-12):
        raise CheckFailed(f"blend_50_50 branch weights = {weights} (expected [0.5, 0.5])")
    # Both branches now use the single Petersen2011PrimaryFD class,
    # distinguished by the `version` parameter (bilinear vs elliptical)
    # rather than by separate subclasses.
    pfd_classes = {b["models"]["primary_surf_displ"] for b in branches}
    if pfd_classes != {"Petersen2011PrimaryFD"}:
        raise CheckFailed(
            f"blend_50_50 primary_surf_displ classes = {pfd_classes} "
            "(expected both branches to be Petersen2011PrimaryFD)"
        )
    print(f"  [PASS] 2 branches, weights = {weights}, PFDs = {sorted(pfd_classes)}")

    # --- Check 5 (optional): GIRS reference regression ---
    _banner("Check 5 | GIRS reference regression (optional)")
    ref_path = HERE / "expected" / "girs_reference.csv"
    if not ref_path.is_file():
        print(f"  [SKIP] no reference file at {ref_path.relative_to(HERE)}")
    else:
        ref = _read_aggregate(ref_path)
        if "annual_rate" not in ref:
            raise CheckFailed(
                f"{ref_path}: missing column 'annual_rate' (expected D0,annual_rate)"
            )
        if not np.allclose(ref["D0"], A["D0"]):
            raise CheckFailed("GIRS reference D0 grid does not match single_bilinear")
        tol = float(os.environ.get("GIRS_TOL", "1e-6"))
        err = np.abs(A["mean"] - ref["annual_rate"])
        if err.max() >= tol:
            raise CheckFailed(
                f"GIRS regression failed. max|diff| = {err.max():.3e} > TOL={tol:.1e}"
            )
        print(f"  [PASS] max|diff| = {err.max():.3e} (TOL={tol:.1e})")

    _banner("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
