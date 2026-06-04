"""Logic Tree Validation (reverse) — verification script.

Same invariants as the strike-slip sibling demo, rebound onto the reverse-fault
configuration (Youngs2003PrimaryFD vs Moss2024PrimaryFD, 50/50 blend).

Passes if:
  - D0 grid is consistent across the three runs.
  - mean column of blend_50_50 equals 0.5 * youngs.mean + 0.5 * moss.mean
    elementwise (tolerance 1e-12).
  - 5/16/50 fractiles equal min(A, B) and 84/95 fractiles match the linear
    interpolation used by openquake.fdha.logic_tree.aggregation.
  - blend_50_50 manifest lists both branches with weight 0.5.

Exit code 0 on success, 1 on any failure.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
RUNS = {
    "single_youngs2003": HERE / "single_youngs2003",
    "single_moss2024":   HERE / "single_moss2024",
    "blend_50_50":       HERE / "blend_50_50",
}
OUT_SUBDIR = "out"
EXPECTED_PFDS = {"Youngs2003PrimaryFD", "Moss2024PrimaryFD"}


class CheckFailed(SystemExit):
    def __init__(self, msg: str) -> None:
        super().__init__(f"[FAIL] {msg}")


def _banner(title: str) -> None:
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def _ensure_runs() -> None:
    from openquake.fdha.logic_tree.driver import FdhaLogicTree

    for name, run_dir in RUNS.items():
        out = run_dir / OUT_SUBDIR / "aggregate_hazard.csv"
        if out.is_file():
            print(f"  [skip run] {name}: found existing {out.relative_to(HERE)}")
            continue
        print(f"  [run]      {name}: computing ...")
        FdhaLogicTree.from_ini(str(run_dir / "fdha.ini")).run(
            outdir=run_dir / OUT_SUBDIR
        )


def _read_aggregate(csv_path: Path) -> dict[str, np.ndarray]:
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        cols = reader.fieldnames or []
        rows = list(reader)
    data = {c: [float(r[c]) for r in rows] for c in cols}
    return {c: np.asarray(v, dtype=float) for c, v in data.items()}


def _load_manifest(run_dir: Path) -> dict:
    return json.loads((run_dir / OUT_SUBDIR / "manifest.json").read_text())


def main() -> int:
    _banner("Logic Tree Validation (reverse) demo — verify.py")

    _ensure_runs()

    A = _read_aggregate(RUNS["single_youngs2003"] / OUT_SUBDIR / "aggregate_hazard.csv")
    B = _read_aggregate(RUNS["single_moss2024"]   / OUT_SUBDIR / "aggregate_hazard.csv")
    X = _read_aggregate(RUNS["blend_50_50"]       / OUT_SUBDIR / "aggregate_hazard.csv")

    required_cols = {"D0", "mean", "p05", "p16", "p50", "p84", "p95"}
    for label, data in (("single_youngs2003", A),
                        ("single_moss2024", B),
                        ("blend_50_50", X)):
        missing = required_cols - set(data.keys())
        if missing:
            raise CheckFailed(f"{label}/aggregate_hazard.csv missing columns: {missing}")

    _banner("Check 1 | D0 grid consistency across runs")
    if not (np.allclose(A["D0"], X["D0"]) and np.allclose(B["D0"], X["D0"])):
        raise CheckFailed("D0 grids differ across runs")
    print(f"  [PASS] D0 grid identical. n_D0 = {len(A['D0'])}")

    _banner("Check 2 | Aggregation arithmetic: mean(X) = 0.5 * mean(A) + 0.5 * mean(B)")
    expected_mean = 0.5 * A["mean"] + 0.5 * B["mean"]
    diff = np.abs(X["mean"] - expected_mean)
    if diff.max() >= 1e-12:
        raise CheckFailed(
            f"Aggregation arithmetic failed. max|X.mean - 0.5*A - 0.5*B| "
            f"= {diff.max():.3e}"
        )
    print(f"  [PASS] max|diff| = {diff.max():.3e}  (tol 1e-12)")

    # Fractiles under linear-interp convention for two equal-weight branches.
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
            raise CheckFailed(f"{q} fractile mismatch. max|err| = {err.max():.3e}")
        print(f"  [PASS] {q}: max|err| = {err.max():.3e}")

    _banner("Check 4 | blend_50_50 manifest: 2 branches, each weight 0.5")
    manifest = _load_manifest(RUNS["blend_50_50"])
    branches = manifest.get("branches", [])
    if len(branches) != 2:
        raise CheckFailed(f"blend_50_50 manifest has {len(branches)} branches (expected 2)")
    weights = sorted(b["weight"] for b in branches)
    if not np.allclose(weights, [0.5, 0.5], atol=1e-12):
        raise CheckFailed(f"blend_50_50 branch weights = {weights} (expected [0.5, 0.5])")
    pfd_classes = {b["models"]["primary_surf_displ"] for b in branches}
    if pfd_classes != EXPECTED_PFDS:
        raise CheckFailed(
            f"blend_50_50 primary_surf_displ classes = {pfd_classes} "
            f"(expected {EXPECTED_PFDS})"
        )
    styles = {b["style"] for b in branches}
    if styles != {"reverse"}:
        raise CheckFailed(f"blend_50_50 styles = {styles} (expected only 'reverse')")
    print(f"  [PASS] 2 branches, weights = {weights}, PFDs = {sorted(pfd_classes)}, "
          f"style = reverse")

    _banner("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
