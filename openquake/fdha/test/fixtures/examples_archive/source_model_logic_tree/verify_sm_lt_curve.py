#!/usr/bin/env python3
"""Numerical + visual check: SMLT weighted mean vs analytic rate aggregation.

Reads ``out_curve/manifest.json`` (hazard_curve + source-model logic tree) and
each per-realisation ``hazard_curves/branch_*.csv``.  Rebuilds the mean as

    analytical_mean(D0) = Σ_i w_i · λ_i(D0)

where ``w_i`` is ``combined_branch_weight`` from the manifest, and compares to
the ``mean`` column of ``aggregate_hazard.csv`` produced by
:class:`openquake.fdha.logic_tree.driver.FdhaLogicTree`.

Also writes ``smlt_validation_hazard_curves.png`` (matplotlib) in this example
directory - same style as ``examples/logic_tree_validation/plot_hazard_curves.py``.

Usage::

    # After generating outputs (recommended: separate dirs)
    python verify_sm_lt_curve.py

    # Or (re)compute curve outputs first
    python verify_sm_lt_curve.py --run
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TOL = 1e-12


def _positive(x: np.ndarray) -> np.ndarray:
    y = x.astype(float).copy()
    y[y <= 0] = np.nan
    return y


def _read_single_site_curve(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Return (d0, annual_rate) from a one-site branch CSV."""
    d0_list: list[float] = []
    r_list: list[float] = []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            d0_list.append(float(row["D0"]))
            r_list.append(float(row["annual_rate"]))
    return np.asarray(d0_list, dtype=float), np.asarray(r_list, dtype=float)


def _read_aggregate(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="") as f:
        r = csv.DictReader(f)
        rows = list(r)
        if not rows:
            raise SystemExit(f"Empty aggregate CSV: {path}")
        keys = rows[0].keys()
        out: dict[str, list[float]] = {k: [] for k in keys}
        for row in rows:
            for k in keys:
                out[k].append(float(row[k]))
        return {k: np.asarray(v, dtype=float) for k, v in out.items()}


def _ensure_out_curve(out_curve: Path, run: bool) -> None:
    if (out_curve / "manifest.json").is_file() and (out_curve / "aggregate_hazard.csv").is_file():
        return
    if not run:
        raise SystemExit(
            f"Missing {out_curve}. Run with --run or execute:\n"
            f"  FdhaLogicTree.from_ini(...).run(outdir={out_curve!r})"
        )
    from openquake.fdha.logic_tree.driver import FdhaLogicTree

    FdhaLogicTree.from_ini(HERE / "fdha_curve.ini").run(outdir=out_curve)
    print(f"Wrote hazard-curve outputs -> {out_curve}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--out-curve",
        type=Path,
        default=HERE / "out_curve",
        help="Directory with manifest + aggregate_hazard.csv (default: ./out_curve)",
    )
    ap.add_argument(
        "--plot",
        type=Path,
        default=HERE / "smlt_validation_hazard_curves.png",
        help="Output PNG path",
    )
    ap.add_argument(
        "--run",
        action="store_true",
        help="If outputs are missing, run fdha_curve.ini into --out-curve",
    )
    args = ap.parse_args()
    out_curve: Path = args.out_curve.resolve()

    _ensure_out_curve(out_curve, args.run)

    manifest = json.loads((out_curve / "manifest.json").read_text())
    if manifest.get("mode") != "hazard_curve":
        raise SystemExit(
            f"Expected manifest mode 'hazard_curve', got {manifest.get('mode')!r}. "
            "Use hazard-curve INI and a separate output directory (see README)."
        )

    branches = sorted(manifest["branches"], key=lambda b: b["global_index"])
    n = len(branches)
    weights = np.asarray([b["combined_branch_weight"] for b in branches], dtype=float)
    if not np.isclose(weights.sum(), 1.0, atol=TOL, rtol=0):
        print(f"WARNING: combined weights sum to {weights.sum()}, expected 1.0", file=sys.stderr)

    curves: list[np.ndarray] = []
    d0_ref: np.ndarray | None = None
    for b in branches:
        cf = b.get("curve_file")
        if not cf:
            raise SystemExit(f"Branch missing curve_file: {b}")
        path = out_curve / cf
        d0, rates = _read_single_site_curve(path)
        if d0_ref is None:
            d0_ref = d0
        elif not np.allclose(d0_ref, d0):
            raise SystemExit("D0 grid mismatch across branch CSVs")
        curves.append(rates)

    rates_stack = np.stack(curves, axis=0)  # (n_branch, n_d0)
    analytical = np.tensordot(weights, rates_stack, axes=(0, 0))

    agg = _read_aggregate(out_curve / "aggregate_hazard.csv")
    if not np.allclose(agg["D0"], d0_ref):
        raise SystemExit("D0 in aggregate_hazard.csv differs from branch CSVs")
    lt_mean = agg["mean"]
    max_err = float(np.max(np.abs(lt_mean - analytical)))

    if max_err > TOL:
        raise SystemExit(
            f"FAIL: max |LT mean - analytical| = {max_err:.3e} "
            f"(tol {TOL:.0e})."
        )
    print(f"PASS: max |LT mean - Σ w_i λ_i| = {max_err:.3e} (tol {TOL:.0e})")

    # ----- plot (optional dependency) -----
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed; skipping plot")
        return

    fig, ax = plt.subplots(figsize=(8.5, 6.0))

    p05, p95 = agg.get("p05"), agg.get("p95")
    p16, p84 = agg.get("p16"), agg.get("p84")
    if p05 is not None and p95 is not None:
        ax.fill_between(
            d0_ref,
            _positive(p05),
            _positive(p95),
            color="tab:red",
            alpha=0.10,
            label="LT 5-95% fractile band",
        )
    if p16 is not None and p84 is not None:
        ax.fill_between(
            d0_ref,
            _positive(p16),
            _positive(p84),
            color="tab:red",
            alpha=0.15,
            label="LT 16-84% fractile band",
        )

    colors = ["tab:blue", "tab:cyan", "tab:orange", "tab:purple"]
    markers = ["o", "v", "s", "^"]
    for i, b in enumerate(branches):
        lab = (
            f'{b["source_model_branch_id"]}  '
            f'(w={b["combined_branch_weight"]:.2f})'
        )
        ax.plot(
            d0_ref,
            _positive(rates_stack[i]),
            color=colors[i % len(colors)],
            lw=1.5,
            marker=markers[i % len(markers)],
            ms=4,
            label=lab,
        )

    ax.plot(
        d0_ref,
        _positive(analytical),
        color="tab:green",
        lw=2.4,
        ls="--",
        label=r"analytical $\sum_i w_i\,\lambda_i(D_0)$",
    )
    ax.plot(
        d0_ref,
        _positive(lt_mean),
        color="tab:red",
        lw=1.2,
        ls=":",
        marker="x",
        ms=7,
        mew=1.6,
        label="LT aggregate mean",
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Displacement threshold $D_0$ (m)")
    ax.set_ylabel(r"Annual rate of exceedance $\lambda(D_0)$ (1/yr)")
    ax.set_title(
        "Source-model logic tree - hazard curve validation\n"
        "(SMLT × bGRRelative × minimal FDHA tree; Pizza2023 PSR + "
        "Petersen2011PrimaryFD bilinear)"
    )
    ax.grid(True, which="both", ls=":", alpha=0.6)
    ax.legend(loc="lower left", fontsize=8, framealpha=0.95)
    ax.annotate(
        (
            r"$\max\,|\,\mathrm{LT\ mean} - \sum_i w_i \lambda_i\,|$"
            f"\n= {max_err:.2e}  (tol $10^{{-12}}$)"
        ),
        xy=(0.98, 0.98),
        xycoords="axes fraction",
        ha="right",
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.6", alpha=0.9),
    )

    fig.tight_layout()
    fig.savefig(args.plot, dpi=150)
    print(f"Wrote {args.plot}")


if __name__ == "__main__":
    main()
