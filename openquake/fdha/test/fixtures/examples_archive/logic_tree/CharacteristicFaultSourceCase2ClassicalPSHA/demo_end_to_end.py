"""End-to-end demonstration of the FDHA logic tree.

Runs ``job_SS_only.ini`` (runnable demo), then ``job.ini`` (template that
must halt at FDLT-007), and checks a battery of invariants on the outputs.

Not part of the automated suite; this is a human-readable transcript the
user can inspect to see the logic tree behaving as the v3 spec demands.
"""
from __future__ import annotations

import json
import math
import sys
import warnings
from pathlib import Path


HERE = Path(__file__).resolve().parent


def _hr(title: str) -> None:
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def _check(label: str, ok: bool, detail: str = "") -> None:
    status = "PASS" if ok else "FAIL"
    line = f"  [{status}] {label}"
    if detail:
        line += f"  -- {detail}"
    print(line)
    if not ok:
        raise AssertionError(f"Invariant failed: {label} ({detail})")


def demo_runnable_ss_only() -> Path:
    ini = HERE / "job_SS_only.ini"
    outdir = HERE / "demo_out"

    _hr("STEP 1 | Parse the runnable demo INI")
    from openquake.fdha.calc.config_loader import load_config
    cfg = load_config(str(ini))
    print("  ini                   :", ini.name)
    print("  [geometry].sites      :", cfg["geometry"]["sites"])
    print(
        "  [calculation].source_model_logic_tree_file :",
        cfg["calculation"]["source_model_logic_tree_file"],
    )
    print(
        "  [calculation].fdha_logic_tree_file         :",
        cfg["calculation"]["fdha_logic_tree_file"],
    )
    has_legacy_models = bool(cfg.get("models"))
    _check("No legacy [models.*] section present (LT-only path)", not has_legacy_models,
           f"cfg['models'] = {cfg.get('models')!r}")

    _hr("STEP 2 | Parse + validate the NRML XML (no calculator yet)")
    from openquake.fdha.logic_tree.nrml_reader import parse as parse_nrml
    from openquake.fdha.logic_tree.validators import validate_spec
    spec = parse_nrml(HERE / "fdha_logic_tree_SS.xml")
    n_bsets = sum(len(lv.branch_sets) for lv in spec.branching_levels)
    print(f"  branching_levels      : {len(spec.branching_levels)}")
    print(f"  branch_sets           : {n_bsets}")
    report = validate_spec(spec, source_ids={"3"})
    print(f"  mandatory errors      : {len(report.errors)}  (expect 0)")
    print(f"  advisory warnings     : {len(report.warnings)}")
    for w in report.warnings:
        print(f"     - {w.code}: {w.message}")
    _check("No FDLT-001..007 errors on the SS spec", len(report.errors) == 0)
    _check("At least one FDLT-101 advisory fires (Petersen D_p_L + Chiou D_sp_Nstar)",
           any(w.code == "FDLT-101" for w in report.warnings))
    _check("FDLT-101 message carries the cited source string",
           any("Sarmiento et al. 2025" in w.message for w in report.warnings))

    _hr("STEP 3 | Run FdhaLogicTree end-to-end")
    from openquake.fdha.logic_tree.driver import FdhaLogicTree
    if outdir.exists():
        for p in sorted(outdir.rglob("*"), reverse=True):
            if p.is_file():
                p.unlink()
            elif p.is_dir():
                p.rmdir()
        outdir.rmdir()
    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter("always")
        lt = FdhaLogicTree.from_ini(str(ini))
        result = lt.run(outdir=outdir)
    dep = [w for w in wl if issubclass(w.category, DeprecationWarning)]
    print(f"  DeprecationWarnings   : {len(dep)}   (expect 0: we're on the LT path)")
    print(f"  outdir                : {result.outdir}")
    print(f"  n_D0                  : {len(result.d0)}")
    print(f"  D0 values             : {[round(x, 4) for x in result.d0]}")
    print(f"  fractile keys         : {sorted(result.fractiles.keys())}")
    print(f"  mean shape            : {len(result.mean_rates)} x {len(result.mean_rates[0])}")
    _check("No DeprecationWarning on LT path", len(dep) == 0)

    return outdir


def inspect_outputs(outdir: Path) -> None:
    _hr("STEP 4 | Inspect output directory layout")
    expected = {
        "aggregate_hazard.csv",
        "manifest.json",
        "validator_report.txt",
        "hazard_curves",
        "branch_configs",
    }
    got = {p.name for p in outdir.iterdir()}
    print(f"  got                   : {sorted(got)}")
    _check("All expected artefacts present", expected.issubset(got),
           f"missing = {expected - got}")

    curves = sorted((outdir / "hazard_curves").glob("branch_*.csv"))
    print(f"  per-branch curve files: {len(curves)}  (one per end-branch)")

    manifest = json.loads((outdir / "manifest.json").read_text())
    n_manifest = len(manifest["branches"])
    print(f"  manifest.branches     : {n_manifest}")
    _check("Per-branch curve count matches manifest", len(curves) == n_manifest,
           f"{len(curves)} vs {n_manifest}")

    required_keys = {"index", "fingerprint", "weight", "source_id", "style", "models", "curve_file"}
    required_slots = {"primary_surf_rup", "primary_surf_displ",
                      "secondary_surf_rup", "secondary_surf_displ"}
    for b in manifest["branches"]:
        missing = required_keys - b.keys()
        _check(f"branch#{b['index']:02d} has all required manifest keys",
               not missing, f"missing={missing}")
        slots_missing = required_slots - set(b["models"].keys())
        _check(f"branch#{b['index']:02d} names all four calculator slots",
               not slots_missing, f"missing={slots_missing}")

    adv = manifest.get("advisory_warnings", [])
    print(f"  advisory warnings     : {len(adv)}")
    for a in adv:
        print(f"     - {a['code']}: {a['message']}")

    _hr("STEP 5 | Show first three end-branches (models + weights)")
    for b in manifest["branches"][:3]:
        print(f"  branch #{b['index']:02d}  fp={b['fingerprint']}  w={b['weight']:.6f}  style={b['style']}")
        for slot, cls in b["models"].items():
            print(f"            {slot:22s} = {cls}")

    _hr("STEP 6 | Weight-sum invariant")
    w_total = sum(b["weight"] for b in manifest["branches"])
    print(f"  Σ weights             : {w_total:.15f}")
    _check("Σ weights == 1.0 ± 1e-12", abs(w_total - 1.0) < 1e-12,
           f"delta={abs(w_total-1.0):.3e}")


def rate_invariants(outdir: Path) -> None:
    import csv
    _hr("STEP 7 | Rate-shape invariants on every per-branch curve")

    def _read_branch_csv(p: Path):
        d0, lam = [], []
        with p.open() as f:
            r = csv.DictReader(f)
            for row in r:
                d0.append(float(row["D0"]))
                lam.append(float(row["annual_rate"]))
        return d0, lam

    curves = sorted((outdir / "hazard_curves").glob("branch_*.csv"))
    worst = 0.0
    for p in curves:
        d0, lam = _read_branch_csv(p)
        for i in range(len(lam) - 1):
            diff = lam[i] - lam[i + 1]
            worst = max(worst, -diff)
            # allow tiny numerical noise
            _check(f"{p.name} monotone non-increasing at D0={d0[i]:.4g}->{d0[i+1]:.4g}",
                   diff >= -1e-15,
                   f"lam[i]={lam[i]:.6e} lam[i+1]={lam[i+1]:.6e}")
    print(f"  worst (lam[i+1]-lam[i]) across all curves = {worst:.3e} "
          "(should be ≤ 0, ≥ -1e-15)")

    _hr("STEP 8 | Aggregation invariants (mean = Σ w·λ, fractiles ordered, bracket mean)")
    manifest = json.loads((outdir / "manifest.json").read_text())
    branches = manifest["branches"]
    # Reconstruct per-branch rate matrix in the same branch order as manifest
    per_branch = []
    weights = []
    d0_ref = None
    for b in branches:
        p = outdir / b["curve_file"]
        d0, lam = _read_branch_csv(p)
        if d0_ref is None:
            d0_ref = d0
        else:
            _check(f"D0 grid consistent across branches for {b['curve_file']}",
                   d0 == d0_ref, f"mismatch")
        per_branch.append(lam)
        weights.append(b["weight"])

    # Aggregate CSV
    agg_d0, mean_from_csv, fr_from_csv = [], [], {}
    qs = ["p05", "p16", "p50", "p84", "p95"]
    with (outdir / "aggregate_hazard.csv").open() as f:
        r = csv.DictReader(f)
        cols = r.fieldnames
        for row in r:
            agg_d0.append(float(row["D0"]))
            mean_from_csv.append(float(row["mean"]))
            for q in qs:
                fr_from_csv.setdefault(q, []).append(float(row[q]))
    _check("aggregate_hazard.csv columns",
           set(["D0", "mean", *qs]).issubset(set(cols)),
           f"cols={cols}")
    _check("Aggregate D0 grid == branch D0 grid", agg_d0 == d0_ref)

    worst_mean_err = 0.0
    for j in range(len(agg_d0)):
        analytic = sum(weights[b] * per_branch[b][j] for b in range(len(branches)))
        err = abs(analytic - mean_from_csv[j])
        worst_mean_err = max(worst_mean_err, err)
        _check(f"mean == Σ w·λ at D0={agg_d0[j]:.4g}",
               err < 1e-12, f"|diff|={err:.3e}")
    print(f"  max |Σw·λ − mean_from_csv| = {worst_mean_err:.3e} (tol 1e-12)")

    for j in range(len(agg_d0)):
        pvals = [fr_from_csv[q][j] for q in qs]
        for a, b in zip(pvals, pvals[1:]):
            _check(f"fractile monotonicity at D0={agg_d0[j]:.4g}",
                   a <= b + 1e-18, f"p{a} > p{b}")
        _check(f"p05 ≤ mean ≤ p95 at D0={agg_d0[j]:.4g}",
               pvals[0] - 1e-18 <= mean_from_csv[j] <= pvals[-1] + 1e-18,
               f"p05={pvals[0]:.3e} mean={mean_from_csv[j]:.3e} p95={pvals[-1]:.3e}")
    print("  fractiles ordered and mean bracketed by [p05, p95] at every D0")


def demo_template_halts() -> None:
    _hr("STEP 9 | Template job.ini must halt at FDLT-007 (no calculator runs)")

    import openquake.fdha.logic_tree.driver as driver_mod
    from openquake.fdha.logic_tree.driver import FdhaLogicTree
    from openquake.fdha.logic_tree.types import LogicTreeValidationError

    calls = {"n": 0}
    orig_run = driver_mod._run_single

    def tripwire(*a, **kw):  # pragma: no cover - should never fire
        calls["n"] += 1
        raise AssertionError("TRIPWIRE: calculator must not run before FDLT-007 halt")

    driver_mod._run_single = tripwire
    try:
        try:
            FdhaLogicTree.from_ini(str(HERE / "job.ini")).run(outdir=HERE / "demo_template_halt")
            raise AssertionError("Expected LogicTreeValidationError; template did not halt")
        except LogicTreeValidationError as exc:
            msg = str(exc)
    finally:
        driver_mod._run_single = orig_run

    print("  LogicTreeValidationError raised as expected.")
    _check("Error mentions FDLT-007", "FDLT-007" in msg, msg.splitlines()[0])
    _check("Error references PLACEHOLDER_WEIGHT", "PLACEHOLDER_WEIGHT" in msg, msg.splitlines()[0])
    _check("Tripwire confirms no calculator ran", calls["n"] == 0,
           f"calculator call count = {calls['n']}")
    # Show only the first FDLT-007 line for brevity.
    first = next((ln for ln in msg.splitlines() if "FDLT-007" in ln), msg.splitlines()[0])
    print(f"  first FDLT-007 line   : {first}")


def main() -> None:
    outdir = demo_runnable_ss_only()
    inspect_outputs(outdir)
    rate_invariants(outdir)
    demo_template_halts()
    print()
    print("=" * 78)
    print("ALL INVARIANTS HELD. Logic tree end-to-end demo complete.")
    print("=" * 78)


if __name__ == "__main__":
    sys.exit(main())
