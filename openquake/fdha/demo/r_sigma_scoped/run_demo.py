"""Demo: scoping the r_sigma_km epistemic branches to SPECIFIC sources.

A ``fdhaCalcRSigma`` branch set accepts the standard ``applyToSources``
filter, so the mapping-accuracy uncertainty can be propagated for one or two
faults while every other fault keeps the default sigma = 0 boxcar
(half-width ``r_threshold_km``). That is the physically natural use: mapping
accuracy is a property of EACH fault's trace (one fault well surveyed,
another concealed), not of the whole model.

Setup: TWO parallel M7.0 characteristic faults ~0.69 km apart
(``source_model_2src.xml``) and two sites, each 40 m from "its" fault:

  site 1 - 40 m from fault A, ~0.65 km from fault B;
  site 2 - 40 m from fault B, ~0.65 km from fault A.

The separation exceeds both near-field scales (boxcar h = 50 m, widest
Gaussian support 2 sigma = 131 m), so each fault's sigma treatment can only
matter at its own site - which is exactly what the assertions check.

Three runs of the SAME two-source job, differing only in the FDHA tree:

  baseline - no sigma branch set: both faults on the boxcar (default);
  scoped   - sigma set (Accurate w=0.5, Concealed w=0.5) with
             applyToSources="src_A": fault A Gaussian, fault B boxcar;
  unscoped - the same sigma set WITHOUT applyToSources: both faults
             Gaussian.

Assertions (machine precision where exact):

  1. scoped == baseline at site 2  -> fault B genuinely stayed on the
     boxcar inside the scoped job (and fault A's far-field contribution
     there is sigma-independent by construction);
  2. scoped == unscoped at site 1  -> fault A's branches are the same
     whether or not the set is scoped;
  3. scoped != baseline at site 1  -> the scoping actually did something;
  4. unscoped != baseline at site 2 -> without the filter, fault B's site
     would have changed too.

Validator FDLT-012 guarantees each source is covered by AT MOST one sigma
branch set (disjoint applyToSources scopes may coexist; an unscoped set
covers every source and cannot coexist with another set).

All sigma values are Petersen et al. (2011, Tables 2-3) two-sided
mapping-accuracy classes, used here as demo inputs, not recommendations.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np

from openquake.fdha.logic_tree.driver import FdhaLogicTree

HERE = Path(__file__).resolve().parent
BASE_DEMO = HERE.parents[3] / "examples"
LOCAL_SRC = HERE / "source_model_2src.xml"
SRC_XML = "source_model.xml"                    # name the SMLT expects
SMLT_XML = "hazard_curve_minimal_source_model_logic_tree.xml"
FDHA_XML = "hazard_curve_minimal_fdha_logic_tree.xml"
JOB_INI = "hazard_curve_minimal.ini"
OUT = HERE / "out"

# Site 1 is 40 m from fault A (same point as site A of the sibling
# r_sigma_epistemic demo); site 2 is the same point shifted east by the
# fault-B offset (0.008 deg), i.e. 40 m from fault B.
SITE_1 = (16.16186114, 39.64824366)             # 40 m from A
SITE_2 = (16.16986114, 39.64824366)             # 40 m from B
H_BOXCAR = 0.05                                 # r_threshold_km in the INI

# Two Petersen mapping-accuracy classes as the demo sigma branches.
SIGMA_BRANCHES = [
    ("RS_ACCURATE",  "0.02689", "0.5"),
    ("RS_CONCEALED", "0.06552", "0.5"),
]


def sigma_level(branches, apply_to_sources=None):
    scope = (f' applyToSources="{apply_to_sources}"'
             if apply_to_sources else "")
    rows = [
        '    <logicTreeBranchingLevel branchingLevelID="bl_5_r_sigma">',
        f'      <logicTreeBranchSet branchSetID="bs_5_r_sigma" '
        f'uncertaintyType="fdhaCalcRSigma"{scope}>',
    ]
    for bid, value, weight in branches:
        rows += [
            f'        <logicTreeBranch branchID="{bid}">',
            f'          <uncertaintyModel>{value}</uncertaintyModel>',
            f'          <uncertaintyWeight>{weight}</uncertaintyWeight>',
            '        </logicTreeBranch>',
        ]
    rows += ['      </logicTreeBranchSet>', '    </logicTreeBranchingLevel>']
    return "\n".join(rows)


def make_job(name: str, branches=None, apply_to_sources=None) -> Path:
    """Clone examples/hazard_curve_minimal into out/<name>/ with the local
    two-source model, the two demo sites, and an optional sigma branch set
    (optionally scoped via applyToSources)."""
    workdir = OUT / name
    workdir.mkdir(parents=True, exist_ok=True)
    shutil.copy(LOCAL_SRC, workdir / SRC_XML)
    shutil.copy(BASE_DEMO / SMLT_XML, workdir)

    lt_text = (BASE_DEMO / FDHA_XML).read_text()
    if branches:
        lt_text = lt_text.replace(
            "  </logicTree>",
            sigma_level(branches, apply_to_sources) + "\n  </logicTree>",
        )
    (workdir / FDHA_XML).write_text(lt_text)

    ini_text = (BASE_DEMO / JOB_INI).read_text()
    ini_text = ini_text.replace(
        "r_threshold_km = 0.1", f"r_threshold_km = {H_BOXCAR}"
    )
    ini_text = ini_text.replace(
        "sites = 16.16573727 39.64704451",
        f"sites = {SITE_1[0]} {SITE_1[1]}, {SITE_2[0]} {SITE_2[1]}",
    )
    ini = workdir / JOB_INI
    ini.write_text(ini_text)
    return ini


def run(name: str, branches=None, apply_to_sources=None):
    ini = make_job(name, branches, apply_to_sources)
    res = FdhaLogicTree.from_ini(str(ini)).run(outdir=OUT / name / "out")
    manifest = json.loads((OUT / name / "out" / "manifest.json").read_text())
    return {
        "d0": np.asarray(res.d0, dtype=float),
        "mean": np.asarray(res.mean_rates, dtype=float),   # (2 sites, D)
        "manifest": manifest,
    }


def rel_diff(a, b):
    return np.abs(a - b) / np.where(b != 0, np.abs(b), 1.0)


def main():
    if OUT.exists():
        shutil.rmtree(OUT)

    print("Running baseline (no sigma set: both faults boxcar)...")
    baseline = run("baseline")
    print("Running scoped   (sigma set applyToSources='src_A')...")
    scoped = run("scoped", SIGMA_BRANCHES, apply_to_sources="src_A")
    print("Running unscoped (same sigma set, no filter: both faults)...")
    unscoped = run("unscoped", SIGMA_BRANCHES)

    d0 = baseline["d0"]

    # 1. Fault B untouched by the scoped set: site 2 identical to baseline.
    e1 = rel_diff(scoped["mean"][1], baseline["mean"][1]).max()
    # 2. Fault A's branches identical scoped vs unscoped: site 1 identical.
    e2 = rel_diff(scoped["mean"][0], unscoped["mean"][0]).max()
    # 3./4. The sigma treatment genuinely changes the near-fault site.
    e3 = rel_diff(scoped["mean"][0], baseline["mean"][0]).max()
    e4 = rel_diff(unscoped["mean"][1], baseline["mean"][1]).max()

    print()
    print(f"scoped   vs baseline @ site 2 (near B): max rel diff = {e1:.3e}"
          "  <- B stayed boxcar")
    print(f"scoped   vs unscoped @ site 1 (near A): max rel diff = {e2:.3e}"
          "  <- A same branches")
    print(f"scoped   vs baseline @ site 1 (near A): max rel diff = {e3:.3e}"
          "  <- scoping acted on A")
    print(f"unscoped vs baseline @ site 2 (near B): max rel diff = {e4:.3e}"
          "  <- unscoped would hit B")
    assert e1 == 0.0, "fault B must be bit-identical to the boxcar baseline"
    # Scoped and unscoped enumerate DIFFERENT realization sets (A-branches x
    # single-B vs branches on both sources), so the weighted aggregation
    # accumulates in a different order; the site-1 rates agree to the last
    # ulp (~2e-16 observed) rather than bit-for-bit.
    assert e2 < 1e-14, "fault A must not care whether the set is scoped"
    assert e3 > 0.1, "sigma branches must genuinely change site 1"
    assert e4 > 0.1, "the unscoped set must genuinely change site 2"

    # The manifest lists realizations PER SOURCE: fault A carries the two
    # sigma branches, fault B a single sigma-less branch. Independent
    # sources are summed, not weighted-averaged, so the weights normalise
    # within each source group (not globally).
    print("\nscoped-run realizations (out/scoped/out/manifest.json):")
    per_source: dict = {}
    for b in scoped["manifest"]["branches"]:
        src = b["fdha_source_id"]
        sig = b.get("fdha_calc_params", {}).get("r_sigma_km", "(boxcar, 0)")
        print(f"  {src:<8} {b['fdha_branch_id'].split('|')[0]:<24.24} "
              f"r_sigma_km={sig!s:<12} w={b['combined_branch_weight']}")
        per_source.setdefault(src, 0.0)
        per_source[src] += b["combined_branch_weight"]
    for src, w_sum in per_source.items():
        assert abs(w_sum - 1.0) < 1e-12, \
            f"weights of source {src} must sum to 1 (got {w_sum})"
    assert set(per_source) == {"src_A", "src_B"}

    # ------------------------------------------------------------- figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), sharey=True)
    for i, (ax, label, fault, note) in enumerate((
        (axes[0], "1", "A", "sigma-scoped fault:\nscoped == unscoped, both differ from boxcar"),
        (axes[1], "2", "B", "excluded fault:\nscoped == boxcar baseline (curves coincide)"),
    )):
        ax.loglog(d0, baseline["mean"][i], color="C0", lw=3.6,
                  solid_capstyle="round", zorder=2,
                  label="baseline: both faults boxcar ($\\sigma$=0)")
        ax.loglog(d0, unscoped["mean"][i], color="#666666", lw=2.0, ls=":",
                  zorder=3, label="unscoped: $\\sigma$ set on BOTH faults")
        ax.loglog(d0, scoped["mean"][i], color="#b2182b", lw=2.0, ls="--",
                  zorder=4,
                  label='scoped: $\\sigma$ set applyToSources="src_A"')
        ax.set_xlabel("Displacement $D_0$ [m]")
        ax.set_title(f"Site {label}: 40 m from fault {fault}", fontsize=10.5)
        ax.text(0.03, 0.05, note, transform=ax.transAxes, fontsize=8.5,
                va="bottom", style="italic", color="#333333")
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(fontsize=8, loc="upper right")
    axes[0].set_ylabel("Annual rate of exceedance")
    fig.suptitle(
        "fdhaCalcRSigma scoped to one source: fault A propagates the "
        "mapping-accuracy branches, fault B keeps the boxcar",
        fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    png = OUT / "r_sigma_scoped_demo.png"
    fig.savefig(png, dpi=150)
    print(f"\nFigure written to {png}")


if __name__ == "__main__":
    main()
