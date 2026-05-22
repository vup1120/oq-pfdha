"""Run the SMLT example end-to-end (curve + map) and print a short summary.

Usage::

    python examples/source_model_logic_tree/run.py

What this script does
---------------------

1. Loads ``fdha_curve.ini``, runs the FDHA logic-tree driver, and prints
   the per-realisation breakdown the driver wrote to ``manifest.json``.
2. Loads ``fdha_map.ini``, runs the same logic tree in hazard-map mode,
   and prints the location of the SMLT-weighted aggregate rate cube.

Both modes use the same ``source_model_logic_tree.xml`` and the same
``fdha_logic_tree.xml``; the only difference between the two INIs is
``[geometry] sites = ...`` vs. ``[geometry] region = ...``.

The driver mirrors OpenQuake's ``get_smlt`` semantics:

* When ``source_model_logic_tree_file`` is set, the SMLT XML is
  authoritative and ``source_model_file`` is ignored.
* When it is absent, ``source_model_file`` is used and the run is
  byte-for-byte identical to a pre-SMLT PFDHA execution.
"""
from __future__ import annotations

import json
from pathlib import Path

from openquake.fdha.logic_tree.driver import FdhaLogicTree


HERE = Path(__file__).resolve().parent


def _print_manifest_summary(out_dir: Path, mode: str) -> None:
    manifest = json.loads((out_dir / "manifest.json").read_text())
    print(f"\n=== {mode} manifest summary ({out_dir}) ===")
    if "source_model_branches" in manifest:
        print(f"  source-model realisations: "
              f"{len(manifest['source_model_branches'])}")
        for sm in manifest["source_model_branches"]:
            print(f"    - {sm['branch_id']:>10s}  weight={sm['weight']:.3f}")
    print(f"  combined PFDHA realisations: {manifest['n_realizations']}")
    print(f"  total combined weight:       "
          f"{manifest['total_combined_weight']:.6f}")
    if manifest.get("branches"):
        print("  per-realisation breakdown:")
        for rec in manifest["branches"]:
            uncs = ", ".join(
                f"{u['uncertainty_type']}={u['value']}"
                for u in rec.get("source_model_uncertainties", [])
            ) or "(none)"
            print(
                f"    [{rec['global_index']:>2d}]  sm={rec['source_model_branch_id']:<12s}"
                f"  fdha={rec['fdha_branch_id']:<10s}"
                f"  weight={rec['combined_branch_weight']:.4f}"
                f"  uncertainties=[{uncs}]"
            )


def run_curve() -> None:
    out = HERE / "out_curve"
    print(f"--- Running hazard-curve mode -> {out}")
    res = FdhaLogicTree.from_ini(HERE / "fdha_curve.ini").run(outdir=out)
    print(f"  mean curve shape: {len(res.mean_rates)} site(s) "
          f"x {len(res.d0)} D0 levels")
    _print_manifest_summary(out, "hazard_curve")
    print(f"  per-realisation CSVs: {out / 'hazard_curves'}")
    print(f"  aggregate CSV:        {out / 'aggregate_hazard.csv'}")


def run_map() -> None:
    out = HERE / "out_map"
    print(f"\n--- Running hazard-map mode -> {out}")
    res = FdhaLogicTree.from_ini(HERE / "fdha_map.ini").run(outdir=out)
    print(f"  mean map shape: {len(res.site_lons or [])} site(s) "
          f"x {len(res.d0)} D0 levels")
    _print_manifest_summary(out, "hazard_map")
    print(f"  per-realisation HDF5: "
          f"{out / 'source_model_branches' / '<idx>_<branch_id>' / 'branches'}")
    print(f"  aggregate HDF5:       {out / 'aggregate' / 'rates_mean.h5'}")
    print(f"  displacement maps:    {out / 'aggregate'}/displacement_map_*.csv")


if __name__ == "__main__":
    run_curve()
    run_map()
