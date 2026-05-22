#!/usr/bin/env python
"""Lightweight consistency checks for the Norcia Case 3 IAEA benchmark."""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path

os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

from openquake.fdha.calc.utils.parsing import parse_source_model_faults
from openquake.fdha.logic_tree.driver import FdhaLogicTree


HERE = Path(__file__).resolve().parent
SOURCE_MODEL = HERE / "source_model_norcia_case3.xml"
CURVE_INI = HERE / "job_norcia_case3_iaea_curve.ini"
OUT_CURVE = HERE / "out_curve"
CURVE_SITES = ((13.278, 42.767), (13.188, 42.749), (13.212, 42.853))
SOURCE_ID = "MVFS,NFS"


def _source_rates() -> dict[str, float]:
    srcs = parse_source_model_faults(
        [str(SOURCE_MODEL)],
        hdf5path="",
        rupture_mesh_spacing=0.5,
        width_of_mfd_bin=0.1,
    )
    return {
        sid: sum(float(rup.occurrence_rate) for rup in src.iter_ruptures())
        for sid, src in srcs.items()
    }


def _ensure_curve() -> None:
    path = OUT_CURVE / "aggregate_hazard.csv"
    manifest = OUT_CURVE / "manifest.json"
    if (
        not path.is_file()
        or not manifest.is_file()
        or not _curve_has_expected_sites(path)
        or not _curve_has_expected_sources(manifest)
    ):
        FdhaLogicTree.from_ini(CURVE_INI).run(outdir=OUT_CURVE)


def _curve_has_expected_sites(path: Path) -> bool:
    try:
        with path.open(newline="") as f:
            reader = csv.DictReader(f)
            if not reader.fieldnames or "site_id" not in reader.fieldnames:
                return False
            seen = {}
            for row in reader:
                sid = int(row["site_id"])
                seen.setdefault(sid, (float(row["lon"]), float(row["lat"])))
        if len(seen) != len(CURVE_SITES):
            return False
        return all(
            math.isclose(seen[i][0], lon) and math.isclose(seen[i][1], lat)
            for i, (lon, lat) in enumerate(CURVE_SITES)
        )
    except (KeyError, OSError, TypeError, ValueError):
        return False


def _load_curve_rates() -> list[float]:
    with (OUT_CURVE / "aggregate_hazard.csv").open(newline="") as f:
        return [float(row["mean"]) for row in csv.DictReader(f)]


def _curve_has_expected_sources(path: Path) -> bool:
    try:
        manifest = json.loads(path.read_text())
        return all(
            branch.get("fdha_source_id") == SOURCE_ID
            for branch in manifest.get("branches", [])
        )
    except (OSError, TypeError, ValueError):
        return False


def main() -> int:
    rates = _source_rates()
    checks = [
        ("MVFS total rate", rates["MVFS"], 4.03281e-4, 5e-10),
        ("NFS total rate", rates["NFS"], 5.407072e-3, 5e-9),
    ]
    for label, got, expected, tol in checks:
        if not math.isclose(got, expected, rel_tol=0.0, abs_tol=tol):
            raise SystemExit(f"{label} mismatch: got {got:.8e}, expected {expected:.8e}")

    _ensure_curve()
    manifest = json.loads((OUT_CURVE / "manifest.json").read_text())
    if len(manifest["branches"]) != 2:
        raise SystemExit(f"Expected 2 FDHA end-branches, found {len(manifest['branches'])}")

    curve_rates = _load_curve_rates()
    if not curve_rates or max(curve_rates) <= 0:
        raise SystemExit("Aggregate curve is empty or non-positive.")

    print("Norcia Case 3 checks passed.")
    print(f"  MVFS rate:          {rates['MVFS']:.6e} /yr")
    print(f"  NFS rate:           {rates['NFS']:.6e} /yr")
    print(f"  FDHA branches:      {len(manifest['branches'])}")
    print(f"  Curve max rate:     {max(curve_rates):.6e} /yr")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
