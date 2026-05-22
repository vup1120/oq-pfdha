#!/usr/bin/env python
"""Run the Norcia Case 3 IAEA curve and map benchmark."""

from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

from openquake.fdha.logic_tree.driver import FdhaLogicTree


HERE = Path(__file__).resolve().parent
CURVE_INI = HERE / "job_norcia_case3_iaea_curve.ini"
MAP_INI = HERE / "job_norcia_case3_iaea_map.ini"


def _run_one(ini: Path, outdir: Path) -> dict:
    result = FdhaLogicTree.from_ini(ini).run(outdir=outdir)
    summary = {
        "mode": result.mode,
        "outdir": result.outdir,
        "n_displ": len(result.d0),
        "n_sites": len(result.site_lons or []),
        "manifest_json": str(outdir / "manifest.json"),
    }
    if result.mode == "hazard_curve":
        summary["aggregate_hazard_csv"] = str(outdir / "aggregate_hazard.csv")
        summary["mean_rates_site0"] = result.mean_rates[0]
    if result.mode == "hazard_map":
        summary["rates_mean_h5"] = str(outdir / "aggregate" / "rates_mean.h5")
        summary["displacement_map_mean_csv"] = str(
            outdir / "aggregate" / "displacement_map_mean.csv"
        )
        summary["target_return_period"] = result.target_return_period
    return summary


def main() -> int:
    curve = _run_one(CURVE_INI, HERE / "out_curve")
    hazard_map = _run_one(MAP_INI, HERE / "out_map")
    results = {"curve": curve, "map": hazard_map}
    results_path = HERE / "results.json"
    results_path.write_text(json.dumps(results, indent=2))
    print(f"Saved benchmark summary to {results_path}")
    print(f"Curve aggregate: {curve['aggregate_hazard_csv']}")
    print(f"Map aggregate: {hazard_map['displacement_map_mean_csv']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
