#!/usr/bin/env python
"""Run the IAEA sensitivity-case demonstration jobs (see sensitivity.py).

Writes ``computed/<case>_<job>.csv`` for every job and compares against
the committed regression snapshots in ``reference_snapshots/``.

Usage (repo root)::

    PYTHONPATH=. python openquake/fdha/test/benchmark/IAEA/run_sensitivity.py \
        [job-filter] [--update-snapshots]

``--update-snapshots`` rewrites the snapshots from the current run (do
this only for intentional model changes, and review the diff).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from sensitivity import SENSITIVITY_JOBS  # noqa: E402


def run_job(case, job, outdir):
    from openquake.fdha.logic_tree.driver import FdhaLogicTree

    ini = HERE / case / f"job_{job}.ini"
    result = FdhaLogicTree.from_ini(ini).run(outdir=outdir)
    d0 = np.asarray(result.d0, dtype=float)
    rates = np.asarray(result.mean_rates, dtype=float)
    if rates.ndim == 2:
        rates = rates[0]
    return d0, rates


def load_snapshot(tag):
    path = HERE / "reference_snapshots" / f"{tag}.csv"
    if not path.exists():
        return None, None
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return data[:, 0], data[:, 1]


def main() -> int:
    update = "--update-snapshots" in sys.argv
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    job_filter = args[0] if args else ""

    (HERE / "computed").mkdir(exist_ok=True)
    (HERE / "reference_snapshots").mkdir(exist_ok=True)

    n_fail = 0
    for sj in SENSITIVITY_JOBS:
        tag = f"{sj.case}_{sj.job}"
        if job_filter and job_filter not in tag:
            continue
        print(f"=== {tag}: {sj.description}")
        d0, rates = run_job(sj.case, sj.job, HERE / "out" / tag)

        with (HERE / "computed" / f"{tag}.csv").open("w") as f:
            f.write("disp_m,rate\n")
            for d, r in zip(d0, rates):
                f.write(f"{d:g},{r:.6e}\n")

        if update:
            with (HERE / "reference_snapshots" / f"{tag}.csv").open("w") as f:
                f.write("disp_m,rate\n")
                for d, r in zip(d0, rates):
                    f.write(f"{d:g},{r:.6e}\n")
            print(f"    snapshot updated (head {rates[0]:.4e})")
            continue

        snap_d, snap_r = load_snapshot(tag)
        if snap_d is None:
            print("    !!! no snapshot - run with --update-snapshots")
            n_fail += 1
            continue
        ok = (np.allclose(snap_d, d0, rtol=1e-9)
              and np.allclose(snap_r, rates, rtol=1e-4, atol=1e-30))
        status = "PASS" if ok else "FAIL"
        n_fail += not ok
        print(f"    head {rates[0]:.4e}  snapshot {snap_r[0]:.4e}  {status}")

    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
