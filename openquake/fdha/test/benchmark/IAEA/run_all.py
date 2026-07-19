#!/usr/bin/env python
"""Run all IAEA PFDHA exercise benchmark jobs and compare with the paper.

For every entry in :mod:`manifest` this script

1. runs the logic-tree job ``<case>/job_<job>.ini``,
2. writes the computed curve to ``computed/<case>_<job>.csv``,
3. compares it (after any documented ``post_factor``) with the digitized
   published curve in ``reference/`` and prints a summary table,
4. writes ``comparison_summary.json``.

Usage (repo root)::

    PYTHONPATH=. python openquake/fdha/test/benchmark/IAEA/run_all.py [job-filter]

An optional substring filter restricts which jobs run, e.g. ``norcia`` or
``distributed_V24``.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from manifest import MANIFEST  # noqa: E402


def load_reference(figure_csv: str, column: str):
    import csv

    with (HERE / "reference" / figure_csv).open() as f:
        rows = list(csv.DictReader(f))
    disp_m = np.array([float(r["disp_cm"]) for r in rows]) / 100.0
    ref = np.array([float(r[column]) for r in rows])
    return disp_m, ref


def run_job(case: str, job: str, outdir: Path):
    from openquake.fdha.logic_tree.driver import FdhaLogicTree

    ini = HERE / case / f"job_{job}.ini"
    result = FdhaLogicTree.from_ini(ini).run(outdir=outdir)
    d0 = np.asarray(result.d0, dtype=float)
    rates = np.asarray(result.mean_rates, dtype=float)
    if rates.ndim == 2:
        rates = rates[0]
    return d0, rates


def compare(entry, d0, rates):
    """Return a comparison dict for one manifest entry."""
    ref_d, ref = load_reference(entry.figure_csv, entry.column)
    computed = rates * entry.post_factor

    if ref_d.shape != d0.shape or not np.allclose(ref_d, d0, rtol=1e-9):
        # different displacement grid (e.g. M11): interpolate computed in log-log
        pos = computed > 0
        computed = np.exp(np.interp(np.log(ref_d), np.log(d0[pos]),
                                    np.log(computed[pos])))

    ok = (ref > 0) & (computed > 0)
    if entry.assert_dmax_m is not None:
        in_range = ok & (ref_d <= entry.assert_dmax_m + 1e-12)
    else:
        in_range = ok
    ratio = np.full_like(ref, np.nan)
    ratio[ok] = computed[ok] / ref[ok]
    with np.errstate(all="ignore"):
        max_re_range = float(np.nanmax(np.abs(ratio[in_range] - 1))) if in_range.any() else float("nan")
        max_re_full = float(np.nanmax(np.abs(ratio[ok] - 1))) if ok.any() else float("nan")

    return {
        "case": entry.case,
        "job": entry.job,
        "model": entry.column,
        "figure": entry.figure_csv,
        "post_factor": entry.post_factor,
        "max_relerr_assert_range": max_re_range,
        "max_relerr_full": max_re_full,
        "assert_max_relerr": entry.assert_max_relerr,
        "assert_dmax_m": entry.assert_dmax_m,
        "ratio": [None if not np.isfinite(x) else round(float(x), 4) for x in ratio],
    }


def main() -> int:
    job_filter = sys.argv[1] if len(sys.argv) > 1 else ""
    computed_dir = HERE / "computed"
    computed_dir.mkdir(exist_ok=True)

    results = []
    for entry in MANIFEST:
        tag = f"{entry.case}_{entry.job}"
        if job_filter and job_filter not in tag:
            continue
        print(f"=== {tag}")
        outdir = HERE / "out" / tag
        d0, rates = run_job(entry.case, entry.job, outdir)
        if not np.any(rates > 0):
            print(f"!!! {tag}: all-zero hazard curve - model chain failed, "
                  f"check logs in {outdir}")
        csv_path = computed_dir / f"{tag}.csv"
        with csv_path.open("w") as f:
            f.write("disp_m,rate,rate_post_factor\n")
            for d, r in zip(d0, rates):
                f.write(f"{d:g},{r:.6e},{r * entry.post_factor:.6e}\n")
        summary = compare(entry, d0, rates)
        results.append(summary)
        print(f"    max|relerr| assert-range={summary['max_relerr_assert_range']:.3%} "
              f"full={summary['max_relerr_full']:.3%} "
              f"(limit={entry.assert_max_relerr}, dmax={entry.assert_dmax_m})")

    # Merge into any existing summary so a filtered rerun updates only the
    # jobs it ran instead of truncating the file.
    summary_path = HERE / "comparison_summary.json"
    merged = {}
    if job_filter and summary_path.exists():
        for r in json.loads(summary_path.read_text()):
            merged[(r["case"], r["job"])] = r
    for r in results:
        merged[(r["case"], r["job"])] = r
    ordered = [merged[(e.case, e.job)] for e in MANIFEST
               if (e.case, e.job) in merged]
    summary_path.write_text(json.dumps(ordered, indent=2))

    print("\n=== SUMMARY " + "=" * 66)
    print(f"{'job':34s} {'model':6s} {'range err':>10s} {'full err':>10s} {'limit':>7s}  status")
    n_fail = 0
    for r in results:
        lim = r["assert_max_relerr"]
        if lim is None:
            status = "qualitative"
        elif r["max_relerr_assert_range"] <= lim:
            status = "PASS"
        else:
            status = "FAIL"
            n_fail += 1
        print(f"{r['case'] + '/' + r['job']:34s} {r['model']:6s} "
              f"{r['max_relerr_assert_range']:>9.1%} {r['max_relerr_full']:>9.1%} "
              f"{('-' if lim is None else format(lim, '.0%')):>7s}  {status}")
    print(f"\nSaved comparison_summary.json ({n_fail} failures)")
    return 1 if n_fail else 0


if __name__ == "__main__":
    # NUMBA_DISABLE_JIT only in the __main__ block (see run_sensitivity.py):
    # setting it at import time breaks any pytest process that imports this
    # module alongside already-jitted hazardlib code.
    os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
    raise SystemExit(main())
