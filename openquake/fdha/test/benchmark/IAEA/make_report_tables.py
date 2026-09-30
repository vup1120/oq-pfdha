#!/usr/bin/env python
"""Emit LaTeX exceedance-rate tables for the SIGMA3 deliverable.

One table per IAEA exercise: published team value against the framework
value at fixed displacement levels, plus the ratio. Values come from the
same files the benchmark asserts on -- ``reference/`` (coordinator
vectors and digitised curves) and ``computed/`` (written by
``run_all.py``) -- so the report cannot drift from the suite.

Usage (repo root)::

    PYTHONPATH=. python openquake/fdha/test/benchmark/IAEA/make_report_tables.py [outdir]
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LEVELS_M = [0.01, 0.05, 0.1, 0.5, 1.0]          # report grid (m)

# exercise -> (caption, reference csv, {model: computed csv})
EXERCISES = [
    ("kumamoto_p1", "Kumamoto principal, option~1 (fixed Uto rupture)",
     "fig4a_kumamoto_principal.csv",
     {"P11": "kumamoto_principal_P11.csv", "C24": "kumamoto_principal_C24.csv",
      "K24": "kumamoto_principal_K24.csv", "T13": "kumamoto_principal_T13.csv",
      "L23": "kumamoto_principal_L23.csv"}),
    ("kumamoto_p2", "Kumamoto principal, option~2 (floating ruptures)",
     "fig4b_kumamoto_principal_floating.csv",
     {"K24": "kumamoto_floating_K24.csv", "T13": "kumamoto_floating_T13.csv",
      "L23": "kumamoto_floating_L23.csv"}),
    ("kumamoto_d1", "Kumamoto 1st distributed exercise ($r = 5.2$~km)",
     "fig6a_kumamoto_distributed.csv",
     {"P11": "kumamoto_distributed_P11.csv",
      "T13": "kumamoto_distributed_T13.csv"}),
    ("kumamoto_d2", "Kumamoto 2nd distributed exercise ($r = 10$~km)",
     "fig15c_kumamoto_distributed_r10.csv",
     {"P11": "kumamoto_distributed_P11_sens4.csv",
      "T13": "kumamoto_distributed_T13_sens4.csv"}),
    ("kumamoto_d3", "Kumamoto 3rd distributed exercise (Suizenji, $r = 0.6$~km)",
     "fig15d_kumamoto_distributed_suizenji.csv",
     {"T13": "kumamoto_distributed_T13_sens1.csv"}),
    ("leteil_p", "Le Teil principal exercise",
     "fig4c_leteil_principal.csv",
     {"K24": "le_teil_principal_K24.csv", "T13": "le_teil_principal_T13.csv",
      "L23": "le_teil_principal_L23.csv"}),
    ("leteil_d1", "Le Teil 1st distributed exercise ($r = 0.6$~km, footwall)",
     "fig6b_leteil_distributed.csv",
     {"T13": "le_teil_distributed_T13.csv"}),
    ("leteil_d2", "Le Teil 2nd distributed exercise (dip $45^\\circ$~NW, hanging wall)",
     "fig19c_leteil_distributed_dipNW.csv",
     {"T13": "le_teil_distributed_T13_dipflip.csv"}),
    ("leteil_d3", "Le Teil 3rd distributed exercise (three faults)",
     "fig19d_leteil_distributed_3faults.csv",
     {"T13": "le_teil_distributed_T13_sens4.csv"}),
    ("norcia_p", "Norcia principal exercise",
     "fig4d_norcia_principal.csv",
     {"Y03": "norcia_principal_Y03.csv", "K24": "norcia_principal_K24.csv",
      "L23": "norcia_principal_L23.csv"}),
    ("norcia_d1", "Norcia 1st distributed exercise (site MS, $r = 7.6$~km)",
     "fig6c_norcia_distributed.csv",
     {"Y03": "norcia_distributed_Y03.csv"}),
    ("norcia_d2", "Norcia 2nd distributed exercise (site SL, $r = 2.4$~km)",
     "fig22b_norcia_distributed_sl.csv",
     {"Y03": "norcia_distributed_Y03_sens2.csv"}),
    ("norcia_d3", "Norcia 3rd distributed exercise (site MS; MVF and NF)",
     "fig22c_norcia_distributed_2sources.csv",
     {"Y03": "norcia_distributed_Y03_sens3.csv"}),
]


def read_reference(fig_csv, column):
    rows = [r for r in csv.DictReader((HERE / "reference" / fig_csv).open())
            if r.get(column, "").strip()]
    d = np.array([float(r["disp_cm"]) for r in rows]) / 100.0
    v = np.array([float(r[column]) for r in rows])
    return d, v


def read_computed(job_csv):
    rows = list(csv.DictReader((HERE / "computed" / job_csv).open()))
    col = "rate_post_factor" if "rate_post_factor" in rows[0] else "rate"
    return (np.array([float(r["disp_m"]) for r in rows]),
            np.array([float(r[col]) for r in rows]))


def fmt(v):
    if v is None or not np.isfinite(v) or v <= 0:
        return "---"
    m, e = f"{v:.2e}".split("e")
    return f"${m}\\!\\times\\!10^{{{int(e)}}}$"


def loglog_at(d, v, x):
    ok = v > 0
    if ok.sum() < 2 or x < d[ok].min() or x > d[ok].max():
        return None
    return float(10.0 ** np.interp(np.log10(x), np.log10(d[ok]),
                                   np.log10(v[ok])))


def make_table(key, caption, fig_csv, jobs):
    lines = [
        "\\begin{table}[H]", "\\centering", "\\small",
        f"\\caption{{{caption}: annual frequency of exceedance (a$^{{-1}}$) "
        "at fixed displacement levels. Published team value (pub.) against "
        "the framework value (comp.), with their ratio.}",
        f"\\label{{tab:rates:{key}}}",
        "\\begin{tabular}{l l " + " ".join(["r"] * len(LEVELS_M)) + "}",
        "\\toprule",
        "\\textbf{Chain} & & "
        + " & ".join(f"\\textbf{{{L:g}~m}}" for L in LEVELS_M) + " \\\\",
        "\\midrule",
    ]
    for model, job_csv in jobs.items():
        rd, rv = read_reference(fig_csv, model)
        cd, cv = read_computed(job_csv)
        pub = [loglog_at(rd, rv, L) for L in LEVELS_M]
        com = [loglog_at(cd, cv, L) for L in LEVELS_M]
        rat = [(c / p if (c and p) else None) for c, p in zip(com, pub)]
        # the model name sits on the first of its three rows; no
        # \multirow, so the tables carry no package dependency beyond
        # booktabs
        lines.append(f"{model} & pub. & "
                     + " & ".join(fmt(x) for x in pub) + " \\\\")
        lines.append("& comp. & " + " & ".join(fmt(x) for x in com) + " \\\\")
        lines.append("& ratio & " + " & ".join(
            f"{x:.2f}" if x else "---" for x in rat) + " \\\\")
        lines.append("\\addlinespace")
    if lines[-1] == "\\addlinespace":
        lines.pop()
    lines += ["\\bottomrule", "\\end{tabular}", "\\end{table}", ""]
    return "\n".join(lines)


def main():
    outdir = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "report_tables"
    outdir.mkdir(parents=True, exist_ok=True)
    allp = outdir / "benchmark_rate_tables.tex"
    with allp.open("w") as f:
        f.write("% Generated by make_report_tables.py -- do not edit by hand.\n")
        f.write("% Regenerate after any change to computed/ or reference/.\n\n")
        for key, caption, fig_csv, jobs in EXERCISES:
            f.write(make_table(key, caption, fig_csv, jobs))
            f.write("\n")
    print(f"wrote {allp} ({len(EXERCISES)} tables)")


if __name__ == "__main__":
    sys.exit(main())
