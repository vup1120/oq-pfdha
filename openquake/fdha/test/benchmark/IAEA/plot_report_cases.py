#!/usr/bin/env python
"""Per-case-study overlay figures for the validation section of the report.

One figure per IAEA case study (Kumamoto, Le Teil, Norcia), one panel per
exercise of TECDOC-2092 Table 8: the principal exercise (two options for
Kumamoto) and the three distributed exercises. Solid lines are oq-pfdha
curves (``computed/`` CSVs written by ``run_all.py`` plus the
supplementary jobs listed below); dashed lines are the published team
curves (``reference/`` CSVs). Every panel carries both.

All twelve exercises now carry reference curves: the coordinator-supplied
vectors for the principal and 1st-distributed exercises, and the curves
digitized from the published figures by
``reference/digitize_tecdoc_figures.py`` for the 2nd and 3rd distributed
exercises. The only supplementary job without a published counterpart is
``le_teil/job_principal_M24.ini`` (the numerical P_sr branch), which no
exercise team ran.

Usage (repo root)::

    PYTHONPATH=. python openquake/fdha/test/benchmark/IAEA/plot_report_cases.py [outdir]
"""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
COMPUTED = HERE / "computed"
REFERENCE = HERE / "reference"

# Models excluded from the exercise comparison. The teams' V24 curves were
# produced with an earlier, under-review revision of Visini et al. (2025):
# plotting them against an implementation of the published model would show a
# vintage mismatch, not a validation result. That implementation is validated
# at model level against the authors' released FDHLab code instead
# (benchmark/visini_et_al_2025); see manifest.py for the full rationale.
EXCLUDED_MODELS = {"V24"}

MODEL_COLORS = {
    "P11": "tab:blue", "T13": "tab:orange", "K24": "tab:green",
    "L23": "tab:red", "C24": "tab:purple", "M11": "tab:brown",
    "Y03": "tab:blue", "V24": "tab:red", "M24": "black",
}


def _computed(job_csv: str):
    rows = list(csv.DictReader((COMPUTED / job_csv).open()))
    d = np.array([float(r["disp_m"]) for r in rows])
    col = "rate_post_factor" if "rate_post_factor" in rows[0] else "rate"
    v = np.array([float(r[col]) for r in rows])
    return d, v


def _supplementary(outdir_name: str):
    """Aggregate curve of one of the supplementary jobs (run into a scratch
    directory or next to the job INI); returns None when not available."""
    for base in (HERE / "computed_supplementary", Path(os.environ.get(
            "IAEA_SUPP_DIR", ""))):
        f = base / outdir_name / "aggregate_hazard.csv" if base else None
        if f and f.exists():
            rows = list(csv.DictReader(f.open()))
            d = np.array([float(r["D0"]) for r in rows])
            v = np.array([float(r["mean"]) for r in rows])
            return d, v
    return None


def _diagnostic(rel_path: str):
    """A curve that is deliberately not a benchmark entry (see the panel
    spec's ``diag`` slot); returns None when it has not been generated."""
    f = HERE / "computed_supplementary" / rel_path
    if not f.exists():
        return None
    rows = list(csv.DictReader(f.open()))
    return (np.array([float(r["disp_m"]) for r in rows]),
            np.array([float(r["rate"]) for r in rows]))


def _reference(fig_csv: str):
    rows = list(csv.DictReader((REFERENCE / fig_csv).open()))
    d = np.array([float(r["disp_cm"]) for r in rows]) / 100.0
    out = {}
    for c in rows[0]:
        if c == "disp_cm":
            continue
        keep = [i for i, r in enumerate(rows) if r[c].strip()]
        if keep:
            out[c] = (d[keep], np.array([float(rows[i][c]) for i in keep]))
    return out


PANELS = {
    "kumamoto": [
        ("Principal, option 1 (fixed rupture)", "fig4a_kumamoto_principal.csv",
         {"P11": "kumamoto_principal_P11.csv", "C24": "kumamoto_principal_C24.csv",
          "K24": "kumamoto_principal_K24.csv", "T13": "kumamoto_principal_T13.csv",
          "L23": "kumamoto_principal_L23.csv"}, {}),
        ("Principal, option 2 (floating ruptures)",
         "fig4b_kumamoto_principal_floating.csv",
         {"K24": "kumamoto_floating_K24.csv", "T13": "kumamoto_floating_T13.csv",
          "L23": "kumamoto_floating_L23.csv"}, {}),
        ("1st distributed (r = 5.2 km)", "fig6a_kumamoto_distributed.csv",
         {"P11": "kumamoto_distributed_P11.csv",
          "T13": "kumamoto_distributed_T13.csv"}, {}),
        ("2nd distributed (r = 10 km)", "fig15c_kumamoto_distributed_r10.csv",
         {"P11": "kumamoto_distributed_P11_sens4.csv",
          "T13": "kumamoto_distributed_T13_sens4.csv"}, {}),
        ("3rd distributed (r = 0.6 km, Suizenji)",
         "fig15d_kumamoto_distributed_suizenji.csv",
         {"T13": "kumamoto_distributed_T13_sens1.csv"}, {},
         {"P11": ("kumamoto_distributed_P11_sens1_extrapolated/curve.csv",
                  "P11 extrapolated below M6\n(diagnostic, 7.6%)")}),
    ],
    "le_teil": [
        ("Principal", "fig4c_leteil_principal.csv",
         {"K24": "le_teil_principal_K24.csv", "T13": "le_teil_principal_T13.csv",
          "L23": "le_teil_principal_L23.csv", "M11": "le_teil_principal_M11.csv"},
         {"M24": "leteil_M24_out"}),
        ("1st distributed (r = 0.6 km, footwall)", "fig6b_leteil_distributed.csv",
         {"T13": "le_teil_distributed_T13.csv",}, {}),
        ("2nd distributed (dip 45$^\\circ$ NW, hanging wall)",
         "fig19c_leteil_distributed_dipNW.csv",
         {"T13": "le_teil_distributed_T13_dipflip.csv",}, {}),
        ("3rd distributed (three faults)",
         "fig19d_leteil_distributed_3faults.csv",
         {"T13": "le_teil_distributed_T13_sens4.csv",}, {}),
    ],
    "norcia": [
        ("Principal", "fig4d_norcia_principal.csv",
         {"Y03": "norcia_principal_Y03.csv", "K24": "norcia_principal_K24.csv",
          "L23": "norcia_principal_L23.csv"}, {}),
        ("1st distributed (MS, r = 7.6 km, MVF)", "fig6c_norcia_distributed.csv",
         {"Y03": "norcia_distributed_Y03.csv",}, {}),
        ("2nd distributed (SL, r = 2.4 km, MVF)",
         "fig22b_norcia_distributed_sl.csv",
         {"Y03": "norcia_distributed_Y03_sens2.csv",}, {}),
        ("3rd distributed (MS; MVF + NF)",
         "fig22c_norcia_distributed_2sources.csv",
         {"Y03": "norcia_distributed_Y03_sens3.csv",}, {}),
    ],
}


def make_figure(case: str, outdir: Path):
    panels = PANELS[case]
    n = len(panels)
    ncols = 2 if n <= 4 else 3
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.4 * nrows),
                             squeeze=False)
    for ax in axes.flat[n:]:
        ax.set_visible(False)
    letters = "abcdef"
    for i, spec in enumerate(panels):
        title, ref_csv, jobs, supp = spec[:4]
        diag = spec[4] if len(spec) > 4 else {}
        ax = axes.flat[i]
        modelled = set(jobs) | set(supp) | set(diag)
        if ref_csv:
            for model, (d, v) in _reference(ref_csv).items():
                if model in EXCLUDED_MODELS:
                    continue
                colour = MODEL_COLORS.get(model, "gray")
                if model in modelled:
                    # Drawn as a wide, pale dashed band *under* the framework
                    # curve: where the two agree closely (e.g. P11 at 1.9%) a
                    # thin dashed line would be hidden by the solid one, which
                    # would make the best agreements look like missing
                    # references. The halo stays visible either way.
                    ax.plot(d, v, "--", color=colour, lw=3.4, alpha=0.32,
                            zorder=1, solid_capstyle="butt")
                else:
                    # A published curve the framework does not reproduce at
                    # all (no computed and no diagnostic counterpart).
                    ax.plot(d, v, "--", color=colour, lw=1.5, alpha=0.9,
                            zorder=1, label=f"{model} published (not run)")
        for model, job_csv in jobs.items():
            d, v = _computed(job_csv)
            ax.plot(d, v, "-", color=MODEL_COLORS.get(model, "gray"),
                    lw=1.5, zorder=3, label=model)
        for model, outdir_name in supp.items():
            got = _supplementary(outdir_name)
            if got is None:
                continue
            d, v = got
            ax.plot(d, v, "-", color=MODEL_COLORS.get(model, "black"),
                    lw=1.5, zorder=3,
                    label=model + (" (dip NW)" if "dipflip" in
                                   outdir_name else ""))
        for model, (rel_path, label) in diag.items():
            got = _diagnostic(rel_path)
            if got is None:
                continue
            d, v = got
            # Dotted, never solid: this curve exists only to show that the
            # published reference is reproducible, and is not an asserted
            # benchmark entry (see kumamoto/diagnose_P11_suizenji.py).
            ax.plot(d, v, ":", color=MODEL_COLORS.get(model, "gray"),
                    lw=1.8, zorder=2, label=label)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1e-3, 10)
        ax.set_ylim(1e-12, 1e-3)
        ax.set_title(f"({letters[i]}) {title}", fontsize=8.5)
        ax.set_xlabel("Displacement (m)", fontsize=8)
        ax.set_ylabel("AFOE (1/yr)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.25, which="both")
        ax.legend(fontsize=6.5, loc="lower left", framealpha=0.6)
    fig.suptitle({"kumamoto": "Kumamoto (strike-slip)",
                  "le_teil": "Le Teil (reverse)",
                  "norcia": "Norcia (normal)"}[case], fontsize=11)
    from matplotlib.lines import Line2D
    conv = [Line2D([], [], color="0.35", lw=1.5, ls="-",
                   label="framework (this study)"),
            Line2D([], [], color="0.35", lw=3.4, ls="--", alpha=0.5,
                   label="published team curve")]
    fig.legend(handles=conv, loc="lower center", ncol=2, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, 0.005))
    fig.tight_layout(rect=(0, 0.045, 1, 0.96))
    out = outdir / f"iaea_{case}_cases.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"), dpi=180)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    outdir = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "Figures"
    outdir.mkdir(parents=True, exist_ok=True)
    for case in PANELS:
        make_figure(case, outdir)
