# -*- coding: utf-8 -*-
"""
Generate Figure 1 of the SRL manuscript: primary (principal) surface rupture
probability models organized by faulting style.

Scientific reference
--------------------
Valentini, A., et al. (2025). Twenty-Five Years of Probabilistic Fault
Displacement Hazard Assessment. Reviews of Geophysics, 63(3).
https://doi.org/10.1029/2024RG000875  (Figure 10; Tables 3-4)

Data provenance
---------------
Every curve is produced by calling the model implementations registered in
``openquake.fdha.primary_surf_rup``. No equations are re-implemented and no
coefficients are hardcoded here.

Validity ranges
---------------
The framework models do NOT carry validity-range metadata and do NOT refuse or
flag out-of-range magnitudes -- each ``get_prob`` returns a logistic/numerical
value for any Mw. The published applicable Mw ranges below are taken from the
project user manual (docs/UserManual_Enhanced/06-Models.md), which tabulates
them from Valentini et al. (2025), Reviews of Geophysics, Table 4. Curves are
plotted ONLY within their published range, clipped to the figure axis
[5.0, 8.5]. No extrapolation is performed.

Known framework gaps (documented, not silently skipped)
-------------------------------------------------------
1. Panel A requests three Youngs et al. (2003) regional subsets (Cordillera,
   Great Basin, Northern Basin and Range). ``Youngs2003PrimarySR`` only exposes
   ``style='all'`` and ``style='normal'`` -- there are no regional subsets.
   Per user direction, the two available styles are plotted instead, and the
   three missing subsets are recorded here and in the caption.
2. ``Mammarella2024PrimarySR.get_prob`` requires MSR, HDD_str and auxiliary
   distribution parameters not given in the task. Per user direction, all three
   MSR variants (0=WC94, 1=Leonard 2010, 2=Thingbaijam 2017) are plotted; the
   auxiliary parameters reuse the project's existing validation-script values.
"""

from __future__ import annotations

import csv
import re
import sys
import traceback
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

# --- Framework model registry (the ONLY source of curve data) ----------------
from openquake.fdha.primary_surf_rup.youngs2003 import Youngs2003PrimarySR
from openquake.fdha.primary_surf_rup.wells_coppersmith1993 import WC1993PrimarySR
from openquake.fdha.primary_surf_rup.moss_ross2011 import MossRoss2011PrimarySR
from openquake.fdha.primary_surf_rup.moss2013 import Moss2013PrimarySR
from openquake.fdha.primary_surf_rup.takao2013 import Takao2013PrimarySR
from openquake.fdha.primary_surf_rup.yang2021 import Yang2021PrimarySR
from openquake.fdha.primary_surf_rup.pizza2023 import Pizza2023PrimarySR
from openquake.fdha.primary_surf_rup.mammarella2024 import Mammarella2024PrimarySR

OUT_DIR = Path(__file__).resolve().parent
AXIS_MIN, AXIS_MAX = 5.0, 8.5

# -----------------------------------------------------------------------------
# Published applicable Mw ranges (Valentini et al. 2025, Rev. Geophys., Table 4,
# as tabulated in docs/UserManual_Enhanced/06-Models.md).
# -----------------------------------------------------------------------------
MW_RANGE = {
    "youngs2003": (4.5, 7.6),
    "mammarella2024": (5.0, 8.0),
    "mossross2011": (5.5, 8.0),
    "moss2013": (4.2, 8.7),
    "pizza2023": (5.5, 7.9),
    "takao2013": (5.5, 7.4),
    "wc1993": (5.0, 8.2),
    "yang2021": (4.7, 6.6),
}

# Bibliographic metadata (APA). DOIs verified against the cited journals;
# Wells & Coppersmith (1993) is a Seismological Research Letters abstract with
# no registered DOI.
CITATION = {
    "youngs2003": (
        "Youngs, R. R., Arabasz, W. J., Anderson, R. E., et al. (2003). "
        "A methodology for probabilistic fault displacement hazard analysis "
        "(PFDHA). Earthquake Spectra, 19(1), 191-219.",
        "10.1193/1.1542891",
    ),
    "mammarella2024": (
        "Mammarella, L., Iervolino, I., & Valentini, A. (2024). Conditional "
        "probability of surface rupture: A numerical approach for principal "
        "faulting. Earthquake Spectra.",
        "10.1177/87552930241293570",
    ),
    "mossross2011": (
        "Moss, R. E. S., & Ross, Z. E. (2011). Probabilistic fault "
        "displacement hazard analysis for reverse faults. Bulletin of the "
        "Seismological Society of America, 101(4), 1542-1553.",
        "10.1785/0120100248",
    ),
    "moss2013": (
        "Moss, R. E. S., Stanton, K. V., & Buelna, M. I. (2013). The impact "
        "of material stiffness on the likelihood of fault rupture propagating "
        "to the ground surface. Seismological Research Letters, 84(3), "
        "485-488.",
        "10.1785/0220120113",
    ),
    "pizza2023": (
        "Pizza, M., Ferrario, M. F., Thomas, F., Tringali, G., & Livio, F. "
        "(2023). Likelihood of primary surface faulting: Updating of "
        "empirical regressions. Bulletin of the Seismological Society of "
        "America, 113(5), 2106-2118.",
        "10.1785/0120230019",
    ),
    "takao2013": (
        "Takao, M., Tsuchiyama, J., Annaka, T., & Kurita, T. (2013). "
        "Application of probabilistic fault displacement hazard analysis in "
        "Japan. Journal of Japan Association for Earthquake Engineering, "
        "13(1), 17-36.",
        "10.5610/jaee.13.1_17",
    ),
    "wc1993": (
        "Wells, D. L., & Coppersmith, K. J. (1993). Likelihood of surface "
        "rupture as a function of magnitude (abstract). Seismological "
        "Research Letters, 64(1), 54.",
        "",  # SRL abstract; no registered DOI
    ),
    "yang2021": (
        "Yang, H., Quigley, M., & King, T. (2021). Surface slip distributions "
        "and geometric complexity of intraplate reverse-faulting earthquakes. "
        "GSA Bulletin, 133(9-10), 1909-1929.",
        "10.1130/B35809.1",
    ),
}

# -----------------------------------------------------------------------------
# Legend naming convention (SRL).
#   >= 3 authors : "FirstAuthor et al. (Year)"
#   2 authors    : "FirstAuthor and SecondAuthor (Year)"
#   1 author     : "Author (Year)"
# No ampersand; "et al." is plain (not italic); model variant appended ", ...".
#
# The framework model classes (Youngs2003PrimarySR, etc.) carry no
# ``display_name`` attribute, so there is no registry string to validate.
# Instead, the citation prefix is *constructed* from the authoritative author
# list below, and every assembled legend label is validated against the
# convention before plotting (see ``validate_label``).
# -----------------------------------------------------------------------------
CITE_AUTHORS = {
    # key -> (first author surname, total number of authors, year)
    "youngs2003": ("Youngs", 20, 2003),
    "mammarella2024": ("Mammarella", 3, 2024),
    "mossross2011": ("Moss", 2, 2011),       # Moss and Ross
    "moss2013": ("Moss", 3, 2013),           # Moss, Stanton, Buelna
    "pizza2023": ("Pizza", 5, 2023),
    "takao2013": ("Takao", 4, 2013),
    "wc1993": ("Wells", 2, 1993),            # Wells and Coppersmith
    "yang2021": ("Yang", 3, 2021),
}
# Second-author surname, required only for two-author papers.
SECOND_AUTHOR = {"mossross2011": "Ross", "wc1993": "Coppersmith"}

# Convention regex: "<cite> (YYYY)" with an optional ", <variant>" suffix.
_LABEL_RE = re.compile(
    r"^[A-Z][\w.'‐-]*"
    r"(?: et al\.| and [A-Z][\w.'‐-]*)?"
    r" \(\d{4}\)(?:, .+)?$"
)


def cite_prefix(key: str) -> str:
    """Build the author-year citation per the SRL naming convention."""
    first, n_authors, year = CITE_AUTHORS[key]
    if n_authors >= 3:
        return f"{first} et al. ({year})"
    if n_authors == 2:
        return f"{first} and {SECOND_AUTHOR[key]} ({year})"
    return f"{first} ({year})"


def make_label(key: str, variant: str = "") -> str:
    """Assemble a legend label: citation prefix plus optional ', variant'."""
    prefix = cite_prefix(key)
    return f"{prefix}, {variant}" if variant else prefix


def validate_label(model_id: str, label: str) -> None:
    """STOP if a legend label violates the SRL naming convention."""
    if "&" in label:
        stop(model_id, ValueError(
            f"legend label uses an ampersand (spell out 'and'): {label!r}"))
    if "$" in label or r"\it" in label:
        stop(model_id, ValueError(
            f"legend label contains italic markup ('et al.' must be "
            f"plain): {label!r}"))
    if not _LABEL_RE.match(label):
        stop(model_id, ValueError(
            f"malformed legend label, does not match the naming "
            f"convention: {label!r}"))


# Mammarella auxiliary parameters: reused verbatim from the project's existing
# validation script openquake/fdha/test/scripts/generate_mammarella_style_comparison.py
MAMM_AUX = dict(dip_sigma=2.0, t_d=2.5, Zs_sigma=2.0, t_z=1.0)
MAMM_MSR = {0: "WC94", 1: "Leonard2010", 2: "Thingbaijam2017"}
MAMM_HDD = {"normal": "AGG_N", "reverse": "AGG_R", "strike-slip": "AGG_S"}
MAMM_DIP = {"normal": 50.0, "reverse": 30.0, "strike-slip": 90.0}
SEISMO_THICKNESS = 9.0  # km (task: Tseis = 9 km)


def clipped_range(key: str) -> tuple[float, float]:
    """Published Mw range intersected with the figure axis [5.0, 8.5]."""
    lo, hi = MW_RANGE[key]
    return max(lo, AXIS_MIN), min(hi, AXIS_MAX)


def fine_grid(lo: float, hi: float, step: float) -> np.ndarray:
    """Magnitude grid on [lo, hi] inclusive at the requested step."""
    n = int(round((hi - lo) / step))
    return lo + step * np.arange(n + 1)


def stop(model_id: str, exc: Exception) -> None:
    """Honor the STOP conditions: report the failing model and abort."""
    print("\n" + "=" * 70, file=sys.stderr)
    print(f"STOP: model '{model_id}' failed on the requested magnitude grid.",
          file=sys.stderr)
    print(f"Error: {type(exc).__name__}: {exc}", file=sys.stderr)
    traceback.print_exc()
    print("=" * 70, file=sys.stderr)
    sys.exit(1)


# -----------------------------------------------------------------------------
# Curve evaluators -- each returns probabilities by calling the framework model.
# -----------------------------------------------------------------------------
def eval_logistic(model_call, key: str, step: float = 0.05):
    """Evaluate a fast logistic model over its clipped, published range."""
    lo, hi = clipped_range(key)
    m = fine_grid(lo, hi, step)
    p = np.asarray(model_call(m), dtype=float)
    return m, p


def eval_mammarella(style: str, msr: int, step: float = 0.25):
    """Evaluate Mammarella et al. (2024) for one style/MSR over its range."""
    lo, hi = clipped_range("mammarella2024")
    m = fine_grid(lo, hi, step)
    model = Mammarella2024PrimarySR()
    p = model.get_prob(
        m,
        MSR=msr,
        HDD_str=MAMM_HDD[style],
        dip_mu=MAMM_DIP[style],
        style=style,
        seismothickness=SEISMO_THICKNESS,
        **MAMM_AUX,
    )
    return m, np.atleast_1d(np.asarray(p, dtype=float))


# -----------------------------------------------------------------------------
# Panel / curve specification.
# Each curve: (model_id, legend label, style key, color, linestyle, marker,
#              evaluator callable -> (m, p))
# -----------------------------------------------------------------------------
def build_curves():
    y = Youngs2003PrimarySR()
    wc = WC1993PrimarySR()
    mr = MossRoss2011PrimarySR()
    m13 = Moss2013PrimarySR()
    tk = Takao2013PrimarySR()
    yg = Yang2021PrimarySR()
    pz = Pizza2023PrimarySR()

    panels = {
        "Normal Faulting": [],
        "Reverse Faulting": [],
        "Strike-Slip Faulting": [],
        "All Slip": [],
    }
    meta = []  # rows for the metadata CSV

    # ---------------- Panel A: Normal ----------------
    panels["Normal Faulting"] = [
        dict(model_id="youngs2003_normal",
             label=make_label("youngs2003", "normal subset"),
             key="youngs2003", color="#1f77b4", ls="-", marker="o",
             eval=lambda: eval_logistic(lambda m: y.get_prob(m, style="normal"),
                                        "youngs2003"),
             cite="youngs2003", fstyle="normal",
             cls="openquake.fdha.primary_surf_rup.youngs2003.Youngs2003PrimarySR"),
        dict(model_id="youngs2003_all",
             label=make_label("youngs2003", "regional 'all' subset"),
             key="youngs2003", color="#17becf", ls="--", marker="v",
             eval=lambda: eval_logistic(lambda m: y.get_prob(m, style="all"),
                                        "youngs2003"),
             cite="youngs2003", fstyle="normal",
             cls="openquake.fdha.primary_surf_rup.youngs2003.Youngs2003PrimarySR"),
        dict(model_id="pizza2023_normal",
             label=make_label("pizza2023", "normal subset"),
             key="pizza2023", color="#2ca02c", ls="-", marker="s",
             eval=lambda: eval_logistic(lambda m: pz.get_prob(m, style="normal"),
                                        "pizza2023"),
             cite="pizza2023", fstyle="normal",
             cls="openquake.fdha.primary_surf_rup.pizza2023.Pizza2023PrimarySR"),
    ]
    for msr, mlbl in MAMM_MSR.items():
        panels["Normal Faulting"].append(
            dict(model_id=f"mammarella2024_normal_msr{msr}",
                 label=make_label("mammarella2024", f"MSR={mlbl}"),
                 key="mammarella2024",
                 color=["#d62728", "#9467bd", "#8c564b"][msr], ls=":",
                 marker=["^", "D", "P"][msr],
                 eval=(lambda s="normal", mm=msr: eval_mammarella(s, mm)),
                 cite="mammarella2024", fstyle="normal",
                 cls="openquake.fdha.primary_surf_rup.mammarella2024.Mammarella2024PrimarySR"))

    # ---------------- Panel B: Reverse ----------------
    panels["Reverse Faulting"] = [
        dict(model_id="mossross2011", label=make_label("mossross2011"),
             key="mossross2011", color="#1f77b4", ls="-", marker="o",
             eval=lambda: eval_logistic(mr.get_prob, "mossross2011"),
             cite="mossross2011", fstyle="reverse",
             cls="openquake.fdha.primary_surf_rup.moss_ross2011.MossRoss2011PrimarySR"),
        dict(model_id="moss2013_reverse_stiff",
             label=make_label("moss2013", "reverse, stiff soil"),
             key="moss2013", color="#2ca02c", ls="-", marker="s",
             eval=lambda: eval_logistic(
                 lambda m: m13.get_prob(m, "reverse", 800.0), "moss2013"),
             cite="moss2013", fstyle="reverse",
             cls="openquake.fdha.primary_surf_rup.moss2013.Moss2013PrimarySR"),
        dict(model_id="moss2013_reverse_soft",
             label=make_label("moss2013", "reverse, soft soil"),
             key="moss2013", color="#98df8a", ls="--", marker="D",
             eval=lambda: eval_logistic(
                 lambda m: m13.get_prob(m, "reverse", 300.0), "moss2013"),
             cite="moss2013", fstyle="reverse",
             cls="openquake.fdha.primary_surf_rup.moss2013.Moss2013PrimarySR"),
        dict(model_id="takao2013_reverse",
             label=make_label("takao2013", "reverse"),
             key="takao2013", color="#ff7f0e", ls="-", marker="^",
             eval=lambda: eval_logistic(tk.get_prob, "takao2013"),
             cite="takao2013", fstyle="reverse",
             cls="openquake.fdha.primary_surf_rup.takao2013.Takao2013PrimarySR"),
        dict(model_id="yang2021", label=make_label("yang2021"),
             key="yang2021", color="#e377c2", ls="-", marker="X",
             eval=lambda: eval_logistic(yg.get_prob, "yang2021"),
             cite="yang2021", fstyle="reverse",
             cls="openquake.fdha.primary_surf_rup.yang2021.Yang2021PrimarySR"),
        dict(model_id="pizza2023_reverse",
             label=make_label("pizza2023", "reverse subset"),
             key="pizza2023", color="#7f7f7f", ls="-", marker="*",
             eval=lambda: eval_logistic(lambda m: pz.get_prob(m, style="reverse"),
                                        "pizza2023"),
             cite="pizza2023", fstyle="reverse",
             cls="openquake.fdha.primary_surf_rup.pizza2023.Pizza2023PrimarySR"),
    ]
    for msr, mlbl in MAMM_MSR.items():
        panels["Reverse Faulting"].append(
            dict(model_id=f"mammarella2024_reverse_msr{msr}",
                 label=make_label("mammarella2024", f"MSR={mlbl}"),
                 key="mammarella2024",
                 color=["#d62728", "#9467bd", "#8c564b"][msr], ls=":",
                 marker=["p", "h", "<"][msr],
                 eval=(lambda s="reverse", mm=msr: eval_mammarella(s, mm)),
                 cite="mammarella2024", fstyle="reverse",
                 cls="openquake.fdha.primary_surf_rup.mammarella2024.Mammarella2024PrimarySR"))

    # ---------------- Panel C: Strike-Slip ----------------
    panels["Strike-Slip Faulting"] = [
        dict(model_id="moss2013_ss_stiff",
             label=make_label("moss2013", "strike-slip, stiff soil"),
             key="moss2013", color="#2ca02c", ls="-", marker="s",
             eval=lambda: eval_logistic(
                 lambda m: m13.get_prob(m, "strike-slip", 800.0), "moss2013"),
             cite="moss2013", fstyle="strike-slip",
             cls="openquake.fdha.primary_surf_rup.moss2013.Moss2013PrimarySR"),
        dict(model_id="moss2013_ss_soft",
             label=make_label("moss2013", "strike-slip, soft soil"),
             key="moss2013", color="#98df8a", ls="--", marker="D",
             eval=lambda: eval_logistic(
                 lambda m: m13.get_prob(m, "strike-slip", 300.0), "moss2013"),
             cite="moss2013", fstyle="strike-slip",
             cls="openquake.fdha.primary_surf_rup.moss2013.Moss2013PrimarySR"),
        dict(model_id="takao2013_strikeslip",
             label=make_label("takao2013", "strike-slip"),
             key="takao2013", color="#ff7f0e", ls="-", marker="^",
             eval=lambda: eval_logistic(tk.get_prob, "takao2013"),
             cite="takao2013", fstyle="strike-slip",
             cls="openquake.fdha.primary_surf_rup.takao2013.Takao2013PrimarySR"),
        dict(model_id="pizza2023_strikeslip",
             label=make_label("pizza2023", "strike-slip subset"),
             key="pizza2023", color="#7f7f7f", ls="-", marker="*",
             eval=lambda: eval_logistic(
                 lambda m: pz.get_prob(m, style="strike-slip"), "pizza2023"),
             cite="pizza2023", fstyle="strike-slip",
             cls="openquake.fdha.primary_surf_rup.pizza2023.Pizza2023PrimarySR"),
    ]
    for msr, mlbl in MAMM_MSR.items():
        panels["Strike-Slip Faulting"].append(
            dict(model_id=f"mammarella2024_strikeslip_msr{msr}",
                 label=make_label("mammarella2024", f"MSR={mlbl}"),
                 key="mammarella2024",
                 color=["#d62728", "#9467bd", "#8c564b"][msr], ls=":",
                 marker=["p", "h", "<"][msr],
                 eval=(lambda s="strike-slip", mm=msr: eval_mammarella(s, mm)),
                 cite="mammarella2024", fstyle="strike-slip",
                 cls="openquake.fdha.primary_surf_rup.mammarella2024.Mammarella2024PrimarySR"))

    # ---------------- Panel D: All Slip ----------------
    panels["All Slip"] = [
        dict(model_id="wc1993", label=make_label("wc1993"),
             key="wc1993", color="#1f77b4", ls="-", marker="o",
             eval=lambda: eval_logistic(wc.get_prob, "wc1993"),
             cite="wc1993", fstyle="all",
             cls="openquake.fdha.primary_surf_rup.wells_coppersmith1993.WC1993PrimarySR"),
        dict(model_id="pizza2023_all",
             label=make_label("pizza2023", "all-style"),
             key="pizza2023", color="#2ca02c", ls="-", marker="s",
             eval=lambda: eval_logistic(lambda m: pz.get_prob(m, style="all"),
                                        "pizza2023"),
             cite="pizza2023", fstyle="all",
             cls="openquake.fdha.primary_surf_rup.pizza2023.Pizza2023PrimarySR"),
    ]
    return panels, meta


# -----------------------------------------------------------------------------
def main():
    rcParams["font.family"] = "serif"
    rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    # Render mathtext (the italic m) in Times New Roman as well.
    rcParams["mathtext.fontset"] = "custom"
    rcParams["mathtext.rm"] = "Times New Roman"
    rcParams["mathtext.it"] = "Times New Roman:italic"
    rcParams["mathtext.bf"] = "Times New Roman:bold"
    rcParams["font.size"] = 12
    rcParams["axes.labelsize"] = 12
    rcParams["axes.titlesize"] = 12
    rcParams["legend.fontsize"] = 12
    rcParams["xtick.labelsize"] = 12
    rcParams["ytick.labelsize"] = 12
    # SSA Art Guidelines: embed fonts as TrueType (Type 42) outlines in PDF/EPS
    # so SSA typesetters receive fully embedded fonts.
    rcParams["pdf.fonttype"] = 42
    rcParams["ps.fonttype"] = 42
    # SSA: rule/axis line weights 0.5-1.0 pt at 100% reproduction scale.
    rcParams["axes.linewidth"] = 0.8

    panels, _ = build_curves()
    panel_titles = list(panels.keys())
    # SRL convention: lowercase parenthesized panel labels (a) (b) (c) (d).
    panel_letters = ["a", "b", "c", "d"]

    fig, axes2d = plt.subplots(2, 2, figsize=(16, 14))
    axes = axes2d.ravel()

    seen_keys: dict[str, str] = {}  # registry collision guard
    meta_rows = []

    for ax, letter, title in zip(axes, panel_letters, panel_titles):
        for spec in panels[title]:
            mid = spec["model_id"]
            if mid in seen_keys:
                print(f"STOP: duplicate registry key '{mid}' "
                      f"(already used by {seen_keys[mid]}).", file=sys.stderr)
                sys.exit(1)
            seen_keys[mid] = title

            # Enforce the SRL legend naming convention (STOP on violation).
            validate_label(mid, spec["label"])

            try:
                m, p = spec["eval"]()
            except Exception as exc:  # noqa: BLE001
                stop(mid, exc)

            if m is None or len(m) == 0 or np.asarray(p).size != len(m):
                stop(mid, ValueError(
                    f"model returned no usable values "
                    f"(m={None if m is None else len(m)}, "
                    f"p={np.asarray(p).size})"))

            p = np.asarray(p, dtype=float)
            if not np.all(np.isfinite(p)):
                stop(mid, ValueError("model returned non-finite probabilities"))

            # marker spacing: every 0.5 Mw. Fine logistic grid uses 0.05 step
            # (markevery=10); Mammarella uses 0.25 step (markevery=2).
            step = round(float(m[1] - m[0]), 4) if len(m) > 1 else 0.05
            markevery = max(1, int(round(0.5 / step)))

            ax.plot(m, p, color=spec["color"], linestyle=spec["ls"],
                    marker=spec["marker"], markevery=markevery, markersize=7,
                    linewidth=2.4, label=spec["label"])

            lo, hi = MW_RANGE[spec["key"]]
            cit, doi = CITATION[spec["cite"]]
            meta_rows.append(dict(
                model_id=mid, faulting_style=spec["fstyle"],
                mw_min=lo, mw_max=hi, source_citation=cit,
                source_doi=doi, framework_class_path=spec["cls"]))

        ax.set_xlim(AXIS_MIN, AXIS_MAX)
        ax.set_ylim(0.0, 1.0)
        ax.set_xlabel(r"Magnitude, $m$")
        ax.set_ylabel("Probability of Primary Surface Rupture")
        ax.set_title(title, fontweight="bold")
        # SRL convention: lowercase parenthesized label, top-left OUTSIDE panel.
        ax.text(-0.10, 1.04, f"({letter})", transform=ax.transAxes,
                fontsize=14, fontweight="bold", ha="left", va="bottom")
        ax.grid(True, which="major", color="0.85", linewidth=0.6)
        ax.grid(False, which="minor")
        # Legend placed below the panel so it never overlaps the curves.
        n_curves = len(panels[title])
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13),
                  ncol=2 if n_curves > 3 else 1, framealpha=0.92,
                  borderaxespad=0.0)

    fig.tight_layout()

    pdf = OUT_DIR / "figure_01_primary_surface_rupture_probability.pdf"
    png = OUT_DIR / "figure_01_primary_surface_rupture_probability.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {pdf}")
    print(f"wrote {png}")

    csv_path = OUT_DIR / "figure_01_metadata.csv"
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            "model_id", "faulting_style", "mw_min", "mw_max",
            "source_citation", "source_doi", "framework_class_path"])
        writer.writeheader()
        for row in meta_rows:
            writer.writerow(row)
    print(f"wrote {csv_path}  ({len(meta_rows)} curves)")


if __name__ == "__main__":
    main()
