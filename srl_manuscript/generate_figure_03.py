# -*- coding: utf-8 -*-
"""
Generate Figure 3 of the SRL manuscript: distributed (secondary) surface
rupture probability models versus distance from the principal fault.

Style reference: Valentini et al. (2025), Reviews of Geophysics, Figure 12.
Three subpanels in a row: (a) normal, (b) reverse, (c) strike-slip. For
dip-slip styles the curve is drawn continuously across r = 0 (negative r =
footwall, positive r = hanging wall); strike-slip is shown for r >= 0 only.
Reference magnitude m = 6.5 for magnitude-dependent models.

All curves are produced by calling the model implementations registered in
``openquake.fdha.secondary_surf_rup``. No equations are re-implemented and
no coefficients are hardcoded here. Validity ranges are taken from
docs/UserManual_Enhanced/06-Models.md (sourced from Valentini et al. (2025),
Reviews of Geophysics, Table 4).

User-directed configuration (decisions from the planning Q&A):
- Visini et al. (2025): plot all three combinations A, B, C per pixel size.
- Rodriguez Padilla and Oskin (2023): include in panel (c).
- Hanging-wall / footwall: continuous curve across r = 0 for dip-slip.
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

# --- Framework model registry --------------------------------------------
from openquake.fdha.secondary_surf_rup.youngs2003 import Youngs2003SecondarySR
from openquake.fdha.secondary_surf_rup.ferrario2021 import (
    FerrarioLivio2021SecondarySR,
)
from openquake.fdha.secondary_surf_rup.visini2025 import Visini2025SecondarySR
from openquake.fdha.secondary_surf_rup.takao2013 import Takao2013SecondarySR
from openquake.fdha.secondary_surf_rup.takao2014 import Takao2014SecondarySR
from openquake.fdha.secondary_surf_rup.petersen2011 import (
    Petersen2011SecondarySR,
)
from openquake.fdha.secondary_surf_rup.rodriguez2023 import (
    Rodriguez2023SecondarySR,
)


OUT_DIR = Path(__file__).resolve().parent
REF_MAG = 6.5

# -----------------------------------------------------------------------------
# Plot styling -- matches Figures 1 and 2.
# -----------------------------------------------------------------------------
rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "font.size": 12,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "axes.linewidth": 0.8,
    "grid.linewidth": 0.6,
    "lines.linewidth": 2.4,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "mathtext.fontset": "custom",
    "mathtext.rm": "Times New Roman",
    "mathtext.it": "Times New Roman:italic",
    "mathtext.bf": "Times New Roman:bold",
    "mathtext.sf": "Times New Roman",
    "mathtext.tt": "Times New Roman",
    "mathtext.cal": "Times New Roman:italic",
    "mathtext.fallback": "cm",
})

# -----------------------------------------------------------------------------
# Validity ranges (Valentini et al. 2025, Rev. Geophys., Table 4, as
# tabulated in docs/UserManual_Enhanced/06-Models.md).
# r ranges are in km. Asymmetric HW/FW limits where the source defines them.
# -----------------------------------------------------------------------------
MW_RANGE = {
    "youngs2003_secondary": (5.5, 7.4),
    "ferrariolivio2021": (6.0, 7.5),
    "visini2025_normal": (5.5, 7.9),
    "visini2025_reverse": (4.9, 7.9),
    "takao2013_secondary": (5.8, 7.4),
    "takao2014": (5.8, 7.4),   # carried over from Takao 2013 calibration
    "petersen2011_secondary": (6.5, 7.5),
    "rodriguez2023": (None, None),  # mag-independent
}
# Per-side r limits in km. (fw_max_km, hw_max_km); negative side uses fw.
R_RANGE_KM = {
    "youngs2003_secondary": (15.0, 15.0),
    "ferrariolivio2021": (12.5, 15.5),
    "visini2025": (8.0, 10.0),
    "takao2013_secondary": (25.0, 25.0),
    "takao2014": (25.0, 25.0),
    "petersen2011_secondary": (2.5, 2.5),
    "rodriguez2023": (3.0, 3.0),
}

# -----------------------------------------------------------------------------
# Legend naming -- SRL convention (mirrors generate_figure_02.py).
# -----------------------------------------------------------------------------
CITE_AUTHORS = {
    "youngs2003": ("Youngs", 20, 2003),
    "ferrariolivio2021": ("Ferrario", 2, 2021),   # Ferrario and Livio
    "visini2025": ("Visini", 7, 2025),
    "takao2013": ("Takao", 4, 2013),
    "takao2014": ("Takao", 5, 2014),
    "petersen2011": ("Petersen", 8, 2011),
    "rodriguez2023": ("Rodriguez Padilla", 2, 2023),
}
SECOND_AUTHOR = {
    "ferrariolivio2021": "Livio",
    "rodriguez2023": "Oskin",
}

CITATION = {
    "youngs2003": (
        "Youngs, R. R., Arabasz, W. J., Anderson, R. E., et al. (2003). "
        "A methodology for probabilistic fault displacement hazard analysis "
        "(PFDHA). Earthquake Spectra, 19(1), 191-219.",
        "10.1193/1.1542891",
    ),
    "ferrariolivio2021": (
        "Ferrario, M. F., & Livio, F. (2021). Conditional probability of "
        "distributed surface rupturing during normal-faulting earthquakes. "
        "Solid Earth, 12(5), 1197-1209.",
        "10.5194/se-12-1197-2021",
    ),
    "visini2025": (
        "Visini, F., Boncio, P., Valentini, A., Scotti, O., Nurminen, F., "
        "Baize, S., & Pace, B. (2025). Empirical regressions for distributed "
        "faulting of dip-slip earthquakes. Earthquake Spectra.",
        "10.1177/87552930241308860",
    ),
    "takao2013": (
        "Takao, M., Tsuchiyama, J., Annaka, T., & Kurita, T. (2013). "
        "Application of probabilistic fault displacement hazard analysis in "
        "Japan. Journal of Japan Association for Earthquake Engineering, "
        "13(1), 17-36.",
        "10.5610/jaee.13.1_17",
    ),
    "takao2014": (
        "Takao, M., Ueta, K., Annaka, T., Kurita, T., Nakase, H., Kyoya, T., "
        "& Kato, J. (2014). Reliability improvement of probabilistic fault "
        "displacement hazard analysis. Journal of Japan Association for "
        "Earthquake Engineering, 14(2), 2-16-2-36.",
        "10.5610/jaee.14.2_16",
    ),
    "petersen2011": (
        "Petersen, M. D., Dawson, T. E., Chen, R., Cao, T., Wills, C. J., "
        "Schwartz, D. P., & Frankel, A. D. (2011). Fault displacement hazard "
        "for strike-slip faults. Bulletin of the Seismological Society of "
        "America, 101(2), 805-825.",
        "10.1785/0120100035",
    ),
    "rodriguez2023": (
        "Rodriguez Padilla, A. M., & Oskin, M. E. (2023). Displacement hazard "
        "from distributed ruptures in strike-slip earthquakes. Bulletin of "
        "the Seismological Society of America, 113(6), 2730-2745.",
        "10.1785/0120230044",
    ),
}

_LABEL_RE = re.compile(
    r"^[A-Z][\w.' ‐-]*"
    r"(?: et al\.| and [A-Z][\w.' ‐-]*)?"
    r" \(\d{4}\)(?:, .+)?$"
)


def stop(model_id: str, exc: Exception) -> None:
    print("\n" + "=" * 70, file=sys.stderr)
    print(f"STOP: model '{model_id}' failed.", file=sys.stderr)
    print(f"Error: {type(exc).__name__}: {exc}", file=sys.stderr)
    traceback.print_exc()
    print("=" * 70, file=sys.stderr)
    sys.exit(1)


def cite_prefix(key: str) -> str:
    first, n_authors, year = CITE_AUTHORS[key]
    if n_authors >= 3:
        return f"{first} et al. ({year})"
    if n_authors == 2:
        return f"{first} and {SECOND_AUTHOR[key]} ({year})"
    return f"{first} ({year})"


def make_label(key: str, variant: str = "") -> str:
    prefix = cite_prefix(key)
    return f"{prefix}, {variant}" if variant else prefix


def validate_label(model_id: str, label: str) -> None:
    if "&" in label:
        stop(model_id, ValueError(f"label uses '&': {label!r}"))
    if "$" in label or r"\it" in label:
        stop(model_id, ValueError(f"label has italic markup: {label!r}"))
    if not _LABEL_RE.match(label):
        stop(model_id, ValueError(f"malformed legend label: {label!r}"))


# -----------------------------------------------------------------------------
# Distance grid utilities.
# -----------------------------------------------------------------------------
# Fine 0.05 km step; markers placed every 1 km (markevery = 20).
R_STEP_KM = 0.05
MARKER_EVERY = int(round(1.0 / R_STEP_KM))


def signed_grid(fw_max_km: float, hw_max_km: float,
                axis_limit_km: float = 10.0) -> np.ndarray:
    """Two-sided distance grid from -min(fw_max, axis_limit) to +min(hw_max,
    axis_limit), in km, with a tiny gap around r = 0 to avoid log(0) when a
    model treats r as the closest distance."""
    fw = min(fw_max_km, axis_limit_km)
    hw = min(hw_max_km, axis_limit_km)
    n_fw = int(round(fw / R_STEP_KM))
    n_hw = int(round(hw / R_STEP_KM))
    fw_pts = -np.linspace(fw, R_STEP_KM, n_fw) if n_fw > 0 else np.array([])
    hw_pts = np.linspace(R_STEP_KM, hw, n_hw) if n_hw > 0 else np.array([])
    return np.concatenate([fw_pts, hw_pts])


def one_sided_grid(r_max_km: float,
                   axis_limit_km: float = 10.0) -> np.ndarray:
    """Positive-only grid for strike-slip panels."""
    r_max = min(r_max_km, axis_limit_km)
    n = int(round(r_max / R_STEP_KM))
    return np.linspace(R_STEP_KM, r_max, n)


# -----------------------------------------------------------------------------
# Curve evaluators.
# -----------------------------------------------------------------------------
def _safe_arr(v):
    a = np.atleast_1d(np.asarray(v, dtype=float))
    return a


def eval_youngs2003(version: str, mag: float):
    """version: '1' (w/o random effects) or '3' (w/ random effects)."""
    m = Youngs2003SecondarySR()
    fw_max, hw_max = R_RANGE_KM["youngs2003_secondary"]
    r = signed_grid(fw_max, hw_max)
    rx = r                       # signed: rx > 0 = HW, rx < 0 = FW
    r_abs = np.abs(r)            # closest distance is always non-negative
    p = _safe_arr(m.get_prob(mag, rx, r_abs, version=version))
    return r, p


def eval_ferrario2021(version: str):
    m = FerrarioLivio2021SecondarySR()
    fw_max, hw_max = R_RANGE_KM["ferrariolivio2021"]
    r = signed_grid(fw_max, hw_max)
    p = _safe_arr(m.get_prob(np.abs(r), r, version=version))
    return r, p


def eval_visini(style: str, mag: float, pixel_size: int, combination: str):
    m = Visini2025SecondarySR()
    fw_max, hw_max = R_RANGE_KM["visini2025"]
    r = signed_grid(fw_max, hw_max)
    # Visini takes r in METERS (closest distance) and rx in meters (signed).
    r_m = np.abs(r) * 1000.0
    rx_m = r * 1000.0
    p = np.array([float(m.get_prob(mag, float(rr), float(rrx),
                                   style=style, pixel_size=pixel_size,
                                   combination=combination))
                  for rr, rrx in zip(r_m, rx_m)])
    return r, p


def eval_takao2013(mag: float, pixel_size: int = 500):
    m = Takao2013SecondarySR()
    fw_max, hw_max = R_RANGE_KM["takao2013_secondary"]
    r = signed_grid(fw_max, hw_max)
    p = _safe_arr(m.get_prob(np.abs(r), mag, pixel_size=pixel_size))
    return r, p


def eval_takao2014(pixel_size: int):
    m = Takao2014SecondarySR()
    fw_max, hw_max = R_RANGE_KM["takao2014"]
    r = signed_grid(fw_max, hw_max)
    p = _safe_arr(m.get_prob(np.abs(r), pixel_size=pixel_size))
    return r, p


def eval_takao2013_ss(mag: float, pixel_size: int = 500):
    m = Takao2013SecondarySR()
    r_max = R_RANGE_KM["takao2013_secondary"][1]
    r = one_sided_grid(r_max)
    p = _safe_arr(m.get_prob(r, mag, pixel_size=pixel_size))
    return r, p


def eval_takao2014_ss(pixel_size: int):
    m = Takao2014SecondarySR()
    r_max = R_RANGE_KM["takao2014"][1]
    r = one_sided_grid(r_max)
    p = _safe_arr(m.get_prob(r, pixel_size=pixel_size))
    return r, p


def eval_petersen2011_ss(cell_size: int):
    m = Petersen2011SecondarySR()
    r_max = R_RANGE_KM["petersen2011_secondary"][1]
    r = one_sided_grid(r_max)
    p = _safe_arr(m.get_prob(r, cell_size=cell_size))
    return r, p


def eval_rodriguez2023_ss():
    m = Rodriguez2023SecondarySR()
    r_max = R_RANGE_KM["rodriguez2023"][1]
    r = one_sided_grid(r_max)
    p = _safe_arr(m.get_prob(r, pixel_size=1))
    return r, p


# -----------------------------------------------------------------------------
# Panel specifications.
# -----------------------------------------------------------------------------
COMB_LS = {"A": "-", "B": "--", "C": ":"}
VISINI_COLOR = {100: "#d62728", 500: "#9467bd"}


def build_panels():
    panels = {"a": [], "b": [], "c": []}

    # ---- Panel (a) Normal ---------------------------------------------------
    panels["a"].extend([
        dict(model_id="youngs2003_normal_norand",
             label=make_label("youngs2003", "500 m, w/o random effects"),
             cite="youngs2003", side="both", pixel=500,
             variant="version=1",
             cls="openquake.fdha.secondary_surf_rup.youngs2003.Youngs2003SecondarySR",
             eval=lambda: eval_youngs2003("1", REF_MAG),
             color="#1f77b4", ls="-", marker="o",
             mw_key="youngs2003_secondary", r_key="youngs2003_secondary"),
        dict(model_id="youngs2003_normal_rand",
             label=make_label("youngs2003", "500 m, w/ random effects"),
             cite="youngs2003", side="both", pixel=500,
             variant="version=3",
             cls="openquake.fdha.secondary_surf_rup.youngs2003.Youngs2003SecondarySR",
             eval=lambda: eval_youngs2003("3", REF_MAG),
             color="#17becf", ls="--", marker="v",
             mw_key="youngs2003_secondary", r_key="youngs2003_secondary"),
        dict(model_id="ferrario2021_regular",
             label=make_label("ferrariolivio2021", "500 m, Regular"),
             cite="ferrariolivio2021", side="both", pixel=500,
             variant="version=regular",
             cls="openquake.fdha.secondary_surf_rup.ferrario2021.FerrarioLivio2021SecondarySR",
             eval=lambda: eval_ferrario2021("regular"),
             color="#2ca02c", ls="-", marker="s",
             mw_key="ferrariolivio2021", r_key="ferrariolivio2021"),
        dict(model_id="ferrario2021_conservative",
             label=make_label("ferrariolivio2021", "500 m, Conservative"),
             cite="ferrariolivio2021", side="both", pixel=500,
             variant="version=conservative",
             cls="openquake.fdha.secondary_surf_rup.ferrario2021.FerrarioLivio2021SecondarySR",
             eval=lambda: eval_ferrario2021("conservative"),
             color="#98df8a", ls="--", marker="D",
             mw_key="ferrariolivio2021", r_key="ferrariolivio2021"),
    ])
    for comb in ("A", "B", "C"):
        panels["a"].append(dict(
            model_id=f"visini2025_normal_500_{comb}",
            label=make_label("visini2025",
                             f"500 m, combination {comb}"),
            cite="visini2025", side="both", pixel=500,
            variant=f"combination={comb}",
            cls="openquake.fdha.secondary_surf_rup.visini2025.Visini2025SecondarySR",
            eval=(lambda c=comb: eval_visini("normal", REF_MAG, 500, c)),
            color=VISINI_COLOR[500], ls=COMB_LS[comb], marker="^",
            mw_key="visini2025_normal", r_key="visini2025"))

    # ---- Panel (b) Reverse --------------------------------------------------
    panels["b"].extend([
        dict(model_id="takao2013_reverse_500",
             label=make_label("takao2013", "500 m"),
             cite="takao2013", side="both", pixel=500,
             variant="pixel_size=500",
             cls="openquake.fdha.secondary_surf_rup.takao2013.Takao2013SecondarySR",
             eval=lambda: eval_takao2013(REF_MAG, pixel_size=500),
             color="#ff7f0e", ls="-", marker="o",
             mw_key="takao2013_secondary", r_key="takao2013_secondary"),
        dict(model_id="takao2014_reverse_500",
             label=make_label("takao2014", "500 m"),
             cite="takao2014", side="both", pixel=500,
             variant="pixel_size=500",
             cls="openquake.fdha.secondary_surf_rup.takao2014.Takao2014SecondarySR",
             eval=lambda: eval_takao2014(500),
             color="#bcbd22", ls="--", marker="s",
             mw_key="takao2014", r_key="takao2014"),
        dict(model_id="takao2014_reverse_100",
             label=make_label("takao2014", "100 m"),
             cite="takao2014", side="both", pixel=100,
             variant="pixel_size=100",
             cls="openquake.fdha.secondary_surf_rup.takao2014.Takao2014SecondarySR",
             eval=lambda: eval_takao2014(100),
             color="#8c564b", ls=":", marker="D",
             mw_key="takao2014", r_key="takao2014"),
    ])
    for pix in (100, 500):
        for comb in ("A", "B", "C"):
            panels["b"].append(dict(
                model_id=f"visini2025_reverse_{pix}_{comb}",
                label=make_label("visini2025",
                                 f"{pix} m, combination {comb}"),
                cite="visini2025", side="both", pixel=pix,
                variant=f"combination={comb}",
                cls="openquake.fdha.secondary_surf_rup.visini2025.Visini2025SecondarySR",
                eval=(lambda p=pix, c=comb:
                      eval_visini("reverse", REF_MAG, p, c)),
                color=VISINI_COLOR[pix], ls=COMB_LS[comb], marker="^",
                mw_key="visini2025_reverse", r_key="visini2025"))

    # ---- Panel (c) Strike-Slip ---------------------------------------------
    panels["c"].extend([
        dict(model_id="petersen2011_ss_100",
             label=make_label("petersen2011", "100 m"),
             cite="petersen2011", side="positive", pixel=100,
             variant="cell_size=100",
             cls="openquake.fdha.secondary_surf_rup.petersen2011.Petersen2011SecondarySR",
             eval=lambda: eval_petersen2011_ss(100),
             color="#1f77b4", ls="-", marker="o",
             mw_key="petersen2011_secondary", r_key="petersen2011_secondary"),
        dict(model_id="takao2013_ss_500",
             label=make_label("takao2013", "500 m"),
             cite="takao2013", side="positive", pixel=500,
             variant="pixel_size=500",
             cls="openquake.fdha.secondary_surf_rup.takao2013.Takao2013SecondarySR",
             eval=lambda: eval_takao2013_ss(REF_MAG, pixel_size=500),
             color="#ff7f0e", ls="-", marker="s",
             mw_key="takao2013_secondary", r_key="takao2013_secondary"),
        dict(model_id="takao2014_ss_100",
             label=make_label("takao2014", "100 m"),
             cite="takao2014", side="positive", pixel=100,
             variant="pixel_size=100",
             cls="openquake.fdha.secondary_surf_rup.takao2014.Takao2014SecondarySR",
             eval=lambda: eval_takao2014_ss(100),
             color="#bcbd22", ls="--", marker="D",
             mw_key="takao2014", r_key="takao2014"),
        dict(model_id="rodriguez2023_ss_1",
             label=make_label("rodriguez2023", "1 m"),
             cite="rodriguez2023", side="positive", pixel=1,
             variant="pixel_size=1",
             cls="openquake.fdha.secondary_surf_rup.rodriguez2023.Rodriguez2023SecondarySR",
             eval=lambda: eval_rodriguez2023_ss(),
             color="#2ca02c", ls=":", marker="^",
             mw_key="rodriguez2023", r_key="rodriguez2023"),
    ])
    return panels


PANEL_TITLES = {
    "a": "Normal Faulting",
    "b": "Reverse Faulting",
    "c": "Strike-Slip Faulting",
}


def _draw_panel(ax, letter, curves, x_signed: bool, ymax: float = 1.0,
                legend_ncol: int = 1):
    seen = set()
    for spec in curves:
        if spec["model_id"] in seen:
            stop(spec["model_id"], ValueError(
                f"duplicate model_id '{spec['model_id']}' in panel"))
        seen.add(spec["model_id"])
        validate_label(spec["model_id"], spec["label"])
        try:
            r, p = spec["eval"]()
        except Exception as exc:  # noqa: BLE001
            stop(spec["model_id"], exc)
        if not np.all(np.isfinite(p)):
            stop(spec["model_id"], ValueError("non-finite probabilities"))
        p = np.clip(p, 0.0, 1.0)
        ax.plot(r, p, color=spec["color"], linestyle=spec["ls"],
                marker=spec["marker"], markevery=MARKER_EVERY,
                markersize=6, linewidth=2.4, label=spec["label"])
    if x_signed:
        ax.set_xlim(-10.0, 10.0)
        ax.axvline(0.0, color="0.5", linewidth=0.6, linestyle="-",
                   alpha=0.5)
    else:
        ax.set_xlim(0.0, 10.0)
    ax.set_yscale("linear")
    ax.set_ylim(0.0, ymax)
    ax.set_xlabel(r"Distance from principal trace, $r$ (km)")
    ax.set_ylabel("Probability of Distributed Surface Rupture")
    ax.set_title(PANEL_TITLES[letter])
    ax.text(-0.18, 1.05, f"({letter})", transform=ax.transAxes,
            fontsize=14, fontweight="bold", ha="left", va="bottom")
    ax.grid(True, which="major", color="0.85", linewidth=0.5)
    ax.grid(False, which="minor")
    # Legend placed to the right of the panel, single column.
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
              framealpha=0.9, borderaxespad=0.0)


def main():
    panels = build_panels()

    # Global duplicate-model_id check.
    all_ids = {}
    for panel, curves in panels.items():
        for spec in curves:
            mid = spec["model_id"]
            if mid in all_ids:
                stop(mid, ValueError(
                    f"duplicate model_id across panels (also in "
                    f"panel {all_ids[mid]})"))
            all_ids[mid] = panel

    fig, axes = plt.subplots(3, 1, figsize=(11, 16))
    _draw_panel(axes[0], "a", panels["a"], x_signed=True, ymax=1.0)
    _draw_panel(axes[1], "b", panels["b"], x_signed=True, ymax=0.5)
    _draw_panel(axes[2], "c", panels["c"], x_signed=False, ymax=0.5)
    fig.tight_layout()

    pdf = OUT_DIR / "figure_03_distributed_surf_rup.pdf"
    png = OUT_DIR / "figure_03_distributed_surf_rup.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=300)
    plt.close(fig)
    print(f"wrote {pdf}")
    print(f"wrote {png}")

    # ---- Metadata CSV -------------------------------------------------------
    csv_path = OUT_DIR / "figure_03_metadata.csv"
    fields = [
        "panel_id", "model_id", "display_name", "faulting_style",
        "side", "pixel_or_cell_size_m", "variant",
        "mw_validity_min", "mw_validity_max",
        "r_validity_min_km", "r_validity_max_km",
        "framework_class_path", "source_citation", "source_doi",
        "exclusion_reason",
    ]
    style_by_panel = {"a": "normal", "b": "reverse", "c": "strike-slip"}
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for letter, curves in panels.items():
            for spec in curves:
                cit, doi = CITATION[spec["cite"]]
                mw_lo, mw_hi = MW_RANGE.get(spec["mw_key"], ("", ""))
                r_fw, r_hw = R_RANGE_KM.get(spec["r_key"], ("", ""))
                w.writerow(dict(
                    panel_id=letter, model_id=spec["model_id"],
                    display_name=spec["label"],
                    faulting_style=style_by_panel[letter],
                    side=spec["side"],
                    pixel_or_cell_size_m=spec["pixel"] or "",
                    variant=spec["variant"],
                    mw_validity_min=mw_lo if mw_lo is not None else "",
                    mw_validity_max=mw_hi if mw_hi is not None else "",
                    r_validity_min_km=-r_fw if spec["side"] == "both" else 0.0,
                    r_validity_max_km=r_hw,
                    framework_class_path=spec["cls"],
                    source_citation=cit, source_doi=doi,
                    exclusion_reason=""))
    print(f"wrote {csv_path}  ({sum(len(c) for c in panels.values())} curves)")


if __name__ == "__main__":
    main()
