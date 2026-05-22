# -*- coding: utf-8 -*-
"""
Generate Figure 4 of the SRL manuscript: median distributed (secondary)
fault displacement versus distance from the principal fault.

Style reference: Valentini et al. (2025), Reviews of Geophysics, Figure 16.
Three subpanels stacked vertically: (a) normal, (b) reverse, (c) strike-slip.
Reference magnitude m = 6.5 for magnitude-dependent models. For dip-slip
panels the curve is drawn continuously across r = 0 (negative r = footwall,
positive r = hanging wall); strike-slip is shown for r >= 0 only.

All curves are produced by calling the model implementations registered in
``openquake.fdha.secondary_surf_displ``. No equations are re-implemented and
no coefficients are hardcoded here. Validity ranges are taken from
docs/UserManual_Enhanced/06-Models.md (sourced from Valentini et al. (2025),
Reviews of Geophysics, Table 4).

User-directed configuration (decisions from the planning Q&A):
- Y-axis: median displacement (m), log scale.
- Default dip for Visini normal/reverse: 60 deg.
- Visini: plot all three combinations A/B/C in both normal and reverse
  panels. Pixel size does not enter the median displacement regression
  (it modulates only the conditional rupture probability, which is plotted
  in Figure 3), so a single curve per combination suffices.
- Moss et al. (2022) included in panel (b); flagged as outside the
  Valentini Fig. 16 inventory.
- Petersen et al. (2011) secondary FD displacement regression (Eq. 18) does
  not depend on cell size, so panel (c) shows a single curve.
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
from scipy.optimize import brentq

# --- Framework model registry --------------------------------------------
from openquake.fdha.secondary_surf_displ.youngs2003 import Youngs2003SecondaryFD
from openquake.fdha.secondary_surf_displ.petersen2011 import (
    Petersen2011SecondaryFD,
)
from openquake.fdha.secondary_surf_displ.visini2025 import Visini2025SecondaryFD
from openquake.fdha.secondary_surf_displ.moss2022 import Moss2022SecondaryFD


OUT_DIR = Path(__file__).resolve().parent
REF_MAG = 6.5
DEFAULT_DIP = 60.0

# -----------------------------------------------------------------------------
# Plot styling -- matches Figures 1-3.
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
# Validity ranges (Valentini et al. 2025 Table 4, via 06-Models.md).
# -----------------------------------------------------------------------------
MW_RANGE = {
    "youngs2003_secondary": (5.5, 7.4),
    "petersen2011_secondary": (6.5, 7.5),
    "visini2025_normal": (5.5, 7.9),
    "visini2025_reverse": (4.9, 7.9),
    "moss2022_secondary": (None, None),  # report-specific
}
# Per-side r limits in km. (fw_max_km, hw_max_km); strike-slip uses hw_max.
R_RANGE_KM = {
    "youngs2003_secondary": (15.0, 15.0),
    "petersen2011_secondary": (2.5, 2.5),
    "visini2025": (8.0, 10.0),
    "moss2022_secondary": (10.0, 10.0),  # report-specific; cap to panel x-limit
}

# -----------------------------------------------------------------------------
# Legend naming -- SRL convention (mirrors Figures 2 and 3).
# -----------------------------------------------------------------------------
CITE_AUTHORS = {
    "youngs2003": ("Youngs", 20, 2003),
    "petersen2011": ("Petersen", 8, 2011),
    "visini2025": ("Visini", 7, 2025),
    "moss2022": ("Moss", 5, 2022),
}

CITATION = {
    "youngs2003": (
        "Youngs, R. R., Arabasz, W. J., Anderson, R. E., et al. (2003). "
        "A methodology for probabilistic fault displacement hazard analysis "
        "(PFDHA). Earthquake Spectra, 19(1), 191-219.",
        "10.1193/1.1542891",
    ),
    "petersen2011": (
        "Petersen, M. D., Dawson, T. E., Chen, R., Cao, T., Wills, C. J., "
        "Schwartz, D. P., & Frankel, A. D. (2011). Fault displacement hazard "
        "for strike-slip faults. Bulletin of the Seismological Society of "
        "America, 101(2), 805-825.",
        "10.1785/0120100035",
    ),
    "visini2025": (
        "Visini, F., Boncio, P., Valentini, A., Scotti, O., Nurminen, F., "
        "Baize, S., & Pace, B. (2025). Empirical regressions for distributed "
        "faulting of dip-slip earthquakes. Earthquake Spectra.",
        "10.1177/87552930241308860",
    ),
    "moss2022": (
        "Moss, R., Thompson, S., Kuo, C.-H., Younesi, K., & Baumont, D. "
        "(2022). Reverse Fault PFDHA. Report GIRS-2022-05.",
        "10.34948/N3F595",
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
        return f"{first} et al. ({year})"  # not used here; placeholder
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
# Distance grids.
# -----------------------------------------------------------------------------
R_STEP_KM = 0.05
MARKER_EVERY = int(round(1.0 / R_STEP_KM))


def signed_grid(fw_max_km: float, hw_max_km: float,
                axis_limit_km: float = 10.0) -> np.ndarray:
    fw = min(fw_max_km, axis_limit_km)
    hw = min(hw_max_km, axis_limit_km)
    n_fw = int(round(fw / R_STEP_KM))
    n_hw = int(round(hw / R_STEP_KM))
    fw_pts = -np.linspace(fw, R_STEP_KM, n_fw) if n_fw > 0 else np.array([])
    hw_pts = np.linspace(R_STEP_KM, hw, n_hw) if n_hw > 0 else np.array([])
    return np.concatenate([fw_pts, hw_pts])


def one_sided_grid(r_max_km: float,
                   axis_limit_km: float = 2.5) -> np.ndarray:
    r_max = min(r_max_km, axis_limit_km)
    n = int(round(r_max / R_STEP_KM))
    return np.linspace(R_STEP_KM, r_max, n)


# -----------------------------------------------------------------------------
# Median accessors (one per model family).
# -----------------------------------------------------------------------------
# Numerical-inversion settings for models that only expose P(D > d).
_INV_D_LO_M = 1e-5     # 0.01 mm
_INV_D_HI_M = 50.0     # 50 m


def _median_from_sf(prob_fn, p_target: float = 0.5):
    """Find d such that P(D > d) = p_target via brentq on log10(d)."""
    def f(log10_d):
        d = 10.0 ** log10_d
        # prob_fn(d) returns scalar
        return float(prob_fn(d)) - p_target

    f_lo = f(np.log10(_INV_D_LO_M))
    f_hi = f(np.log10(_INV_D_HI_M))
    if f_lo < 0:                            # already below target at lo bound
        return float("nan")
    if f_hi > 0:                            # still above target at hi bound
        return float("nan")
    log10_med = brentq(f, np.log10(_INV_D_LO_M), np.log10(_INV_D_HI_M),
                       xtol=1e-5, rtol=1e-6)
    return 10.0 ** log10_med


def eval_youngs2003(percentile: str, mag: float):
    """Median secondary displacement (m) vs signed r (km), normal style."""
    m = Youngs2003SecondaryFD()
    fw_max, hw_max = R_RANGE_KM["youngs2003_secondary"]
    r_km = signed_grid(fw_max, hw_max)
    rx_km = r_km
    r_abs_km = np.abs(r_km)
    med = np.empty_like(r_km)
    for i, (rxi, ri) in enumerate(zip(rx_km, r_abs_km)):
        def sf(d_m, rxi=rxi, ri=ri):
            p = m.get_prob(np.array([d_m]), mag, np.array([rxi]),
                           np.array([ri]), percentile=percentile)
            return float(np.asarray(p).reshape(-1)[0])
        med[i] = _median_from_sf(sf)
    return r_km, med


def eval_petersen2011_ss(mag: float):
    """Median secondary displacement (m) vs positive r (km), strike-slip."""
    m = Petersen2011SecondaryFD()
    r_max = R_RANGE_KM["petersen2011_secondary"][1]
    r_km = one_sided_grid(r_max)
    med = np.empty_like(r_km)
    for i, ri in enumerate(r_km):
        def sf(d_m, ri=ri):
            p = m.get_prob(np.array([d_m]), mag, np.array([ri]))
            return float(np.asarray(p).reshape(-1)[0])
        med[i] = _median_from_sf(sf)
    return r_km, med


def eval_visini(style: str, mag: float, pixel_size: int, combination: str,
                dip: float = DEFAULT_DIP):
    """Median secondary displacement (m) vs signed r (km), normal or reverse.

    Uses ``Visini2025SecondaryFD.get_median_displacement`` directly. TPFm is
    pre-computed from the calling convention's WC1994 scaling, smoothed per
    site over a window proportional to s_km, as in the model's own
    ``get_prob`` path.
    """
    m = Visini2025SecondaryFD()
    fw_max, hw_max = R_RANGE_KM["visini2025"]
    r_km = signed_grid(fw_max, hw_max)
    s_m = np.abs(r_km) * 1000.0          # closest distance, meters
    rx_m = r_km * 1000.0                  # signed cross-fault distance, meters
    med = np.empty_like(r_km)
    for i, (si, rxi, ski) in enumerate(zip(s_m, rx_m, np.abs(r_km))):
        tpfm_i, _ = m.compute_tpfm_from_scaling(
            mag=mag, style=style, model="WC1994",
            norm_pos=0.5, distance=float(ski), dip=dip,
        )
        med[i] = float(m.get_median_displacement(
            mag=mag, s=np.array([si]), rx=np.array([rxi]),
            tpfm=np.array([tpfm_i]), style=style, combination=combination,
        ).reshape(-1)[0])
    # Visini's pixel-size dependence is consumed elsewhere (the conditional
    # rupture probability), not in the median displacement regression. We
    # still tag curves by pixel size for plot-legend bookkeeping but the
    # numeric values for a given combination are identical across pixel
    # sizes; that fact is documented in the caption.
    del pixel_size
    return r_km, med


def eval_moss2022(mag: float, method: str, percentile: str = "50"):
    """Median secondary displacement (m) vs signed r (km), reverse."""
    m = Moss2022SecondaryFD()
    fw_max, hw_max = R_RANGE_KM["moss2022_secondary"]
    r_km = signed_grid(fw_max, hw_max)
    rx_km = r_km
    r_abs_km = np.abs(r_km)
    med = np.empty_like(r_km)
    for i, (rxi, ri) in enumerate(zip(rx_km, r_abs_km)):
        def sf(d_m, rxi=rxi, ri=ri):
            p = m.get_prob(np.array([d_m]), mag, np.array([ri]),
                           np.array([rxi]),
                           method=method, percentile=percentile,
                           version="MD", completeness="complete",
                           sigma_type="recommended", faulting="simple")
            return float(np.asarray(p).reshape(-1)[0])
        med[i] = _median_from_sf(sf)
    return r_km, med


# -----------------------------------------------------------------------------
# Panel specifications.
# -----------------------------------------------------------------------------
COMB_LS = {"A": "-", "B": "--", "C": ":"}
# One colour per combination; pixel-size dimension dropped (see module
# docstring: median regression has no pixel-size term).
VISINI_COMB_COLOR = {"A": "#d62728", "B": "#9467bd", "C": "#e377c2"}


def build_panels():
    panels = {"a": [], "b": [], "c": []}

    # ---- Panel (a) Normal ---------------------------------------------------
    panels["a"].extend([
        dict(model_id="youngs2003_normal_85",
             label=make_label("youngs2003", "85th percentile"),
             cite="youngs2003", side="both", pixel=None,
             variant="percentile=85",
             cls="openquake.fdha.secondary_surf_displ.youngs2003.Youngs2003SecondaryFD",
             eval=lambda: eval_youngs2003("85", REF_MAG),
             color="#1f77b4", ls="-", marker="o",
             mw_key="youngs2003_secondary", r_key="youngs2003_secondary"),
        dict(model_id="youngs2003_normal_95",
             label=make_label("youngs2003", "95th percentile"),
             cite="youngs2003", side="both", pixel=None,
             variant="percentile=95",
             cls="openquake.fdha.secondary_surf_displ.youngs2003.Youngs2003SecondaryFD",
             eval=lambda: eval_youngs2003("95", REF_MAG),
             color="#17becf", ls="--", marker="v",
             mw_key="youngs2003_secondary", r_key="youngs2003_secondary"),
    ])
    for comb in ("A", "B", "C"):
        panels["a"].append(dict(
            model_id=f"visini2025_normal_{comb}",
            label=make_label("visini2025", f"combination {comb}"),
            cite="visini2025", side="both", pixel=None,
            variant=f"combination={comb}",
            cls="openquake.fdha.secondary_surf_displ.visini2025.Visini2025SecondaryFD",
            eval=(lambda c=comb: eval_visini("normal", REF_MAG, 500, c)),
            color=VISINI_COMB_COLOR[comb], ls=COMB_LS[comb], marker="^",
            mw_key="visini2025_normal", r_key="visini2025"))

    # ---- Panel (b) Reverse --------------------------------------------------
    for comb in ("A", "B", "C"):
        panels["b"].append(dict(
            model_id=f"visini2025_reverse_{comb}",
            label=make_label("visini2025", f"combination {comb}"),
            cite="visini2025", side="both", pixel=None,
            variant=f"combination={comb}",
            cls="openquake.fdha.secondary_surf_displ.visini2025.Visini2025SecondaryFD",
            eval=(lambda c=comb: eval_visini("reverse", REF_MAG, 500, c)),
            color=VISINI_COMB_COLOR[comb], ls=COMB_LS[comb], marker="^",
            mw_key="visini2025_reverse", r_key="visini2025"))
    panels["b"].extend([
        dict(model_id="moss2022_reverse_gamma",
             label=make_label("moss2022", "gamma, 50th percentile"),
             cite="moss2022", side="both", pixel=None,
             variant="method=gamma; percentile=50",
             cls="openquake.fdha.secondary_surf_displ.moss2022.Moss2022SecondaryFD",
             eval=lambda: eval_moss2022(REF_MAG, "gamma", "50"),
             color="#2ca02c", ls="-", marker="s",
             mw_key="moss2022_secondary", r_key="moss2022_secondary"),
        dict(model_id="moss2022_reverse_envelope",
             label=make_label("moss2022", "envelope, 50th percentile"),
             cite="moss2022", side="both", pixel=None,
             variant="method=envelope; percentile=50",
             cls="openquake.fdha.secondary_surf_displ.moss2022.Moss2022SecondaryFD",
             eval=lambda: eval_moss2022(REF_MAG, "envelope", "50"),
             color="#98df8a", ls="--", marker="D",
             mw_key="moss2022_secondary", r_key="moss2022_secondary"),
    ])

    # ---- Panel (c) Strike-Slip ---------------------------------------------
    panels["c"].append(dict(
        model_id="petersen2011_ss",
        label=make_label("petersen2011"),
        cite="petersen2011", side="positive", pixel=None,
        variant="cell_size-independent (Eq. 18)",
        cls="openquake.fdha.secondary_surf_displ.petersen2011.Petersen2011SecondaryFD",
        eval=lambda: eval_petersen2011_ss(REF_MAG),
        color="#1f77b4", ls="-", marker="o",
        mw_key="petersen2011_secondary", r_key="petersen2011_secondary"))
    return panels


PANEL_TITLES = {
    "a": "Normal Faulting",
    "b": "Reverse Faulting",
    "c": "Strike-Slip Faulting",
}


def _draw_panel(ax, letter, curves, x_signed: bool, x_max: float,
                y_lo: float = 1e-3, y_hi: float = 1.0):
    seen = set()
    for spec in curves:
        if spec["model_id"] in seen:
            stop(spec["model_id"], ValueError(
                f"duplicate model_id '{spec['model_id']}' in panel"))
        seen.add(spec["model_id"])
        validate_label(spec["model_id"], spec["label"])
        try:
            r, y = spec["eval"]()
        except Exception as exc:  # noqa: BLE001
            stop(spec["model_id"], exc)
        # Filter non-finite points individually so a curve survives any
        # numerical-inversion gaps.
        mask = np.isfinite(y) & (y > 0)
        if not np.any(mask):
            stop(spec["model_id"], ValueError("no finite/positive medians"))
        ax.plot(r[mask], y[mask], color=spec["color"], linestyle=spec["ls"],
                marker=spec["marker"], markevery=MARKER_EVERY,
                markersize=6, linewidth=2.4, label=spec["label"])
    if x_signed:
        ax.set_xlim(-x_max, x_max)
        ax.axvline(0.0, color="0.5", linewidth=0.6, linestyle="-",
                   alpha=0.5)
    else:
        ax.set_xlim(0.0, x_max)
    ax.set_yscale("log")
    ax.set_ylim(y_lo, y_hi)
    ax.set_xlabel(r"Distance from principal trace, $r$ (km)")
    ax.set_ylabel("Distributed displacement (m)")
    ax.set_title(PANEL_TITLES[letter])
    ax.text(-0.18, 1.05, f"({letter})", transform=ax.transAxes,
            fontsize=14, fontweight="bold", ha="left", va="bottom")
    ax.grid(True, which="major", color="0.85", linewidth=0.5)
    ax.grid(True, which="minor", color="0.94", linewidth=0.4)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
              framealpha=0.9, borderaxespad=0.0)


def main():
    panels = build_panels()

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
    _draw_panel(axes[0], "a", panels["a"], x_signed=True, x_max=10.0)
    _draw_panel(axes[1], "b", panels["b"], x_signed=True, x_max=10.0)
    _draw_panel(axes[2], "c", panels["c"], x_signed=False, x_max=2.5)
    fig.tight_layout()

    pdf = OUT_DIR / "figure_04_distributed_surf_displ.pdf"
    png = OUT_DIR / "figure_04_distributed_surf_displ.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=300)
    plt.close(fig)
    print(f"wrote {pdf}")
    print(f"wrote {png}")

    csv_path = OUT_DIR / "figure_04_metadata.csv"
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
