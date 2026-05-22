# -*- coding: utf-8 -*-
"""
Generate Figure 2 of the SRL manuscript: principal fault displacement models.

Subfigure A (3x2): normalized displacement profiles D/MD (row 1) and D/AD
(row 2) vs x/L on [0, 0.5], by faulting style (normal, reverse, strike-slip).
Subfigure B (3x1): absolute median displacement profiles at m = 6.5 on a log
y-axis.

Scientific reference
--------------------
Valentini, A., et al. (2025). Twenty-Five Years of Probabilistic Fault
Displacement Hazard Assessment. Reviews of Geophysics, 63(3).
https://doi.org/10.1029/2024RG000875  (Figures 13 and 14; Tables 3-4)

Data provenance
---------------
Every curve is produced by calling model implementations in
``openquake.fdha.primary_surf_displ``. No equations or coefficients are
reimplemented here. For Moss et al. (2024) the source-defined normalized
helpers (``get_prob_D_AD`` / ``get_prob_D_MD``) were added directly to the
framework class in this branch; the median is then obtained by numerical
inversion of those public methods. For Takao et al. (2013) the framework
already exposes the same helpers. For Subfigure B, the median absolute
displacement at m = 6.5 is obtained by numerical inversion of each model's
public ``get_prob(d, X_L_ratio, mag)``.

Known framework gaps (documented in figure_02_metadata.csv, not silently
skipped)
-------------------------------------------------------------------------
1. Moss and Ross (2011) primary fault displacement model is not implemented
   in the framework (import commented out in
   ``openquake.fdha.primary_surf_displ.__init__``). Panels (b) and (e)
   exclude the three Moss-and-Ross entries.
2. Petersen et al. (2011) is, in its source publication, an absolute log-
   normal model -- the source does not define a D/MD or D/AD distribution.
   Subfigure A panels (c) and (f) therefore exclude the three Petersen
   functional forms; they appear in Subfigure B panel (i) where they belong.
3. Framework model classes carry no ``display_name`` attribute and no
   validity-range metadata. Legend labels are assembled in-script from an
   authoritative author table (CITE_AUTHORS) and validated against the SRL
   naming convention; validity ranges are taken from the project user manual
   docs/UserManual_Enhanced/06-Models.md (sourced from Valentini et al.
   (2025), Reviews of Geophysics, Table 4).
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

# Framework model registry --------------------------------------------------
from openquake.fdha.primary_surf_displ.youngs2003 import Youngs2003PrimaryFD
from openquake.fdha.primary_surf_displ.petersen2011 import (
    Petersen2011PrimaryFD_bilinear,
    Petersen2011PrimaryFD_elliptical,
    Petersen2011PrimaryFD_quadratic,
)
from openquake.fdha.primary_surf_displ.moss2024 import Moss2024PrimaryFD
from openquake.fdha.primary_surf_displ.moss_ross2011 import MossRoss2011PrimaryFD
from openquake.fdha.primary_surf_displ.takao2013 import Takao2013PrimaryFD
from openquake.fdha.primary_surf_displ.lavrentiadis2023 import (
    Lavrentiadis2023PrimaryFD,
)
from openquake.fdha.primary_surf_displ.kuehn2024.kuehn2024 import (
    Kuehn2024PrimaryFD,
)
from openquake.fdha.primary_surf_displ.chiou2025 import Chiou2025PrimaryFD


OUT_DIR = Path(__file__).resolve().parent
REF_MAG_B = 6.5  # reference magnitude for Subfigure B

# -----------------------------------------------------------------------------
# SRL Art Guidelines compliant matplotlib initialization (per task spec).
# -----------------------------------------------------------------------------
rcParams.update({
    # Match Figure 1: Times New Roman serif, size 12 baseline, italics for
    # mathvariables only. Times New Roman is installed at
    # ~/.local/share/fonts/Times.TTF (plus bold/italic faces).
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "font.size": 12,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    # SRL rule rules/axis lines 0.5-1.0 pt; data curves are made thicker
    # for visibility (matches Figure 1).
    "axes.linewidth": 0.8,
    "grid.linewidth": 0.6,
    "lines.linewidth": 2.4,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    # SSA Art Guidelines: embed fonts as TrueType (Type 42) in PDF/EPS.
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    # Mathtext renders italic math variables (m, x, L, D, MD, AD) in
    # Times New Roman to match the body text (custom mathtext).
    "mathtext.fontset": "custom",
    "mathtext.rm": "Times New Roman",
    "mathtext.it": "Times New Roman:italic",
    "mathtext.bf": "Times New Roman:bold",
})

# -----------------------------------------------------------------------------
# Validity ranges (Valentini et al. 2025, Rev. Geophys., Table 4 -- as
# tabulated in docs/UserManual_Enhanced/06-Models.md). Used to gate Subfigure
# B by magnitude and to populate the metadata CSV.
# -----------------------------------------------------------------------------
MW_RANGE = {
    "youngs2003": (4.5, 7.6),
    "petersen2011": (6.0, 8.0),
    "takao2013": (5.8, 7.4),
    "moss2024": (4.7, 8.0),
    "lavrentiadis2023": (5.0, 8.5),
    "kuehn2024_reverse": (5.0, 8.0),
    "kuehn2024_normal": (6.0, 8.0),
    "kuehn2024_strikeslip": (6.0, 8.0),
    "chiou2025": (6.0, 8.3),
    "mossross2011": (5.5, 8.0),  # not implemented but listed for the CSV
}

# All implemented models accept x/L on [0, 1] (folded symmetrically about 0.5).
XL_RANGE = (0.0, 0.5)

# -----------------------------------------------------------------------------
# Citation table and legend naming convention.
# -----------------------------------------------------------------------------
CITE_AUTHORS = {
    # key -> (first author surname, total author count, year)
    "youngs2003": ("Youngs", 20, 2003),
    "petersen2011": ("Petersen", 8, 2011),
    "takao2013": ("Takao", 4, 2013),
    "mossross2011": ("Moss", 2, 2011),         # Moss and Ross
    "moss2024": ("Moss", 5, 2024),
    "lavrentiadis2023": ("Lavrentiadis", 2, 2023),  # Lavrentiadis and Abrahamson
    "kuehn2024": ("Kuehn", 5, 2024),
    "chiou2025": ("Chiou", 6, 2025),
}
SECOND_AUTHOR = {"mossross2011": "Ross", "lavrentiadis2023": "Abrahamson"}

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
    "takao2013": (
        "Takao, M., Tsuchiyama, J., Annaka, T., & Kurita, T. (2013). "
        "Application of probabilistic fault displacement hazard analysis in "
        "Japan. Journal of Japan Association for Earthquake Engineering, "
        "13(1), 17-36.",
        "10.5610/jaee.13.1_17",
    ),
    "mossross2011": (
        "Moss, R. E. S., & Ross, Z. E. (2011). Probabilistic fault "
        "displacement hazard analysis for reverse faults. Bulletin of the "
        "Seismological Society of America, 101(4), 1542-1553.",
        "10.1785/0120100248",
    ),
    "moss2024": (
        "Moss, R. E. S., Thompson, S. C., Kuo, C.-H., Younesi, K., & "
        "Baumont, D. (2024). New probabilistic fault displacement hazard "
        "models for reverse faulting. Earthquake Spectra.",
        "10.1177/87552930241288560",
    ),
    "lavrentiadis2023": (
        "Lavrentiadis, G., & Abrahamson, N. A. (2023). A non-ergodic "
        "spectral acceleration ground motion model for California developed "
        "with random vibration theory. Bulletin of Earthquake Engineering, "
        "21, 5265-5291.",
        "10.1007/s10518-023-01773-0",
    ),
    "kuehn2024": (
        "Kuehn, N. M., Kottke, A. R., Sarmiento, A. C., Madugo, C. M., "
        "& Bozorgnia, Y. (2024). A fault displacement model based on the "
        "FDHI database. Earthquake Spectra.",
        "10.1177/87552930241291077",
    ),
    "chiou2025": (
        "Chiou, B. S. J., Chen, R., Thomas, K., Milliner, C., Dawson, T., "
        "& Petersen, M. (2025). Fault displacement model for surface "
        "principal rupture of strike-slip faults. Earthquake Spectra.",
        "10.1177/87552930251337703",
    ),
}

# SRL legend regex: "Surname[ et al.| and Surname] (YYYY)[, variant]?"
_LABEL_RE = re.compile(
    r"^[A-Z][\w.'‐-]*"
    r"(?: et al\.| and [A-Z][\w.'‐-]*)?"
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
        stop(model_id, ValueError(
            f"legend label uses '&' (spell out 'and'): {label!r}"))
    if "$" in label or r"\it" in label:
        stop(model_id, ValueError(
            f"legend label contains italic markup: {label!r}"))
    if not _LABEL_RE.match(label):
        stop(model_id, ValueError(
            f"malformed legend label: {label!r}"))


# -----------------------------------------------------------------------------
# Median utilities -- numerical inversion of the framework's get_prob.
# -----------------------------------------------------------------------------
def _bracketed_brentq(f, lo, hi, expand_hi_to=1e3, expand_factor=10.0):
    """Find root of f within [lo, hi]; expand hi geometrically if needed."""
    f_lo = f(lo)
    if not np.isfinite(f_lo):
        raise ValueError("non-finite f(lo) bracketing median")
    h = hi
    while h <= expand_hi_to:
        f_hi = f(h)
        if not np.isfinite(f_hi):
            raise ValueError(f"non-finite f({h}) bracketing median")
        if f_lo * f_hi < 0:
            return brentq(f, lo, h, xtol=1e-8)
        h *= expand_factor
    raise ValueError(
        f"could not bracket median: f({lo})={f_lo}, f({h/expand_factor})={f_hi}")


def median_from_sf(sf_fn, lo=1e-9, hi=1.0):
    """Median of a positive RV via root of (SF - 0.5)."""
    return _bracketed_brentq(lambda v: float(sf_fn(v)) - 0.5, lo, hi)


def median_absolute(model_get_prob, x_L, mag, lo=1e-5, hi=2.0, **kw):
    """Median absolute displacement (meters) at a single x/L via inversion."""
    def f(d):
        v = model_get_prob(d, np.atleast_1d(x_L), mag, **kw)
        v = np.atleast_1d(np.asarray(v, dtype=float)).ravel()
        return float(v[0]) - 0.5
    return _bracketed_brentq(f, lo, hi, expand_hi_to=500.0)


# -----------------------------------------------------------------------------
# Per-curve evaluators -- return (xL grid, median y-array).
# -----------------------------------------------------------------------------
XL_GRID = np.linspace(0.0, 0.5, 21)  # 0.025 step -> markers every 0.05


def eval_normalized(sf_at_xl, *, lo, hi):
    """Median normalized D/XD at each x/L via numerical inversion."""
    out = np.empty_like(XL_GRID)
    for i, x in enumerate(XL_GRID):
        out[i] = median_from_sf(lambda v: sf_at_xl(v, x), lo=lo, hi=hi)
    return XL_GRID, out


def eval_absolute(model_get_prob, mag, **kw):
    """Median absolute displacement (cm) at each x/L for Subfigure B."""
    out = np.empty_like(XL_GRID)
    for i, x in enumerate(XL_GRID):
        d_m = median_absolute(model_get_prob, x, mag, **kw)
        out[i] = d_m * 100.0  # meters -> centimeters
    return XL_GRID, out


# -----------------------------------------------------------------------------
# Curve specifications.
# -----------------------------------------------------------------------------
def build_subfigure_A():
    y = Youngs2003PrimaryFD()
    m24 = Moss2024PrimaryFD()
    mr = MossRoss2011PrimaryFD()
    tk = Takao2013PrimaryFD()
    # Takao srl regimes
    SRL_LONG, SRL_SHORT = 20.0, 5.0  # km, representative of >= 10 / < 10

    panels = {
        # row 1: D/MD
        "a": [],  # Normal D/MD
        "b": [],  # Reverse D/MD
        "c": [],  # Strike-Slip D/MD
        # row 2: D/AD
        "d": [],
        "e": [],
        "f": [],
    }

    # Panel (a) Normal D/MD  -- Youngs et al. (2003), Pezzopane & Dawson data.
    # The framework's get_prob_D_MD always uses the source-defined beta
    # distribution; it is not a selectable option, so it is not in the label.
    panels["a"].append(dict(
        model_id="youngs2003_md",
        label=make_label("youngs2003", "Pezzopane and Dawson data"),
        cite="youngs2003", style="normal", dist="Beta", pred="D/MD",
        cls="openquake.fdha.primary_surf_displ.youngs2003.Youngs2003PrimaryFD",
        eval=lambda: eval_normalized(
            lambda v, x: y.get_prob_D_MD(v, x).reshape(-1)[0],
            lo=1e-9, hi=1.0 - 1e-9),
        color="#1f77b4", ls="-", marker="o",
    ))

    # Panel (b) Reverse D/MD  -- Moss & Ross (2011) Beta + Moss et al. (2024)
    panels["b"].extend([
        dict(model_id="mossross2011_md_beta",
             label=make_label("mossross2011", "Beta on D/MD"),
             cite="mossross2011", style="reverse", dist="Beta", pred="D/MD",
             cls="openquake.fdha.primary_surf_displ.moss_ross2011.MossRoss2011PrimaryFD",
             eval=lambda: eval_normalized(
                 lambda v, x: float(mr.get_prob_D_MD(v, x).reshape(-1)[0]),
                 lo=1e-9, hi=1.0 - 1e-9),
             color="#9467bd", ls="--", marker="D"),
        dict(model_id="moss2024_md_gamma",
             label=make_label("moss2024"),
             cite="moss2024", style="reverse", dist="Gamma", pred="D/MD",
             cls="openquake.fdha.primary_surf_displ.moss2024.Moss2024PrimaryFD",
             eval=lambda: eval_normalized(
                 lambda v, x: float(m24.get_prob_D_MD(v, x).reshape(-1)[0]),
                 lo=1e-9, hi=1.0 - 1e-9),
             color="#1f77b4", ls="-", marker="o"),
    ])

    # Panel (c) Strike-Slip D/MD  -- Takao 2 regimes (Petersen excluded)
    panels["c"].extend([
        dict(model_id="takao2013_md_long",
             label=make_label("takao2013", "L >= 10 km"),
             cite="takao2013", style="strike-slip", dist="Beta", pred="D/MD",
             cls="openquake.fdha.primary_surf_displ.takao2013.Takao2013PrimaryFD",
             eval=lambda: eval_normalized(
                 lambda v, x: float(tk.get_prob_D_MD(v, x, SRL_LONG)),
                 lo=1e-9, hi=1.0 - 1e-9),
             color="#1f77b4", ls="-", marker="o"),
        dict(model_id="takao2013_md_short",
             label=make_label("takao2013", "L < 10 km"),
             cite="takao2013", style="strike-slip", dist="Beta", pred="D/MD",
             cls="openquake.fdha.primary_surf_displ.takao2013.Takao2013PrimaryFD",
             eval=lambda: eval_normalized(
                 lambda v, x: float(tk.get_prob_D_MD(v, x, SRL_SHORT)),
                 lo=1e-9, hi=1.0 - 1e-9),
             color="#2ca02c", ls="--", marker="s"),
    ])

    # Panel (d) Normal D/AD  -- Youngs et al. (2003), Pezzopane & Dawson data.
    # get_prob_D_AD always uses the source-defined gamma distribution.
    panels["d"].append(dict(
        model_id="youngs2003_ad",
        label=make_label("youngs2003", "Pezzopane and Dawson data"),
        cite="youngs2003", style="normal", dist="Gamma", pred="D/AD",
        cls="openquake.fdha.primary_surf_displ.youngs2003.Youngs2003PrimaryFD",
        eval=lambda: eval_normalized(
            lambda v, x: float(y.get_prob_D_AD(v, x).reshape(-1)[0]),
            lo=1e-9, hi=5.0),
        color="#1f77b4", ls="-", marker="o",
    ))

    # Panel (e) Reverse D/AD  -- Moss & Ross (2011) Gamma + Weibull + Moss 2024 + Takao.
    panels["e"].extend([
        dict(model_id="mossross2011_ad_gamma",
             label=make_label("mossross2011", "Gamma on D/AD"),
             cite="mossross2011", style="reverse", dist="Gamma", pred="D/AD",
             cls="openquake.fdha.primary_surf_displ.moss_ross2011.MossRoss2011PrimaryFD",
             eval=lambda: eval_normalized(
                 lambda v, x: float(
                     mr.get_prob_D_AD(v, x, variant="gamma").reshape(-1)[0]),
                 lo=1e-9, hi=5.0),
             color="#9467bd", ls="--", marker="D"),
        dict(model_id="mossross2011_ad_weibull",
             label=make_label("mossross2011", "Weibull on D/AD"),
             cite="mossross2011", style="reverse", dist="Weibull", pred="D/AD",
             cls="openquake.fdha.primary_surf_displ.moss_ross2011.MossRoss2011PrimaryFD",
             eval=lambda: eval_normalized(
                 lambda v, x: float(
                     mr.get_prob_D_AD(v, x, variant="weibull").reshape(-1)[0]),
                 lo=1e-9, hi=5.0),
             color="#8c564b", ls=":", marker="X"),
        dict(model_id="moss2024_ad_gamma",
             label=make_label("moss2024"),
             cite="moss2024", style="reverse", dist="Gamma", pred="D/AD",
             cls="openquake.fdha.primary_surf_displ.moss2024.Moss2024PrimaryFD",
             eval=lambda: eval_normalized(
                 lambda v, x: float(m24.get_prob_D_AD(v, x).reshape(-1)[0]),
                 lo=1e-9, hi=5.0),
             color="#1f77b4", ls="-", marker="o"),
        dict(model_id="takao2013_ad_long",
             label=make_label("takao2013", "L >= 10 km"),
             cite="takao2013", style="reverse", dist="Gamma", pred="D/AD",
             cls="openquake.fdha.primary_surf_displ.takao2013.Takao2013PrimaryFD",
             eval=lambda: eval_normalized(
                 lambda v, x: float(tk.get_prob_D_AD(v, x, SRL_LONG)),
                 lo=1e-9, hi=5.0),
             color="#ff7f0e", ls="-", marker="^"),
        dict(model_id="takao2013_ad_short",
             label=make_label("takao2013", "L < 10 km"),
             cite="takao2013", style="reverse", dist="Gamma", pred="D/AD",
             cls="openquake.fdha.primary_surf_displ.takao2013.Takao2013PrimaryFD",
             eval=lambda: eval_normalized(
                 lambda v, x: float(tk.get_prob_D_AD(v, x, SRL_SHORT)),
                 lo=1e-9, hi=5.0),
             color="#2ca02c", ls="--", marker="s"),
    ])

    # Panel (f) Strike-Slip D/AD  -- Petersen 2011 (3 variants) + Takao
    # Petersen et al. (2011) Eqs. 14-17 publish source-defined lognormal
    # ln(D/AD) distributions for bilinear, quadratic, and elliptical forms;
    # the framework helpers were added to petersen2011.py in this branch.
    pet_bil_f = Petersen2011PrimaryFD_bilinear()
    pet_ell_f = Petersen2011PrimaryFD_elliptical()
    pet_quad_f = Petersen2011PrimaryFD_quadratic()
    panels["f"].extend([
        dict(model_id="petersen2011_bilinear_ad",
             label=make_label("petersen2011", "Bilinear"),
             cite="petersen2011", style="strike-slip", dist="Lognormal",
             pred="D/AD",
             cls="openquake.fdha.primary_surf_displ.petersen2011.Petersen2011PrimaryFD_bilinear",
             eval=lambda: eval_normalized(
                 lambda v, x: float(pet_bil_f.get_prob_D_AD(v, x)),
                 lo=1e-9, hi=5.0),
             color="#9467bd", ls="-", marker="o"),
        dict(model_id="petersen2011_elliptical_ad",
             label=make_label("petersen2011", "Elliptical"),
             cite="petersen2011", style="strike-slip", dist="Lognormal",
             pred="D/AD",
             cls="openquake.fdha.primary_surf_displ.petersen2011.Petersen2011PrimaryFD_elliptical",
             eval=lambda: eval_normalized(
                 lambda v, x: float(pet_ell_f.get_prob_D_AD(v, x)),
                 lo=1e-9, hi=5.0),
             color="#8c564b", ls="--", marker="D"),
        dict(model_id="petersen2011_quadratic_ad",
             label=make_label("petersen2011", "Quadratic"),
             cite="petersen2011", style="strike-slip", dist="Lognormal",
             pred="D/AD",
             cls="openquake.fdha.primary_surf_displ.petersen2011.Petersen2011PrimaryFD_quadratic",
             eval=lambda: eval_normalized(
                 lambda v, x: float(pet_quad_f.get_prob_D_AD(v, x)),
                 lo=1e-9, hi=5.0),
             color="#e377c2", ls=":", marker="X"),
        dict(model_id="takao2013_ad_long_ss",
             label=make_label("takao2013", "L >= 10 km"),
             cite="takao2013", style="strike-slip", dist="Gamma", pred="D/AD",
             cls="openquake.fdha.primary_surf_displ.takao2013.Takao2013PrimaryFD",
             eval=lambda: eval_normalized(
                 lambda v, x: float(tk.get_prob_D_AD(v, x, SRL_LONG)),
                 lo=1e-9, hi=5.0),
             color="#ff7f0e", ls="-", marker="^"),
        dict(model_id="takao2013_ad_short_ss",
             label=make_label("takao2013", "L < 10 km"),
             cite="takao2013", style="strike-slip", dist="Gamma", pred="D/AD",
             cls="openquake.fdha.primary_surf_displ.takao2013.Takao2013PrimaryFD",
             eval=lambda: eval_normalized(
                 lambda v, x: float(tk.get_prob_D_AD(v, x, SRL_SHORT)),
                 lo=1e-9, hi=5.0),
             color="#2ca02c", ls="--", marker="s"),
    ])

    return panels


def build_subfigure_B():
    lv = Lavrentiadis2023PrimaryFD()
    ku = Kuehn2024PrimaryFD()
    ch = Chiou2025PrimaryFD()
    pet_bil = Petersen2011PrimaryFD_bilinear()
    pet_ell = Petersen2011PrimaryFD_elliptical()
    pet_quad = Petersen2011PrimaryFD_quadratic()

    panels = {"g": [], "h": [], "i": []}

    # (g) Normal
    panels["g"].extend([
        dict(model_id="lav2023_normal",
             label=make_label("lavrentiadis2023"),
             cite="lavrentiadis2023", style="normal",
             dist="aggregate-primary", pred="absolute_cm",
             cls="openquake.fdha.primary_surf_displ.lavrentiadis2023.Lavrentiadis2023PrimaryFD",
             eval=lambda: eval_absolute(
                 lambda d, x, m: lv.get_prob(d, x, m, style="normal"),
                 mag=REF_MAG_B),
             color="#1f77b4", ls="-", marker="o"),
        dict(model_id="kuehn2024_normal",
             label=make_label("kuehn2024"),
             cite="kuehn2024", style="normal",
             dist="aggregate", pred="absolute_cm",
             cls="openquake.fdha.primary_surf_displ.kuehn2024.kuehn2024.Kuehn2024PrimaryFD",
             eval=lambda: eval_absolute(
                 lambda d, x, m: ku.get_prob(d, x, m, style="normal",
                                             epistemic_uncertainty=False),
                 mag=REF_MAG_B),
             color="#d62728", ls="--", marker="s"),
    ])

    # (h) Reverse
    panels["h"].extend([
        dict(model_id="lav2023_reverse",
             label=make_label("lavrentiadis2023"),
             cite="lavrentiadis2023", style="reverse",
             dist="aggregate-primary", pred="absolute_cm",
             cls="openquake.fdha.primary_surf_displ.lavrentiadis2023.Lavrentiadis2023PrimaryFD",
             eval=lambda: eval_absolute(
                 lambda d, x, m: lv.get_prob(d, x, m, style="reverse"),
                 mag=REF_MAG_B),
             color="#1f77b4", ls="-", marker="o"),
        dict(model_id="kuehn2024_reverse",
             label=make_label("kuehn2024"),
             cite="kuehn2024", style="reverse",
             dist="aggregate", pred="absolute_cm",
             cls="openquake.fdha.primary_surf_displ.kuehn2024.kuehn2024.Kuehn2024PrimaryFD",
             eval=lambda: eval_absolute(
                 lambda d, x, m: ku.get_prob(d, x, m, style="reverse",
                                             epistemic_uncertainty=False),
                 mag=REF_MAG_B),
             color="#d62728", ls="--", marker="s"),
        # Moss et al. (2024) is intentionally omitted from panel (h): the
        # source publication predicts only normalized D/MD and D/AD (no
        # absolute D form), matching Valentini et al. (2025) Figure 14 which
        # does not include Moss et al. (2024) in the reverse absolute panel.
    ])

    # (i) Strike-Slip
    panels["i"].extend([
        dict(model_id="petersen2011_bilinear",
             label=make_label("petersen2011", "Bilinear"),
             cite="petersen2011", style="strike-slip",
             dist="single-principal", pred="absolute_cm",
             cls="openquake.fdha.primary_surf_displ.petersen2011.Petersen2011PrimaryFD_bilinear",
             eval=lambda: eval_absolute(
                 lambda d, x, m: pet_bil.get_prob(d, x, m, version="bilinear"),
                 mag=REF_MAG_B),
             color="#9467bd", ls="-", marker="o"),
        dict(model_id="petersen2011_elliptical",
             label=make_label("petersen2011", "Elliptical"),
             cite="petersen2011", style="strike-slip",
             dist="single-principal", pred="absolute_cm",
             cls="openquake.fdha.primary_surf_displ.petersen2011.Petersen2011PrimaryFD_elliptical",
             eval=lambda: eval_absolute(
                 lambda d, x, m: pet_ell.get_prob(d, x, m, version="elliptical"),
                 mag=REF_MAG_B),
             color="#8c564b", ls="--", marker="s"),
        dict(model_id="petersen2011_quadratic",
             label=make_label("petersen2011", "Quadratic"),
             cite="petersen2011", style="strike-slip",
             dist="single-principal", pred="absolute_cm",
             cls="openquake.fdha.primary_surf_displ.petersen2011.Petersen2011PrimaryFD_quadratic",
             eval=lambda: eval_absolute(
                 lambda d, x, m: pet_quad.get_prob(d, x, m, version="quadratic"),
                 mag=REF_MAG_B),
             color="#e377c2", ls=":", marker="D"),
        dict(model_id="lav2023_strikeslip",
             label=make_label("lavrentiadis2023"),
             cite="lavrentiadis2023", style="strike-slip",
             dist="aggregate-primary", pred="absolute_cm",
             cls="openquake.fdha.primary_surf_displ.lavrentiadis2023.Lavrentiadis2023PrimaryFD",
             eval=lambda: eval_absolute(
                 lambda d, x, m: lv.get_prob(d, x, m, style="strike-slip"),
                 mag=REF_MAG_B),
             color="#1f77b4", ls="-", marker="^"),
        dict(model_id="chiou2025_ss",
             label=make_label("chiou2025"),
             cite="chiou2025", style="strike-slip",
             dist="sum-of-principal", pred="absolute_cm",
             cls="openquake.fdha.primary_surf_displ.chiou2025.Chiou2025PrimaryFD",
             eval=lambda: eval_absolute(
                 lambda d, x, m: ch.get_prob(d, x, m, style="strike-slip"),
                 mag=REF_MAG_B),
             color="#ff7f0e", ls="-.", marker="X"),
        dict(model_id="kuehn2024_strikeslip",
             label=make_label("kuehn2024"),
             cite="kuehn2024", style="strike-slip",
             dist="aggregate", pred="absolute_cm",
             cls="openquake.fdha.primary_surf_displ.kuehn2024.kuehn2024.Kuehn2024PrimaryFD",
             eval=lambda: eval_absolute(
                 lambda d, x, m: ku.get_prob(d, x, m, style="strike-slip",
                                             epistemic_uncertainty=False),
                 mag=REF_MAG_B),
             color="#d62728", ls="--", marker="P"),
    ])
    return panels


# -----------------------------------------------------------------------------
# Documented exclusions -- recorded in the metadata CSV.
# -----------------------------------------------------------------------------
EXCLUSIONS = [
    *[dict(panel_id="c", subfigure="A",
           model_id=f"petersen2011_{variant.lower()}_md",
           display_name=make_label("petersen2011", variant),
           faulting_style="strike-slip", prediction_type="D/MD",
           statistical_distribution="(not source-defined)", percentile=50,
           framework_class_path=("openquake.fdha.primary_surf_displ."
                                 f"petersen2011.Petersen2011PrimaryFD_"
                                 f"{variant.lower()}"),
           exclusion_reason="Petersen et al. (2011) source publication "
                            "publishes ln(D) and ln(D/AD) formulations "
                            "(Eqs. 7-17) but NOT ln(D/MD); no source-"
                            "defined D/MD distribution exists.",
           cite="petersen2011")
      for variant in ["Bilinear", "Elliptical", "Quadratic"]],
]


# -----------------------------------------------------------------------------
# Plot helpers.
# -----------------------------------------------------------------------------
TITLE_STYLE = {"normal": "Normal Faulting",
               "reverse": "Reverse Faulting",
               "strike-slip": "Strike-Slip Faulting"}


def _draw_panel(ax, letter, style_key, curves, y_label, y_lim,
                y_log=False, legend_loc="lower right"):
    seen = set()
    for spec in curves:
        if spec["model_id"] in seen:
            stop(spec["model_id"], ValueError(
                f"duplicate registry key '{spec['model_id']}' in panel"))
        seen.add(spec["model_id"])
        validate_label(spec["model_id"], spec["label"])
        try:
            x, yv = spec["eval"]()
        except Exception as exc:  # noqa: BLE001
            stop(spec["model_id"], exc)
        if not np.all(np.isfinite(yv)):
            stop(spec["model_id"],
                 ValueError("non-finite median values returned"))
        ax.plot(x, yv, color=spec["color"], linestyle=spec["ls"],
                marker=spec["marker"], markevery=2, markersize=7,
                linewidth=2.4, label=spec["label"])
    ax.set_xlim(0.0, 0.5)
    if y_log:
        ax.set_yscale("log")
        ax.set_ylim(*y_lim)
    else:
        ax.set_ylim(*y_lim)
    ax.set_xlabel(r"$x/L$")
    ax.set_ylabel(y_label)
    ax.set_title(f"{TITLE_STYLE[style_key]}")
    ax.text(-0.18, 1.05, f"({letter})", transform=ax.transAxes,
            fontsize=14, fontweight="bold", ha="left", va="bottom")
    ax.grid(True, which="major", color="0.85", linewidth=0.5)
    ax.grid(False, which="minor")
    ax.legend(loc=legend_loc, framealpha=0.9)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    panels_A = build_subfigure_A()
    panels_B = build_subfigure_B()

    # Detect duplicate model_id across the whole figure (STOP condition #2).
    all_ids = {}
    for panel, curves in {**panels_A, **panels_B}.items():
        for spec in curves:
            mid = spec["model_id"]
            if mid in all_ids:
                stop(mid, ValueError(
                    f"duplicate registry key across panels: also in "
                    f"panel {all_ids[mid]}"))
            all_ids[mid] = panel

    fig = plt.figure(figsize=(18, 13))
    gs = fig.add_gridspec(
        3, 3, height_ratios=[1, 1, 1.05],
        hspace=0.55, wspace=0.30,
        left=0.06, right=0.98, top=0.96, bottom=0.06)

    # Subfigure A row 1: D/MD
    for col, (letter, style) in enumerate(
            zip(["a", "b", "c"], ["normal", "reverse", "strike-slip"])):
        ax = fig.add_subplot(gs[0, col])
        _draw_panel(ax, letter, style, panels_A[letter],
                    y_label=r"$D/\mathrm{MD}$", y_lim=(0.0, 1.0))

    # Subfigure A row 2: D/AD
    for col, (letter, style) in enumerate(
            zip(["d", "e", "f"], ["normal", "reverse", "strike-slip"])):
        ax = fig.add_subplot(gs[1, col])
        # Panels (e) and (f) are crowded near the lower-right; place the
        # legend in the upper-left corner where the curves are still below
        # ~D/AD = 0.5 at small x/L.
        loc = "upper left" if letter in ("e", "f") else "lower right"
        _draw_panel(ax, letter, style, panels_A[letter],
                    y_label=r"$D/\mathrm{AD}$", y_lim=(0.0, 2.0),
                    legend_loc=loc)

    # Subfigure B row 3: absolute median displacement (cm) at m = 6.5
    for col, (letter, style) in enumerate(
            zip(["g", "h", "i"], ["normal", "reverse", "strike-slip"])):
        ax = fig.add_subplot(gs[2, col])
        _draw_panel(ax, letter, style, panels_B[letter],
                    y_label="Displacement (cm)", y_lim=(1.0, 100.0),
                    y_log=True, legend_loc="lower right")

    pdf = OUT_DIR / "figure_02_principal_fdm.pdf"
    png = OUT_DIR / "figure_02_principal_fdm.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=300)
    plt.close(fig)
    print(f"wrote {pdf}")
    print(f"wrote {png}")

    # Metadata CSV ------------------------------------------------------------
    csv_path = OUT_DIR / "figure_02_metadata.csv"
    fields = [
        "panel_id", "subfigure", "model_id", "display_name",
        "faulting_style", "prediction_type", "statistical_distribution",
        "percentile", "mw_validity_min", "mw_validity_max",
        "xL_validity_min", "xL_validity_max", "framework_class_path",
        "source_citation", "source_doi", "exclusion_reason",
    ]
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()

        def _mw_range(spec_cite, spec_style=None):
            if spec_cite == "kuehn2024":
                key = {"normal": "kuehn2024_normal",
                       "reverse": "kuehn2024_reverse",
                       "strike-slip": "kuehn2024_strikeslip"}[spec_style]
                return MW_RANGE[key]
            return MW_RANGE.get(spec_cite, ("", ""))

        # Plotted curves
        for sub_label, panels in [("A", panels_A), ("B", panels_B)]:
            for letter, curves in panels.items():
                for spec in curves:
                    cite_str, doi = CITATION[spec["cite"]]
                    mw_lo, mw_hi = _mw_range(spec["cite"], spec["style"])
                    w.writerow(dict(
                        panel_id=letter, subfigure=sub_label,
                        model_id=spec["model_id"],
                        display_name=spec["label"],
                        faulting_style=spec["style"],
                        prediction_type=spec["pred"],
                        statistical_distribution=spec["dist"],
                        percentile=50,
                        mw_validity_min=mw_lo, mw_validity_max=mw_hi,
                        xL_validity_min=XL_RANGE[0],
                        xL_validity_max=XL_RANGE[1],
                        framework_class_path=spec["cls"],
                        source_citation=cite_str, source_doi=doi,
                        exclusion_reason=""))

        # Documented exclusions
        for exc in EXCLUSIONS:
            cite_str, doi = CITATION[exc["cite"]]
            mw_lo, mw_hi = MW_RANGE.get(exc["cite"], ("", ""))
            w.writerow(dict(
                panel_id=exc["panel_id"], subfigure=exc["subfigure"],
                model_id=exc["model_id"],
                display_name=exc["display_name"],
                faulting_style=exc["faulting_style"],
                prediction_type=exc["prediction_type"],
                statistical_distribution=exc["statistical_distribution"],
                percentile=exc["percentile"],
                mw_validity_min=mw_lo, mw_validity_max=mw_hi,
                xL_validity_min=XL_RANGE[0],
                xL_validity_max=XL_RANGE[1],
                framework_class_path=exc["framework_class_path"],
                source_citation=cite_str, source_doi=doi,
                exclusion_reason=exc["exclusion_reason"]))

    print(f"wrote {csv_path}")


if __name__ == "__main__":
    main()
