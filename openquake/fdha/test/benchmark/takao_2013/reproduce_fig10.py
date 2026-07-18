#!/usr/bin/env python
"""Reproduce the Takao et al. (2013) probability-density example (their Fig. 10).

Fig. 10 shows, for the case-(a) worked example (Mw 6.8), the three
probability densities entering their Eq. (12):

    solid - AD:    lognormal of Eq. (10), log10(AD) = -4.80 + 0.69 Mw
    dotted - D/AD:  the conditional gamma distribution
    dashed - D:     the marginal of D = AD x (D/AD)

Empirical identification (see README): the paper's y axis is probability
MASS on a 0.01-decade log grid - the digitized AD peak (x = 0.7793,
h = 0.01108) equals the Eq. 10 lognormal mass exactly (theory x = 0.7798,
h = 0.01108). The dotted/dashed curves correspond to the SHORT-segment
gamma of their Eq. (8) (a = 1.53, b = 0.58, the L < 10 km branch), not the
long-segment Eq. (6) at x/L = 0 that the figure caption suggests: the
digitized D/AD peak sits at 0.894 (Eq. 8 log-mass: 0.887; Eq. 6 at
x/L = 0: 0.497) and the digitized D peak at (0.584, 7.57e-3) matches the
Eq. 8 convolution (0.586, 7.57e-3) with no free scale, while Eq. 6 gives
(0.336, 8.17e-3). The dotted curve alone carries an unexplained uniform
scale factor (~1/7.9, fitted and reported, shape compared).

The D curve is the strongest single test: it validates the AD-grid mass x
gamma composition exactly as ``Takao2013PrimaryFD`` implements it,
against the paper's own plotted numbers, with no adjustable parameter.

Writes ``Figures/takao2013_fig10_reproduction.png`` and
``fig10_agreement.json``.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from openquake.fdha.primary_surf_displ import Takao2013PrimaryFD

HERE = Path(__file__).resolve().parent

MAG = 6.8
GRID_DLOG10 = 0.01  # the paper's plotting grid (empirical, see docstring)
N_SIGMA = 5.0
SHORT_SRL_KM = 1.0  # any value < 10 selects the Eq. 8 coefficients

# digitized-peak anchors measured from reference/fig10_points.csv
PAPER_AD_PEAK = (0.7793, 0.01108)
PAPER_D_PEAK = (0.5836, 0.007573)
PAPER_DAD_PEAK_X = 0.8942


def theoretical_curves():
    """(x, mass) log10-grid mass curves for AD, D/AD (Eq. 8) and D."""
    model = Takao2013PrimaryFD(n_sigma=N_SIGMA)
    log_mean = -4.80 + 0.69 * MAG
    sigma = 0.36

    lz = np.arange(-3, 3 + 1e-9, GRID_DLOG10)
    z = 10 ** lz

    # AD: the model's normalized mass on its own integration grid,
    # rescaled to mass per GRID_DLOG10
    model_dlog = 2 * N_SIGMA * sigma / 999  # get_prob_* use 1000-point grids
    ad_mass = np.array([model.get_prob_avg_displacement(v, MAG) for v in z])
    ad_mass *= GRID_DLOG10 / model_dlog

    # D/AD: mass from the Eq. 8 survival function (short-segment branch)
    sf = np.array([float(model.get_prob_D_AD(v, 0.0, SHORT_SRL_KM)) for v in z])
    dad_mass = -np.gradient(sf, lz)* GRID_DLOG10

    # D = AD x (D/AD): exceedance composed exactly as the model composes it
    # (AD-grid masses x gamma survival), differentiated on the log grid
    lower = 10 ** (log_mean - N_SIGMA * sigma)
    upper = 10 ** (log_mean + N_SIGMA * sigma)
    ad_grid = np.logspace(np.log10(lower), np.log10(upper), 1000)
    ad_pdf = np.array([model.get_prob_avg_displacement(v, MAG) for v in ad_grid])
    exceed = np.zeros_like(z)
    for ad, p in zip(ad_grid, ad_pdf):
        exceed += p * np.asarray(
            [float(model.get_prob_D_AD(v / ad, 0.0, SHORT_SRL_KM)) for v in z])
    d_mass = -np.gradient(exceed, lz) * GRID_DLOG10

    return z, ad_mass, dad_mass, d_mass


def digitized_curves():
    """Split the extracted Fig 10 points into (AD, D, D/AD) point sets."""
    import csv

    rows = list(csv.DictReader(
        (HERE / "reference" / "fig10_points.csv").open()))
    x = np.array([float(r["x"]) for r in rows])
    y = np.array([float(r["y"]) for r in rows])
    src = np.array([r["source"] for r in rows])
    sz = np.array([float(r["size_pt"]) for r in rows])

    ad = np.c_[x[src == "stroke"], y[src == "stroke"]]

    # D backbone: the large-dash elements
    m = (src == "dash") & (sz >= 0.85) & (y > 1e-4)
    order = np.argsort(x[m])
    d_pts = np.c_[x[m][order], y[m][order]]

    # dotted D/AD: remaining dash elements clearly below the D backbone
    m_all = (src == "dash") & (y > 2e-4)
    y_on_d = np.interp(np.log10(x[m_all]), np.log10(d_pts[:, 0]), d_pts[:, 1])
    keep = y[m_all] < 0.55 * y_on_d
    dad_pts = np.c_[x[m_all][keep], y[m_all][keep]]
    return ad, d_pts, dad_pts


def _peak(xs, ys):
    i = int(np.argmax(ys))
    s = slice(max(0, i - 6), i + 7)
    c = np.polyfit(np.log10(xs[s]), ys[s], 2)
    lx = -c[1] / (2 * c[0])
    return float(10 ** lx), float(np.polyval(c, lx))


def _pointwise(pts, xg, yg, threshold=0.05):
    """max |relerr| of digitized points vs curve, where curve > threshold*peak."""
    yi = np.interp(np.log10(pts[:, 0]), np.log10(xg), yg)
    m = yi > threshold * yg.max()
    if not m.any():
        return None
    return float(np.max(np.abs(pts[m, 1] / yi[m] - 1)))


def agreement():
    z, ad_mass, dad_mass, d_mass = theoretical_curves()
    ad_pts, d_pts, dad_pts = digitized_curves()

    res = {}
    for name, pts, curve in [("AD", ad_pts, ad_mass), ("D", d_pts, d_mass)]:
        px, ph = _peak(pts[:, 0], pts[:, 1])
        tx, th = _peak(z, curve)
        res[name] = {
            "digitized_peak": [px, ph], "computed_peak": [tx, th],
            "peak_x_ratio": px / tx, "peak_h_ratio": ph / th,
            "max_relerr_above_5pct_peak": _pointwise(pts, z, curve),
        }

    # dotted D/AD: shape comparison with one fitted scale (see docstring)
    px, ph = _peak(dad_pts[:, 0], dad_pts[:, 1])
    tx, th = _peak(z, dad_mass)
    yi = np.interp(np.log10(dad_pts[:, 0]), np.log10(z), dad_mass)
    m = yi > 0.1 * dad_mass.max()
    scale = float(np.sum(dad_pts[m, 1] * yi[m]) / np.sum(yi[m] ** 2))
    res["D_AD"] = {
        "digitized_peak_x": px, "computed_peak_x": tx, "peak_x_ratio": px / tx,
        "fitted_scale": scale,
        "max_shape_relerr_above_10pct_peak": float(
            np.max(np.abs(dad_pts[m, 1] / (scale * yi[m]) - 1))),
    }
    return res, (z, ad_mass, dad_mass, d_mass), (ad_pts, d_pts, dad_pts)


def make_figure(curves, points, res):
    z, ad_mass, dad_mass, d_mass = curves
    ad_pts, d_pts, dad_pts = points
    scale = res["D_AD"]["fitted_scale"]

    fig, ax = plt.subplots(figsize=(8.5, 5))
    ax.plot(z, ad_mass, "-", color="k", lw=1.5, label="AD (Eq. 10 lognormal)")
    ax.plot(z, d_mass, "--", color="k", lw=1.5, label="D (Eq. 12 convolution)")
    ax.plot(z, dad_mass * scale, ":", color="k", lw=1.5,
            label=f"D/AD (Eq. 8) x fitted {scale:.3g}")
    ax.plot(ad_pts[:, 0], ad_pts[:, 1], ".", ms=2.5, color="tab:blue",
            label="paper solid (digitized)")
    ax.plot(d_pts[:, 0], d_pts[:, 1], ".", ms=2.5, color="tab:red",
            label="paper dashed (digitized)")
    ax.plot(dad_pts[:, 0], dad_pts[:, 1], ".", ms=2.5, color="tab:green",
            label="paper dotted (digitized)")
    ax.set_xscale("log")
    ax.set_xlim(1e-3, 1e3)
    ax.set_ylim(0, 1.2e-2)
    ax.set_xlabel("AD (m), D/AD (unitless), D (m)")
    ax.set_ylabel(f"Probability mass per {GRID_DLOG10} decade")
    ax.set_title("Takao et al. (2013) Fig. 10 - computed vs digitized")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    figdir = HERE / "Figures"
    figdir.mkdir(exist_ok=True)
    out = figdir / "takao2013_fig10_reproduction.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main():
    res, curves, points = agreement()
    out_png = make_figure(curves, points, res)
    (HERE / "fig10_agreement.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))
    print(f"figure: {out_png}")


if __name__ == "__main__":
    import os

    os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
    main()
