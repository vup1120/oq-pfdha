"""Schematic of the two rupture-location-weight treatments (manual figure).

Illustrates the two mutually exclusive W_p paths of
``openquake.fdha.calc.location_weight`` -- NOT two settings of one function:

  sigma = 0  distance-threshold treatment: boxcar W_p = 1{|r| <= h} with the
             COMPLEMENTARY distributed weight G = 1 - W_p (Youngs 2003 /
             Takao 2013 either/or);
  sigma > 0  distance-weighted treatment: Petersen's pure Gaussian
             exp(-r^2 / 2 sigma^2), pinned to 1 on the trace and truncated
             at +-2 sigma, with complementary weight G = 1 - W_p.
             h plays NO role on this path.

The four sigma values are the Petersen (2011) Tables 2-3 two-sided
mapping-accuracy classes; their largest ("complex", 0.116 km) is left out of
the curves for legibility and quoted in the caption instead.

Every curve is evaluated through the shipped kernel, so the figure cannot drift
from the implementation.

Run with the worktree on PYTHONPATH:
  PYTHONPATH=<repo> python openquake/fdha/demo/r_sigma_epistemic/wp_schematic.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec

from openquake.fdha.calc.location_weight import location_weight

HERE = Path(__file__).resolve().parent
PNG = HERE / "wp_two_treatments.png"
PDF = HERE / "wp_two_treatments.pdf"

H_KM = 0.1            # r_threshold_km default: the boxcar half-width
R_MAX_M = 300.0
C_BOX = "#1f77b4"     # C0, the repo's sigma = 0 colour

# Petersen (2011) Tables 2-3 two-sided mapping-accuracy classes (km), in the
# repo's Fig.-9c-style greys, dark -> light with worsening accuracy.
CLASSES = [
    ("Accurate",    0.02689, "#1a1a1a"),
    ("Approximate", 0.04382, "#636363"),
    ("Concealed",   0.06552, "#9e9e9e"),
    ("Inferred",    0.07269, "#c4c4c4"),
]
COMPLEX_SIGMA = 0.116  # Table 3, drawn as an extent marker only
MAP_CLASS = 2          # index of the class shown in the map strip (Concealed)

r_m = np.linspace(-R_MAX_M, R_MAX_M, 2001)
r_km = r_m / 1000.0


def wp(sigma_km):
    return location_weight(r_km, r_threshold_km=H_KM, r_sigma_km=sigma_km)


def strip(ax, weight, color, title):
    """Map view: across-strike distance on x (aligned with the panels above),
    along-strike on y, shaded by the principal weight."""
    cmap = LinearSegmentedColormap.from_list("w", ["#ffffff", color])
    ax.imshow(np.tile(weight, (2, 1)), aspect="auto", cmap=cmap,
              vmin=0.0, vmax=1.0, origin="lower",
              extent=[-R_MAX_M, R_MAX_M, 0.0, 1.0], interpolation="bilinear")
    ax.axvline(0.0, color="#b22222", lw=1.6, solid_capstyle="round")
    ax.set_yticks([])
    ax.set_xlim(-R_MAX_M, R_MAX_M)
    ax.set_xlabel("Across-strike distance from mapped trace, $r$ (m)",
                  fontsize=10)
    ax.set_title(title, fontsize=10, loc="left", pad=6)
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)


fig = plt.figure(figsize=(11.5, 6.8), dpi=300)
gs = GridSpec(2, 2, height_ratios=[3.0, 0.8], hspace=0.34, wspace=0.20,
              figure=fig)
ax_wp, ax_g = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
ax_m0 = fig.add_subplot(gs[1, 0], sharex=ax_wp)
ax_m1 = fig.add_subplot(gs[1, 1], sharex=ax_g)

# ---------------------------------------------------------------- (a) W_p
box = wp(0.0)
ax_wp.plot(r_m, box, color=C_BOX, lw=2.4, solid_capstyle="round", zorder=5,
           label=fr"threshold, $h$ = {H_KM*1000:.0f} m  ($\sigma$ = 0)")
for name, s_km, col in CLASSES:
    ax_wp.plot(r_m, wp(s_km), color=col, lw=2.0, solid_capstyle="round",
               label=fr"{name}, $\sigma$ = {s_km*1000:.0f} m")

ax_wp.annotate("hard cutoff\nat $\\pm2\\sigma$\n(54 - 145 m)",
               xy=(150, 0.02), xytext=(186, 0.22),
               fontsize=8.6, color="#444444", ha="left", va="bottom",
               arrowprops=dict(arrowstyle="->", color="#444444", lw=1.0,
                               shrinkA=2, shrinkB=2))
ax_wp.set_ylabel(r"Principal weight  $W_p(r)$", fontsize=10.5)
ax_wp.set_title("(a)  Principal contribution", fontsize=11.5, loc="left")
ax_wp.legend(fontsize=8.4, frameon=False, loc="upper center", ncol=2,
             handlelength=1.7, borderaxespad=0.15, labelspacing=0.35,
             columnspacing=1.4)

# ---------------------------------------------------------------- (b) G
ax_g.plot(r_m, 1.0 - box, color=C_BOX, lw=2.4, solid_capstyle="round",
          zorder=5, label=r"threshold:  $G = 1 - W_p$")
for name, s_km, col in CLASSES:
    ax_g.plot(r_m, 1.0 - wp(s_km), color=col, lw=2.0,
              label=fr"{name}, $\sigma$ = {s_km*1000:.0f} m")
ax_g.set_ylabel(r"Distributed weight  $G(r)$", fontsize=10.5)
ax_g.set_title("(b)  Distributed contribution", fontsize=11.5, loc="left")
ax_g.legend(fontsize=8.4, frameon=False, loc="upper center", ncol=2,
            handlelength=1.7, borderaxespad=0.15, labelspacing=0.35,
            columnspacing=1.4)

for ax in (ax_wp, ax_g):
    ax.set_xlim(-R_MAX_M, R_MAX_M)
    ax.set_ylim(-0.04, 1.42)
    ax.set_yticks([0.0, 0.5, 1.0])
    ax.grid(True, axis="y", color="#d8d8d8", lw=0.8)
    ax.set_axisbelow(True)
    ax.axvline(0.0, color="#b22222", lw=1.0, alpha=0.55, zorder=1)
    ax.tick_params(labelbottom=False)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

# ------------------------------------------------------- (c)/(d) map views
strip(ax_m0, box, C_BOX,
      fr"(c)  Map view: hard-edged band, $|r| \leq$ {H_KM*1000:.0f} m")
name, s_km, col = CLASSES[MAP_CLASS]
strip(ax_m1, wp(s_km), "#1a1a1a",
      fr"(d)  Map view: {name}, $\sigma$ = {s_km*1000:.0f} m, "
      fr"cut at $\pm${2*s_km*1000:.0f} m")

fig.suptitle("Two treatments of rupture-location uncertainty "
             "(not limiting cases of one another)",
             fontsize=12.5, y=0.97)
fig.savefig(PNG, bbox_inches="tight", facecolor="white")
fig.savefig(PDF, bbox_inches="tight", facecolor="white")
print(f"wrote {PNG}\n      {PDF}")
