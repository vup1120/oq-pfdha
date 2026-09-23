"""Petersen-(2011)-Fig.-9-style displacement profile across the fault.

x axis: fault distance (m, transect perpendicular to the demo M7 trace);
y axis: displacement with a 10,000-yr return period (lambda* = 1e-4/yr).

Uses the same job chain as run_demo.py (Youngs2003 models, boxcar
h = r_threshold_km = 0.05 km) but with the single M7.0 characteristic event
at a Petersen-style rate of 1/140 yr (the demo's 1e-4/yr event can never
reach a 1e-4/yr exceedance rate, so the RP-10,000 profile would be zero).

Curves, exercising the two W_p paths:

  sigma = 0      boxcar + COMPLEMENTARY split -> box profile: full principal
                 inside +-50 m, distributed-only shoulders outside (blue);
  four Gaussians the Petersen (2011) Table 2-3 two-sided mapping-accuracy
                 classes - Accurate 26.89 m, Approximate 43.82 m, Concealed
                 65.52 m, Inferred 72.69 m - pure Gaussian W_p SUMMED with
                 distributed, in Fig.-9c-style greys (dark -> light).

Run with the worktree on PYTHONPATH:
  PYTHONPATH=<repo> python openquake/fdha/demo/r_sigma_epistemic/profile_figure.py
"""
from __future__ import annotations

import shutil
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

from openquake.fdha.logic_tree.driver import FdhaLogicTree

HERE = Path(__file__).resolve().parent
BASE_DEMO = HERE.parents[3] / "examples"
LOCAL_SRC = HERE / "source_model_M7.xml"
SRC_XML = "source_model.xml"
SMLT_XML = "hazard_curve_minimal_source_model_logic_tree.xml"
FDHA_XML = "hazard_curve_minimal_fdha_logic_tree.xml"
JOB_INI = "hazard_curve_minimal.ini"
OUT = HERE / "out_profile"

RP_YEARS = 10_000.0
LAM_T = 1.0 / RP_YEARS                 # 1e-4 /yr
EVENT_RATE = 1.0 / 140.0               # Petersen's M7 example recurrence
H_BOXCAR = 0.05                        # km
S_ZERO = "0.0"
# Petersen (2011) Tables 2-3 two-sided mapping-accuracy classes (km), with
# Fig.-9c-style greys, dark -> light with worsening accuracy.
CLASSES = [
    ("Accurate",    "0.02689", "#1a1a1a"),
    ("Approximate", "0.04382", "#636363"),
    ("Concealed",   "0.06552", "#9e9e9e"),
    ("Inferred",    "0.07269", "#cfcfcf"),
]
R_MAX_M = 300.0                        # transect half-extent (m)
# Reference point of the original example (used only to anchor the transect).
S_REF = np.array([16.16573727, 39.64704451])
# Dense displacement grid for a smooth return-period inversion.
D_LEVELS = [round(float(v), 6) for v in np.logspace(-3, np.log10(20.0), 40)]


def read_trace():
    ns = {"gml": "http://www.opengis.net/gml"}
    root = ET.parse(LOCAL_SRC).getroot()
    pos = root.find(".//gml:posList", ns).text.split()
    return np.array([float(v) for v in pos]).reshape(-1, 2)


def transect_sites():
    """Sites on the perpendicular through the trace point nearest S_REF,
    at signed fault distances -R_MAX..+R_MAX (m). Returns (r_m, lonlats)."""
    trace = read_trace()
    kx, ky = 111.0 * np.cos(np.radians(S_REF[1])), 111.0

    def to_km(p):
        return np.array([(p[0] - S_REF[0]) * kx, (p[1] - S_REF[1]) * ky])

    s = np.array([0.0, 0.0])           # S_REF in km space
    best = (np.inf, None)
    for i in range(len(trace) - 1):
        a, b = to_km(trace[i]), to_km(trace[i + 1])
        ab = b - a
        t = np.clip(np.dot(s - a, ab) / np.dot(ab, ab), 0.0, 1.0)
        f = a + t * ab
        d = np.hypot(*(s - f))
        if d < best[0]:
            best = (d, f)
    d0, foot = best
    u = (s - foot) / d0                # unit vector, trace -> S_REF side
    # Composite sampling: 2 m inside +-160 m (so the boxcar cliff at +-h and
    # every class's Gaussian toe at +-2 sigma, up to 145 m for Inferred, are
    # resolved - with 10 m steps a vertical step is drawn as a 10 m ramp),
    # 20 m in the flat far field.
    fine = np.arange(-160.0, 160.0 + 0.1, 2.0)
    coarse = np.concatenate([
        np.arange(-R_MAX_M, -160.0, 20.0), np.arange(180.0, R_MAX_M + 0.1, 20.0)
    ])
    r_m = np.unique(np.concatenate([fine, coarse]))
    lonlats = []
    for r in r_m:
        p = foot + (r / 1000.0) * u
        lonlats.append((S_REF[0] + p[0] / kx, S_REF[1] + p[1] / ky))
    return r_m, lonlats


def sigma_level(branches):
    rows = [
        '    <logicTreeBranchSet branchSetID="bs_5_r_sigma" '
        'uncertaintyType="fdhaCalcRSigma">',
    ]
    for bid, value, weight in branches:
        rows += [
            f'      <logicTreeBranch branchID="{bid}">',
            f'        <uncertaintyModel>{value}</uncertaintyModel>',
            f'        <uncertaintyWeight>{weight}</uncertaintyWeight>',
            '      </logicTreeBranch>',
        ]
    rows += ['    </logicTreeBranchSet>']
    return "\n".join(rows)


def make_job(name, branches, sites_str):
    workdir = OUT / name
    workdir.mkdir(parents=True, exist_ok=True)
    src = LOCAL_SRC.read_text().replace(
        "<occurRates>1.0000000E-04</occurRates>",
        f"<occurRates>{EVENT_RATE:.7E}</occurRates>",
    )
    (workdir / SRC_XML).write_text(src)
    shutil.copy(BASE_DEMO / SMLT_XML, workdir)

    lt_text = (BASE_DEMO / FDHA_XML).read_text()
    lt_text = lt_text.replace(
        "  </logicTree>", sigma_level(branches) + "\n  </logicTree>"
    )
    (workdir / FDHA_XML).write_text(lt_text)

    ini_text = (BASE_DEMO / JOB_INI).read_text()
    ini_text = ini_text.replace(
        "r_threshold_km = 0.1", f"r_threshold_km = {H_BOXCAR}"
    )
    ini_text = ini_text.replace(
        "sites = 16.16573727 39.64704451", f"sites = {sites_str}"
    )
    import re
    ini_text = re.sub(
        r"displacement_measure_levels = .*",
        'displacement_measure_levels = {"FD": %s}' % D_LEVELS,
        ini_text,
    )
    ini = workdir / "job.ini"
    ini.write_text(ini_text)
    return ini


def run_rates(name, branches, sites_str):
    ini = make_job(name, branches, sites_str)
    res = FdhaLogicTree.from_ini(str(ini)).run(outdir=OUT / name / "out")
    return np.asarray(res.d0, float), np.asarray(res.mean_rates, float)


def displacement_at_rp(d0, rates):
    """Log-log inversion of each site's hazard curve at LAM_T (0 below)."""
    out = np.zeros(rates.shape[0])
    for i, row in enumerate(rates):
        ok = row > 0
        if not ok.any() or row[ok].max() < LAM_T:
            continue
        x = np.log(row[ok][::-1])           # increasing
        y = np.log(d0[ok][::-1])
        out[i] = float(np.exp(np.interp(np.log(LAM_T), x, y)))
    return out


def main():
    if OUT.exists():
        shutil.rmtree(OUT)
    r_m, lonlats = transect_sites()
    sites_str = ", ".join(f"{lo:.8f} {la:.8f}" for lo, la in lonlats)
    print(f"{len(r_m)} transect sites, r = {r_m[0]:.0f}..{r_m[-1]:.0f} m")

    print("running sigma = 0 branch (boxcar + complementary)...")
    d0, rates_box = run_rates("s_zero", [("RS_ZERO", S_ZERO, "1.0")], sites_str)
    prof_box = displacement_at_rp(d0, rates_box)
    print(f"  sigma=0 peak {prof_box.max():.2f} m")

    profs = {}
    for name, sig, _col in CLASSES:
        print(f"running {name} (sigma = {sig} km, Gaussian + sum)...")
        _, rates = run_rates(f"s_{name.lower()}",
                             [(f"RS_{name.upper()}", sig, "1.0")], sites_str)
        profs[name] = displacement_at_rp(d0, rates)
        print(f"  {name} peak {profs[name].max():.2f} m")

    # ---------------------------------------------------------------- figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9.5, 5.6), dpi=150)
    for name, sig, col in CLASSES:
        ax.plot(r_m, profs[name], color=col, lw=2.2, solid_capstyle="round",
                label=f"{name} ($\\sigma$ = {float(sig)*1000:.2f} m)")
    ax.plot(r_m, prof_box, "C0-", lw=2.4, solid_capstyle="round",
            label=f"$\\sigma$ = 0 (boxcar h = {int(H_BOXCAR*1000)} m)")

    for x in (-H_BOXCAR * 1000, H_BOXCAR * 1000):
        ax.axvline(x, color="C0", lw=0.9, ls=(0, (4, 3)), alpha=0.6)
    ax.annotate("$\\pm h$ = 50 m",
                xy=(-H_BOXCAR * 1000, prof_box.max() * 0.55),
                xytext=(-H_BOXCAR * 1000 - 105, prof_box.max() * 0.62),
                color="C0", fontsize=10,
                arrowprops=dict(arrowstyle="->", color="C0", lw=1.0))

    ax.set_xlabel("Fault distance (m)", fontsize=11)
    ax.set_ylabel(f"Displacement (m), RP = {int(RP_YEARS):,} yr", fontsize=11)
    ax.set_xlim(-R_MAX_M, R_MAX_M)
    ax.set_ylim(bottom=0)
    ax.grid(True, axis="y", color="#d8d8d8", lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.text(
        0.02, 0.97,
        f"RP = {int(RP_YEARS):,} yr ($\\lambda^*$ = {LAM_T:.0e}/yr)\n"
        f"single M 7.0, rate 1/140 yr\nYoungs2003 principal + distributed",
        transform=ax.transAxes, va="top", fontsize=9.5,
        bbox=dict(fc="white", ec="0.6", boxstyle="round,pad=0.4"),
    )
    ax.legend(fontsize=8.5, frameon=False, loc="upper right",
              bbox_to_anchor=(1.0, 1.0), handlelength=1.3,
              labelspacing=0.4, borderaxespad=0.0)
    ax.set_title("Displacement across the fault", fontsize=12)
    fig.tight_layout()
    png = OUT / "r_sigma_profile_RP10000.png"
    fig.savefig(png, bbox_inches="tight", facecolor="white")
    print(f"figure written to {png}")


if __name__ == "__main__":
    main()
