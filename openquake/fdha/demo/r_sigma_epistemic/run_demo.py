"""Demo: r_sigma_km (mapping accuracy) as epistemic logic-tree branches.

Derives job variants from the canonical ``examples/hazard_curve_minimal``
job (Youngs2003 normal models), runs them at TWO sites, checks the
MODE A/MODE B equivalence, and renders a four-panel figure:

  A. map of the fault trace, the two sites, the 50 m boxcar band and the
     Gaussian W_p bands (drawn for the widest class, Inferred);
  B. hazard curves at site A (40 m from the trace);
  C. hazard curves at site B (80 m from the trace);
  D. machine-precision consistency check (both sites).

The rupture-location weight W_p has two separate paths, each with its own
combination rule (docs/design/rupture_location_uncertainty.md, section 2):

  sigma = 0    -> boxcar 1{|r| <= r_threshold_km}, COMPLEMENTARY split
                  (inside h only principal, outside only distributed);
  sigma > 0    -> Petersen's pure Gaussian exp(-r^2/2 sigma^2), pinned,
                  truncated at +-2 sigma, SUMMED with the full distributed
                  term (Petersen eq. 1 + eq. 2).

The sigma > 0 branches are the four Petersen (2011) Table 2-3 two-sided
mapping-accuracy classes - Accurate 26.89 m, Approximate 43.82 m, Concealed
65.52 m, Inferred 72.69 m - drawn in Fig.-9c-style greys (dark -> light).

The source is a demo-local **single M7.0 characteristic event** at a low
annual rate (1e-4/yr; see ``source_model_M7.xml``) so the curves plateau
near ~1e-4/yr. The two sites bracket the 50 m boxcar edge:

  site A, r = 0.04 km (inside the boxcar):
      sigma=0: full principal, no distributed (W_p = 1, either/or);
      classes: W_p = 0.33 / 0.66 / 0.83 / 0.86 x principal + distributed.
  site B, r = 0.08 km (outside the boxcar):
      sigma=0: distributed only (W_p = 0);
      classes: W_p = 0 (Accurate, beyond its 2 sigma) / 0.19 / 0.48 / 0.55
      x principal + distributed.

PROPAGATION (the point of the epistemic mechanism): the ``weighted`` run
carries ALL FOUR classes in ONE ``fdhaCalcRSigma`` branch set (weight 0.25
each). The driver enumerates the branches, runs each end branch, and
aggregates: the weighted-mean curve plus the canonical fractiles
(5/16/50/84/95%). Because the mean is linear in the branch rates, the
aggregated mean must equal 0.25 x the sum of the four single-class runs at
machine precision (verification V7 of the design doc) - asserted below.
A second figure (``out/r_sigma_epistemic_propagation.png``) shows the four
branch curves, the weighted mean and the 5-95% fractile band per site.

All values are demo inputs, not recommendations.
"""
from __future__ import annotations

import json
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

from openquake.fdha.logic_tree.driver import FdhaLogicTree

HERE = Path(__file__).resolve().parent
BASE_DEMO = HERE.parents[3] / "examples"
LOCAL_SRC = HERE / "source_model_M7.xml"       # single M7.0 scenario, 1e-4/yr
SRC_XML = "source_model.xml"                    # name the SMLT expects
SMLT_XML = "hazard_curve_minimal_source_model_logic_tree.xml"
FDHA_XML = "hazard_curve_minimal_fdha_logic_tree.xml"
JOB_INI = "hazard_curve_minimal.ini"
OUT = HERE / "out"

# Two demo sites bracketing the 50 m boxcar edge (0.04 / 0.08 km from the
# trace) and the boxcar half-width written into the INI.
SITE_A = (16.16186114, 39.64824366, 0.04)      # lon, lat, r (km) - inside
SITE_B = (16.16229540, 39.64810931, 0.08)      # lon, lat, r (km) - outside
H_BOXCAR = 0.05                                 # r_threshold_km in the INI
S_ZERO = "0.0"
# Petersen (2011) Tables 2-3 two-sided mapping-accuracy classes (km), with
# Fig.-9c-style greys, dark -> light with worsening accuracy.
CLASSES = [
    ("Accurate",    "0.02689", "#1a1a1a"),
    ("Approximate", "0.04382", "#636363"),
    ("Concealed",   "0.06552", "#9e9e9e"),
    ("Inferred",    "0.07269", "#cfcfcf"),
]
COL_A, COL_B = "#2166ac", "#2e8b57"             # site markers / note colours


def wp_gauss(r_km: float, sig_km: float) -> float:
    if r_km > 2.0 * sig_km:
        return 0.0
    return float(np.exp(-r_km**2 / (2.0 * sig_km**2)))


def sigma_level(branches):
    rows = [
        '    <logicTreeBranchingLevel branchingLevelID="bl_5_r_sigma">',
        '      <logicTreeBranchSet branchSetID="bs_5_r_sigma" '
        'uncertaintyType="fdhaCalcRSigma">',
    ]
    for bid, value, weight in branches:
        rows += [
            f'        <logicTreeBranch branchID="{bid}">',
            f'          <uncertaintyModel>{value}</uncertaintyModel>',
            f'          <uncertaintyWeight>{weight}</uncertaintyWeight>',
            '        </logicTreeBranch>',
        ]
    rows += ['      </logicTreeBranchSet>', '    </logicTreeBranchingLevel>']
    return "\n".join(rows)


def make_job(name: str, branches=None) -> Path:
    """Clone the examples/hazard_curve_minimal job into out/<name>/, with the
    demo-local single-M7.0 source swapped in and the two demo sites.

    baseline (branches=None): MODE A, no r_sigma_km anywhere == sigma 0
    (boxcar path), the tool default.

    MODE B cases: the ``fdhaCalcRSigma`` branch set is appended to the FDHA
    logic tree. The INI never carries an ``r_sigma_km`` scalar here (the two
    mechanisms are mutually exclusive). The boxcar half-width
    ``r_threshold_km`` stays in the INI - a fixed calculation parameter that
    only matters on the sigma = 0 branch.
    """
    workdir = OUT / name
    workdir.mkdir(parents=True, exist_ok=True)
    shutil.copy(LOCAL_SRC, workdir / SRC_XML)
    shutil.copy(BASE_DEMO / SMLT_XML, workdir)

    lt_text = (BASE_DEMO / FDHA_XML).read_text()
    if branches:
        lt_text = lt_text.replace(
            "  </logicTree>", sigma_level(branches) + "\n  </logicTree>"
        )
    (workdir / FDHA_XML).write_text(lt_text)
    ini_text = (BASE_DEMO / JOB_INI).read_text()
    ini_text = ini_text.replace(
        "r_threshold_km = 0.1", f"r_threshold_km = {H_BOXCAR}"
    )
    ini_text = ini_text.replace(
        "sites = 16.16573727 39.64704451",
        f"sites = {SITE_A[0]} {SITE_A[1]}, {SITE_B[0]} {SITE_B[1]}",
    )
    (workdir / JOB_INI).write_text(ini_text)
    return workdir / JOB_INI


def run(name: str, branches=None):
    ini = make_job(name, branches)
    res = FdhaLogicTree.from_ini(str(ini)).run(outdir=OUT / name / "out")
    d0 = np.asarray(res.d0, dtype=float)
    mean = np.asarray(res.mean_rates, dtype=float)          # (n_sites, D)
    fractiles = {
        float(q): np.asarray(v, dtype=float)                # (n_sites, D)
        for q, v in (res.fractiles or {}).items()
    }
    manifest = json.loads((OUT / name / "out" / "manifest.json").read_text())
    return {"d0": d0, "mean": mean, "fractiles": fractiles,
            "manifest": manifest}


def rel_diff(a, b):
    return np.abs(a - b) / np.where(b != 0, np.abs(b), 1.0)


def read_trace():
    """Return the fault-trace lon/lat vertices from the demo source model."""
    ns = {"gml": "http://www.opengis.net/gml"}
    root = ET.parse(LOCAL_SRC).getroot()
    pos = root.find(".//gml:posList", ns).text.split()
    vals = [float(v) for v in pos]
    return np.array(vals).reshape(-1, 2)  # (n, 2) lon, lat


def dist_to_polyline_km(lon_grid, lat_grid, trace, ref_lon, ref_lat):
    """Min distance (km) from each grid point to the trace polyline."""
    kx = 111.0 * np.cos(np.radians(ref_lat))
    ky = 111.0
    gx = (lon_grid - ref_lon) * kx
    gy = (lat_grid - ref_lat) * ky
    tx = (trace[:, 0] - ref_lon) * kx
    ty = (trace[:, 1] - ref_lat) * ky
    dmin = np.full(gx.shape, np.inf)
    for i in range(len(trace) - 1):
        ax, ay = tx[i], ty[i]
        bx, by = tx[i + 1], ty[i + 1]
        abx, aby = bx - ax, by - ay
        denom = abx * abx + aby * aby
        t = np.clip(((gx - ax) * abx + (gy - ay) * aby) / denom, 0.0, 1.0)
        px, py = ax + t * abx, ay + t * aby
        dmin = np.minimum(dmin, np.hypot(gx - px, gy - py))
    return dmin


def main():
    if OUT.exists():
        shutil.rmtree(OUT)

    print("Running baseline (MODE A, single M7.0 source, no r_sigma_km == sigma 0)...")
    baseline = run("baseline")
    print("Running case1 (MODE B: r_sigma_km = 0.0, w=1)...")
    case1 = run("case1", [("RS_ZERO", S_ZERO, "1.0")])
    class_res = {}
    for name, sig, _col in CLASSES:
        print(f"Running {name} (MODE B: r_sigma_km = {sig}, w=1)...")
        class_res[name] = run(
            name.lower(), [(f"RS_{name.upper()}", sig, "1.0")])

    # PROPAGATION: all four classes in ONE weighted branch set. The driver
    # enumerates the end branches, runs each one, and aggregates the weighted
    # mean and the canonical fractiles across them.
    print("Running weighted (MODE B: four classes, w=0.25 each, ONE tree)...")
    weighted = run("weighted", [
        (f"RS_{name.upper()}", sig, "0.25") for name, sig, _col in CLASSES
    ])

    d0 = baseline["d0"]
    e_equiv = rel_diff(case1["mean"], baseline["mean"]).max()
    print()
    print(f"case1 (branch sigma=0, w=1) vs baseline: max rel diff = {e_equiv:.3e}")
    for i, label in enumerate(("A", "B")):
        row = " ".join(
            f"{n[:3]}={class_res[n]['mean'][i][0]:.2e}" for n, _, _ in CLASSES)
        print(f"site {label}: plateau sigma0={case1['mean'][i][0]:.2e}  {row}")

    assert e_equiv == 0.0, "MODE B (sigma=0) must be bit-identical to MODE A default"
    # The widest class must genuinely differ from the boxcar at both sites.
    for i in range(2):
        spread = rel_diff(
            class_res["Inferred"]["mean"][i], case1["mean"][i]).max()
        assert spread > 0.1, f"Inferred must differ from sigma=0 at site {i}"

    # V7 linearity: the aggregated mean of the weighted tree must equal the
    # weighted sum of the four independent single-class runs.
    lin_expected = 0.25 * sum(
        class_res[name]["mean"] for name, _sig, _col in CLASSES)
    e_lin = rel_diff(weighted["mean"], lin_expected).max()
    print(f"weighted tree mean vs 0.25*sum(single runs): "
          f"max rel diff = {e_lin:.3e}")
    assert e_lin < 1e-12, "aggregated mean must be linear in the branch rates"

    # The manifest documents the propagated branches: composed branch ids,
    # the sigma value each realization ran with, and the combined weights.
    print("weighted tree realizations (from out/weighted/out/manifest.json):")
    for b in weighted["manifest"]["branches"]:
        sig_val = b["fdha_calc_params"]["r_sigma_km"]
        print(f"  {b['fdha_branch_id'].split('|')[0]:<15} "
              f"r_sigma_km={sig_val:<8} w={b['combined_branch_weight']}")
    w_sum = sum(b["combined_branch_weight"]
                for b in weighted["manifest"]["branches"])
    assert abs(w_sum - 1.0) < 1e-12, "combined weights must sum to 1"

    # ------------------------------------------------------------- figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    fig, (axm, axA, axB) = plt.subplots(
        1, 3, figsize=(17, 5.4), gridspec_kw={"width_ratios": [1.9, 2.5, 2.5]}
    )

    # -- Panel A: map of trace, sites, boxcar band + Gaussian W_p bands ----
    # Bands drawn for the WIDEST class (Inferred): its 2 sigma support is the
    # envelope of all four classes.
    trace = read_trace()
    c_lon = 0.5 * (SITE_A[0] + SITE_B[0])
    c_lat = 0.5 * (SITE_A[1] + SITE_B[1])
    dlon, dlat = 0.014, 0.011
    lon_min, lon_max = c_lon - dlon, c_lon + dlon
    lat_min, lat_max = c_lat - dlat, c_lat + dlat
    LON, LAT = np.meshgrid(
        np.linspace(lon_min, lon_max, 500), np.linspace(lat_min, lat_max, 500)
    )
    dist = dist_to_polyline_km(LON, LAT, trace, c_lon, c_lat)
    sig_env = float(CLASSES[-1][1])                # Inferred, widest
    r_half = sig_env * np.sqrt(2.0 * np.log(2.0))  # W_p = 0.5 radius
    axm.contourf(
        LON, LAT, dist, levels=[0, H_BOXCAR, r_half, 2 * sig_env, 1e9],
        colors=["#b2182b", "#ef8a62", "#fddbc7", "#e8e8e8"], alpha=0.8,
    )
    axm.plot(trace[:, 0], trace[:, 1], "k-", lw=2.4, marker="o", ms=4)
    axm.plot(SITE_A[0], SITE_A[1], marker="*", ms=8, mfc=COL_A,
             mec="white", mew=0.5, ls="none", zorder=6)
    axm.plot(SITE_B[0], SITE_B[1], marker="*", ms=8, mfc=COL_B,
             mec="white", mew=0.5, ls="none", zorder=6)
    axm.text(0.04, 0.32, f"site A: r = {int(SITE_A[2]*1000)} m\n(inside boxcar)",
             transform=axm.transAxes, fontsize=9.5, color=COL_A, va="top")
    axm.text(0.96, 0.04, f"site B: r = {int(SITE_B[2]*1000)} m\n(outside boxcar)",
             transform=axm.transAxes, fontsize=9.5, color=COL_B,
             ha="right", va="bottom")
    axm.set_xlim(lon_min, lon_max)
    axm.set_ylim(lat_min, lat_max)
    axm.set_aspect(1.0 / np.cos(np.radians(c_lat)))
    axm.set_xlabel("Longitude"); axm.set_ylabel("Latitude")
    axm.set_title("Two sites bracket the 50 m boxcar edge")
    axm.legend(
        handles=[
            Line2D([], [], color="k", lw=2.4, marker="o", ms=4, label="fault trace"),
            Patch(fc="#b2182b", alpha=0.8,
                  label=f"$\\leq${int(H_BOXCAR*1000)} m: boxcar ($\\sigma$=0)"),
            Patch(fc="#ef8a62", alpha=0.8,
                  label="$W_p\\geq0.5$ (Inferred $\\sigma$=72.69 m)"),
            Patch(fc="#fddbc7", alpha=0.8,
                  label=f"$W_p>0$ to $2\\sigma$={int(2*sig_env*1000)} m"),
        ],
        fontsize=7.5, loc="upper left",
    )

    # -- Panels B/C: hazard curves per site ---------------------------------
    for i, (ax, label, site) in enumerate((
        (axA, "A", SITE_A), (axB, "B", SITE_B),
    )):
        r = site[2]
        box = 1.0 if r <= H_BOXCAR else 0.0
        ax.loglog(d0, case1["mean"][i], color="C0", lw=3.6,
                  solid_capstyle="round", zorder=2,
                  label=f"$\\sigma$ = 0 (boxcar, $W_p$={box:.0f})")
        for name, sig, col in CLASSES:
            ax.loglog(d0, class_res[name]["mean"][i], color=col, lw=2.0,
                      label=f"{name} ($\\sigma$={float(sig)*1000:.2f} m, "
                            f"$W_p$≈{wp_gauss(r, float(sig)):.2f})")
        ax.set_xlabel("Displacement $D_0$ [m]")
        ax.set_title(
            f"Site {label} at {int(r*1000)} m "
            f"({'inside' if box else 'outside'} the boxcar)", fontsize=10.5)
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(fontsize=8, loc="lower left")
    axA.set_ylabel("Annual rate of exceedance")

    fig.tight_layout()
    png = OUT / "r_sigma_epistemic_demo.png"
    fig.savefig(png, dpi=150)
    print(f"\nFigure written to {png}")

    # -- Figure 2: propagation through the weighted tree --------------------
    # Per site: the four branch curves (the enumerated end branches of the
    # ONE weighted tree - identical to the single-class runs by V7), the
    # aggregated weighted mean, and the 5-95% fractile band.
    q_lo, q_hi = 0.05, 0.95
    fig2, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), sharey=True)
    for i, (ax, label, site) in enumerate((
        (axes[0], "A", SITE_A), (axes[1], "B", SITE_B),
    )):
        for name, sig, col in CLASSES:
            ax.loglog(d0, class_res[name]["mean"][i], color=col, lw=1.4,
                      zorder=2, label=f"branch {name} "
                      f"($\\sigma$={float(sig)*1000:.2f} m, w=0.25)")
        if q_lo in weighted["fractiles"] and q_hi in weighted["fractiles"]:
            ax.fill_between(
                d0, weighted["fractiles"][q_lo][i],
                weighted["fractiles"][q_hi][i],
                color="#c6dbef", alpha=0.6, zorder=1,
                label=f"{int(q_lo*100)}-{int(q_hi*100)}% fractile band")
        ax.loglog(d0, weighted["mean"][i], color="#b2182b", lw=3.2,
                  solid_capstyle="round", zorder=3, label="weighted mean")
        ax.set_xlabel("Displacement $D_0$ [m]")
        ax.set_title(
            f"Site {label} at {int(site[2]*1000)} m: one tree, "
            "four $\\sigma$ branches", fontsize=10.5)
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(fontsize=8, loc="lower left")
    axes[0].set_ylabel("Annual rate of exceedance")
    fig2.suptitle(
        "Propagating the mapping-accuracy class as fdhaCalcRSigma branches",
        fontsize=12)
    fig2.tight_layout(rect=(0, 0, 1, 0.96))
    png2 = OUT / "r_sigma_epistemic_propagation.png"
    fig2.savefig(png2, dpi=150)
    print(f"Figure written to {png2}")


if __name__ == "__main__":
    main()
