"""Demo: r_threshold_km as epistemic logic-tree branches.

Derives four job variants from the canonical ``examples/hazard_curve_minimal``
job (Youngs2003 normal models, same site), runs them, checks the MODE A/MODE B
equivalence and the weighted-mean identity, and renders a three-panel figure:

  A. map of the fault trace, the site, and the two threshold buffers;
  B. hazard curves for the four variants;
  C. machine-precision consistency checks.

The source is swapped to a demo-local **single M7.0 characteristic event** at
a low annual rate (1e-4/yr; see ``source_model_M7.xml``) so the curves plateau
near ~1e-4/yr, comparable to a published single-scenario figure, instead of
reflecting the very active Gutenberg-Richter fault (~0.02/yr). The site is
~0.39 km from the fault, so the two demo thresholds 0.1 km and 0.5 km straddle
it: the single rupture is routed distributed under 0.1 km and principal under
0.5 km (a clean binary flip). Values are demo inputs, not recommendations.
"""
from __future__ import annotations

import json
import shutil
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

from openquake.fdha.logic_tree.driver import FdhaLogicTree

HERE = Path(__file__).resolve().parent
# Reuse the canonical examples/hazard_curve_minimal job (Youngs2003 normal
# models, same site) but swap in a demo-local single M7.0 characteristic
# source at a low annual rate, so the curves plateau near ~1e-4/yr instead of
# reflecting the very active Gutenberg-Richter fault (~0.02/yr).
BASE_DEMO = HERE.parents[3] / "examples"
LOCAL_SRC = HERE / "source_model_M7.xml"       # single M7.0 scenario, 1e-4/yr
SRC_XML = "source_model.xml"                    # name the SMLT expects
SMLT_XML = "hazard_curve_minimal_source_model_logic_tree.xml"
FDHA_XML = "hazard_curve_minimal_fdha_logic_tree.xml"
JOB_INI = "hazard_curve_minimal.ini"
OUT = HERE / "out"

# Site from the examples job (unmodified) and the two demo thresholds.
SITE_LON, SITE_LAT = 16.16573727, 39.64704451
T_LOW, T_HIGH = "0.1", "0.5"
CASES = {
    "case1": [("RT_LOW", T_LOW, "1.0")],
    "case2": [("RT_HIGH", T_HIGH, "1.0")],
    "case3": [("RT_LOW", T_LOW, "0.5"), ("RT_HIGH", T_HIGH, "0.5")],
}


def threshold_level(branches):
    rows = [
        '    <logicTreeBranchingLevel branchingLevelID="bl_5_r_threshold">',
        '      <logicTreeBranchSet branchSetID="bs_5_r_threshold" '
        'uncertaintyType="fdhaCalcRThreshold">',
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
    demo-local single-M7.0 source swapped in.

    baseline (branches=None): MODE A, INI scalar r_threshold_km = 0.1.

    MODE B cases: the ``fdhaCalcRThreshold`` branch set is appended to the FDHA
    logic tree, and the INI scalar ``r_threshold_km`` is removed (the two
    mechanisms are mutually exclusive; keeping both is a configuration error).
    """
    workdir = OUT / name
    workdir.mkdir(parents=True, exist_ok=True)
    # Demo-local single-M7.0 source, named as the SMLT wrapper expects.
    shutil.copy(LOCAL_SRC, workdir / SRC_XML)
    shutil.copy(BASE_DEMO / SMLT_XML, workdir)

    lt_text = (BASE_DEMO / FDHA_XML).read_text()
    ini_text = (BASE_DEMO / JOB_INI).read_text()
    if branches:
        lt_text = lt_text.replace(
            "  </logicTree>", threshold_level(branches) + "\n  </logicTree>"
        )
        ini_text = "".join(
            line for line in ini_text.splitlines(keepends=True)
            if not line.lstrip().startswith("r_threshold_km")
        )
    (workdir / FDHA_XML).write_text(lt_text)
    (workdir / JOB_INI).write_text(ini_text)
    return workdir / JOB_INI


def run(name: str, branches=None):
    ini = make_job(name, branches)
    res = FdhaLogicTree.from_ini(str(ini)).run(outdir=OUT / name / "out")
    d0 = np.asarray(res.d0, dtype=float)
    mean = np.asarray(res.mean_rates, dtype=float)[0]
    fractiles = {q: np.asarray(v, dtype=float)[0] for q, v in res.fractiles.items()}
    manifest = json.loads((OUT / name / "out" / "manifest.json").read_text())
    return {"d0": d0, "mean": mean, "fractiles": fractiles, "manifest": manifest}


def rel_diff(a, b):
    return np.abs(a - b) / np.where(b != 0, np.abs(b), 1.0)


def read_trace():
    """Return the fault-trace lon/lat vertices from the demo source model."""
    ns = {"gml": "http://www.opengis.net/gml"}
    root = ET.parse(LOCAL_SRC).getroot()
    pos = root.find(".//gml:posList", ns).text.split()
    vals = [float(v) for v in pos]
    return np.array(vals).reshape(-1, 2)  # (n, 2) lon, lat


def dist_to_polyline_km(lon_grid, lat_grid, trace):
    """Min distance (km) from each grid point to the trace polyline."""
    kx = 111.0 * np.cos(np.radians(SITE_LAT))
    ky = 111.0
    gx = (lon_grid - SITE_LON) * kx
    gy = (lat_grid - SITE_LAT) * ky
    tx = (trace[:, 0] - SITE_LON) * kx
    ty = (trace[:, 1] - SITE_LAT) * ky
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

    print("Running baseline (MODE A, single M7.0 source, scalar r_threshold_km = 0.1)...")
    baseline = run("baseline")
    results = {}
    for name, branches in CASES.items():
        label = " + ".join(f"{v} (w={w})" for _, v, w in branches)
        print(f"Running {name} (MODE B: r_threshold_km = {label})...")
        results[name] = run(name, branches)
    case1, case2, case3 = results["case1"], results["case2"], results["case3"]

    d0 = baseline["d0"]
    average = 0.5 * case1["mean"] + 0.5 * case2["mean"]
    e_equiv = rel_diff(case1["mean"], baseline["mean"]).max()
    e_mean = rel_diff(case3["mean"], average).max()
    spread = rel_diff(case2["mean"], case1["mean"]).max()

    print()
    print(f"case1 (branch {T_LOW}, w=1) vs baseline : max rel diff = {e_equiv:.3e}")
    print(f"case3 mean vs 0.5*case1 + 0.5*case2     : max rel diff = {e_mean:.3e}")
    print(f"case1 vs case2 (sanity: must differ)    : max rel diff = {spread:.3e}")
    n3 = case3["manifest"]["n_realizations"]
    w3 = sum(b["combined_branch_weight"] for b in case3["manifest"]["branches"])
    print(f"case3 manifest: {n3} realizations, total weight = {w3}")

    print(f"low-displacement plateau: case1 (distributed) = {case1['mean'][0]:.2e}/yr, "
          f"case2 (principal) = {case2['mean'][0]:.2e}/yr  (source rate 1e-4/yr)")

    assert e_equiv == 0.0, "MODE B (0.1) must be bit-identical to MODE A scalar"
    assert e_mean < 1e-12, "case3 mean must equal the weighted average"
    assert spread > 0.1, "thresholds must actually change the routing"

    # ------------------------------------------------------------- figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    fig, (axm, ax1, ax2) = plt.subplots(
        1, 3, figsize=(17, 5.4), gridspec_kw={"width_ratios": [2.1, 3, 2]}
    )

    # -- Panel A: map of trace, site, threshold buffers --------------------
    trace = read_trace()
    margin = 0.02
    lon_min, lon_max = trace[:, 0].min() - margin, trace[:, 0].max() + margin
    lat_min, lat_max = trace[:, 1].min() - margin, trace[:, 1].max() + margin
    LON, LAT = np.meshgrid(
        np.linspace(lon_min, lon_max, 400), np.linspace(lat_min, lat_max, 400)
    )
    dist = dist_to_polyline_km(LON, LAT, trace)
    # Bands: <=0.1 principal under both; 0.1-0.5 principal only under 0.5.
    axm.contourf(
        LON, LAT, dist, levels=[0, 0.1, 0.5, 1e9],
        colors=["#b2182b", "#f4a582", "#e8e8e8"], alpha=0.75,
    )
    axm.plot(trace[:, 0], trace[:, 1], "k-", lw=2.4, marker="o", ms=4,
             label="fault trace")
    axm.plot(SITE_LON, SITE_LAT, marker="*", ms=20, mfc="#2166ac",
             mec="white", mew=1.2, ls="none", label="site")
    axm.annotate(
        "nearest rupture\n$R_{JB}\\approx0.39$ km", xy=(SITE_LON, SITE_LAT),
        xytext=(SITE_LON - 0.017, SITE_LAT - 0.028), fontsize=9,
        ha="center", arrowprops=dict(arrowstyle="->", lw=1.1),
    )
    axm.set_aspect(1.0 / np.cos(np.radians(SITE_LAT)))
    axm.set_xlabel("Longitude"); axm.set_ylabel("Latitude")
    axm.set_title("Where is the site?")
    axm.legend(
        handles=[
            Line2D([], [], color="k", lw=2.4, marker="o", ms=4, label="fault trace"),
            Line2D([], [], color="#2166ac", marker="*", ms=14, ls="none",
                   mec="white", label="site (~0.39 km out)"),
            Patch(fc="#b2182b", alpha=0.75, label="$\\leq$0.1 km: principal (both)"),
            Patch(fc="#f4a582", alpha=0.75, label="0.1-0.5 km: principal only if thr=0.5"),
            Patch(fc="#e8e8e8", alpha=0.75, label="$>$0.5 km: distributed (both)"),
        ],
        fontsize=8, loc="upper right",
    )

    # -- Panel B: hazard curves --------------------------------------------
    ax1.loglog(d0, baseline["mean"], color="0.55", lw=6, alpha=0.6,
               label="baseline (MODE A scalar 0.1 km)")
    ax1.loglog(d0, case1["mean"], "C0-", lw=1.8,
               label=f"case 1: branch {T_LOW} km (w=1.0)")
    ax1.loglog(d0, case2["mean"], "C3-", lw=1.8,
               label=f"case 2: branch {T_HIGH} km (w=1.0)")
    ax1.loglog(d0, case3["mean"], "C2-", lw=2.4,
               label=f"case 3 mean: {{{T_LOW}, {T_HIGH}}} (w=0.5 each)")
    ax1.loglog(d0, average, "k--", lw=1.2, marker="o", ms=4, mfc="none",
               label="0.5*case1 + 0.5*case2 (check)")
    if 0.05 in case3["fractiles"] and 0.95 in case3["fractiles"]:
        ax1.fill_between(d0, case3["fractiles"][0.05], case3["fractiles"][0.95],
                         color="C2", alpha=0.12, label="case 3 p05-p95")
    ax1.set_xlabel("Displacement $D_0$ [m]")
    ax1.set_ylabel("Annual rate of exceedance")
    ax1.set_title("Hazard curves (single M7.0 scenario, rate 1e-4/yr)")
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(fontsize=8.5, loc="lower left")

    # -- Panel C: consistency ----------------------------------------------
    ax2.semilogx(d0, rel_diff(case3["mean"], average), "ko-", ms=4,
                 label="|case3 $-$ avg| / avg")
    ax2.semilogx(d0, rel_diff(case1["mean"], baseline["mean"]), "C0s-", ms=4,
                 label="|case1 $-$ baseline| / baseline")
    ax2.set_ylim(-1e-16, max(1e-15, ax2.get_ylim()[1]))
    ax2.set_xlabel("Displacement $D_0$ [m]")
    ax2.set_ylabel("Relative difference")
    ax2.set_title("Consistency checks (0 = exact)")
    ax2.grid(True, which="both", alpha=0.25)
    ax2.legend(fontsize=9)

    fig.tight_layout()
    png = OUT / "r_threshold_epistemic_demo.png"
    fig.savefig(png, dpi=150)
    print(f"\nFigure written to {png}")


if __name__ == "__main__":
    main()
