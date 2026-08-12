#!/usr/bin/env python
"""Paper figures for the Norcia epistemic-anatomy application section.

Reads the outputs of job_anatomy_curve.ini (45 end branches, 3 sites),
job_profile_transect.ini (101-site cross-trace profile), and
job_anatomy_map.ini (45 branches, 0.01 deg grid), and renders:

  fig_ms_tunnel        distributed-model comparison at the Mount Serra
                       site with the observed 2016 San Benedetto tunnel
                       offset (Galli et al. 2019 WTC: 20 cm vertical,
                       13 cm left-lateral). One curve per distributed
                       CHAIN (model x its published displacement
                       variant), each the weighted mean over the primary
                       nodes that gate it. No aggregate mean is drawn -
                       the figure is a model comparison, not a hazard
                       result.
  fig_three_site_fans  per-branch spaghetti + weighted mean + 5-95%
                       band at PF / MS / SL
  fig_profile          displacement at the 1e5 yr return period across
                       the MVF trace under the Petersen Gaussian W_p
                       (sigma = 26.89 m, "Accurately located")
The near-trace corridor is characterised by fig_profile. The regional
map appears only as panel (a) of fig_ms_tunnel, and shows the PRINCIPAL
component alone.

Colors are the Okabe & Ito (2008) colorblind-safe palette: one hue
family per model, distinct hues for its published variants; solid lines
throughout except the 100 m cell-size sensitivity.
"""

from __future__ import annotations

import csv
import json
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
FIGDIR = HERE / "figures"
FIGDIR.mkdir(exist_ok=True)

# ----------------------------------------------------------------- style
plt.rcParams.update({
    "font.size": 8.5,
    "axes.labelsize": 9,
    "axes.titlesize": 9.5,
    "legend.fontsize": 7.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.linewidth": 0.8,
    "grid.linewidth": 0.4,
    "grid.alpha": 0.35,
    "lines.linewidth": 1.4,
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
})

# Okabe-Ito hues: one family hue per model, lightness steps within a
# family for its published variants
C_Y03 = "#0072B2"      # blue      - Youngs et al. (2003), 85th
C_Y03_95 = "#56B4E9"   # sky       - Youngs et al. (2003), 95th
C_V25 = "#D55E00"      # vermilion - Visini et al. (2025), WC1994
C_V25_TB = "#E69F00"   # orange    - Visini et al. (2025), Thingbaijam
C_V25_LE = "#CC79A7"   # purple    - Visini et al. (2025), Leonard
C_MEAN = "#000000"

# Distributed cell convention. Youngs et al. (2003) and Ferrario &
# Livio (2021) are both developed for a 500 x 500 m cell and expose no
# cell parameter, so the tree pins Visini et al. (2025) - the only
# cell-parameterised model - to 500 m. The 100 m variant (the IAEA site
# dimension) is carried as an explicit sensitivity from out_cellsize.
CELL_TREE_M = 500
CELL_ALT_M = 100

SITES = {0: "PF", 1: "MS", 2: "SL"}
SITE_TITLES = {
    0: "PF (on the MVF trace)",
    1: "MS ($r$ = 7.6 km)",
    2: "SL ($r$ = 2.4 km)",
}

# observed 2016 offset at the San Benedetto tunnel (Galli et al. 2019 WTC
# proceedings / Galli 2020): 20 cm vertical, 13 cm left-lateral
D_OBS_V = 0.20
PF_SITE, MS_SITE, SL_SITE = 0, 1, 2
# what was seen at each antithetic site in 2016
OBSERVED = {
    MS_SITE: (0.20, "Observed 2016 San Benedetto\ntunnel offset, 0.20 m vertical"),
    SL_SITE: (1.00, "Observed 2016 road\noffset, ~1 m"),
}
D_OBS_NET = float(np.hypot(0.20, 0.13))

RP_TARGET = 1.0e5  # IAEA reference return period (a)

# ------------------------------------------------- node weights (tree)
W_PSR = {"PIZZA23": 0.4, "Y03EC": 0.2, "Y03NBR": 0.2, "Y03GB": 0.2}
W_PFD = {"Y03_AD": 0.25, "Y03_MD": 0.25, "LAV23": 0.5}
W_SSR = {"Y03": 0.5, "V25": 0.5}
W_SFD = {"Y03_P85": 0.5, "Y03_P95": 0.5,
         "V25_WC94": 0.33, "V25_TB17": 0.33, "V25_LEO10": 0.34}
W_RS = {"ACCURATE": 1.0}

#: The three figures the paper carries. The other analyses below still
#: run, because their numbers are quoted in `discussion_notes.md` and in
#: the section text, but they are not rendered to disk.
PAPER_FIGURES = {"fig_logic_tree", "fig_ms_tunnel", "fig_principal_fractiles"}


def save(fig, name):
    """Write a figure only if the paper carries it."""
    if name in PAPER_FIGURES:
        for ext in ("pdf", "png"):
            fig.savefig(FIGDIR / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)


# every distributed chain in the tree, as (S_sr, S_fd) -> label, colour
# Only the WC1994 Visini variant is drawn: its Thingbaijam2017 and
# Leonard2010 siblings differ by x1.1 and overplot it. Youngs is shown at
# its 85th-percentile fit. The tree still carries every variant; this is
# a plotting selection, not a change to the weights.
CHAINS = [
    (("Y03", "Y03_P85"), "Youngs et al. (2003) - 85th, 500 m", C_Y03),
    (("Y03", "Y03_P95"), "Youngs et al. (2003) - 95th, 500 m", C_Y03_95),
    (("V25", "V25_WC94"), "Visini, 500 m cell (pre-existing fault)", C_V25),
]
#: displacement variant whose 100 m arm is drawn as the cell-size
#: sensitivity - must match a CHAINS entry to be a controlled comparison
CELL_ALT_SFD = "V25_WC94"


def parse_branch_id(bid: str, r_sigma: float) -> dict:
    """Decompose a composed fdha_branch_id into anatomy node levels."""
    tok = {t.split("_", 1)[0]: t for t in bid.split("|")}
    psr = tok["PSR"].replace("PSR_", "")
    psr = psr.replace("_YDEF", "").replace("_LAV", "")
    node = {
        "PSR": psr,
        "PFD": tok["PFD"].replace("PFD_", ""),
        "SSR": tok["SSR"].replace("SSR_", ""),
        "SFD": tok["SFD"].replace("SFD_", ""),
        "RS": tok["RS"].replace("RS_", ""),
    }
    assert (node["RS"] == "BOXCAR") == (r_sigma == 0.0), (bid, r_sigma)
    return node


def branch_weight(node: dict) -> float:
    return (W_PSR[node["PSR"]] * W_PFD[node["PFD"]] * W_SSR[node["SSR"]]
            * W_SFD[node["SFD"]] * W_RS[node["RS"]])


def load_curves(outdir: Path):
    """Return (d0, branches) with per-branch total/principal/distributed
    rate arrays of shape (nsites, nlevels)."""
    manifest = json.loads((outdir / "manifest.json").read_text())
    branches = []
    d0 = None
    for b in manifest["branches"]:
        rows = list(csv.DictReader(open(outdir / b["curve_file"])))
        sids = sorted({int(r["site_id"]) for r in rows})
        if d0 is None:
            d0 = np.array(sorted({float(r["D0"]) for r in rows}))
        nd = len(d0)
        tot = np.zeros((len(sids), nd))
        pri = np.zeros((len(sids), nd))
        dis = np.zeros((len(sids), nd))
        for r in rows:
            i = int(r["site_id"])
            j = int(np.searchsorted(d0, float(r["D0"])))
            tot[i, j] = float(r["annual_rate"])
            pri[i, j] = float(r["annual_rate_principal"])
            dis[i, j] = float(r["annual_rate_distributed"])
        rs = float(b.get("fdha_calc_params", {}).get("r_sigma_km", 0.0))
        entry = {
            "id": b["fdha_branch_id"], "rates": tot,
            "principal": pri, "distributed": dis, "r_sigma": rs,
        }
        if "PSR_" in b["fdha_branch_id"]:
            entry["node"] = parse_branch_id(b["fdha_branch_id"], rs)
            entry["w"] = branch_weight(entry["node"])
        branches.append(entry)
    return d0, branches


def rate_at(d0, rates, d):
    """Log-log interpolated annual rate at displacement d (per site row)."""
    r = np.asarray(rates, dtype=float)
    out = np.full(r.shape[:-1], np.nan)
    it = np.ndindex(r.shape[:-1])
    for idx in it:
        y = r[idx]
        pos = y > 0
        if pos.sum() < 2 or d > d0[pos].max() or d < d0[pos].min():
            out[idx] = 0.0
            continue
        out[idx] = np.exp(np.interp(np.log(d), np.log(d0[pos]),
                                    np.log(y[pos])))
    return out


def displ_at_rate(d0, rates, target):
    """Displacement whose exceedance rate equals target (log-log interp);
    0 where the whole curve is below target."""
    y = np.asarray(rates, dtype=float)
    pos = y > 0
    if pos.sum() < 2 or y.max() < target:
        return 0.0
    lx, ly = np.log(d0[pos]), np.log(y[pos])
    # rates decrease with displacement: reverse for np.interp
    return float(np.exp(np.interp(np.log(target), ly[::-1], lx[::-1])))


def fault_traces():
    """(name, trace xy) for every source with a simpleFaultGeometry."""
    gml = "{http://www.opengis.net/gml}"
    root = ET.parse(HERE / "source_model_norcia_case3.xml").getroot()
    traces = []
    for src in root.iter():
        if not src.tag.endswith("Source"):
            continue
        pls = src.findall(f".//{gml}posList")
        if not pls:
            continue
        xy = np.array(pls[0].text.split(), dtype=float).reshape(-1, 2)
        traces.append((src.get("name") or src.get("id"), xy))
    return traces


SITE_XY = {0: (13.278, 42.767), 1: (13.188, 42.749), 2: (13.212, 42.853)}
# marker fills matching the submission figure (PF black, MS blue, SL orange)
SITE_COLORS = {0: "0.15", 1: "#0072B2", 2: "#E69F00"}


def add_return_period_axis(ax):
    sec = ax.secondary_yaxis(
        "right", functions=(lambda r: 1.0 / np.maximum(r, 1e-300),
                            lambda t: 1.0 / np.maximum(t, 1e-300)))
    sec.set_ylabel("return period (a)")
    return sec


def load_cellsize_curves(d0_ref, subdir="out_cellsize", site_id=1):
    """Weighted-mean distributed curve at one site per Visini cell,
    from a companion cell-size run. Returns {} if not run."""
    outdir = HERE / subdir
    if not (outdir / "manifest.json").exists():
        return {}
    manifest = json.loads((outdir / "manifest.json").read_text())
    acc = {}
    for b in manifest["branches"]:
        cell = 100 if "SSR_V25_C100" in b["fdha_branch_id"] else 500
        node = parse_branch_id(
            b["fdha_branch_id"].replace("SSR_V25_C100", "SSR_V25")
                              .replace("SSR_V25_C500", "SSR_V25"),
            float(b.get("fdha_calc_params", {}).get("r_sigma_km", 0.0)))
        # weight within one cell: primary x displacement pairing only
        if node["SFD"] != CELL_ALT_SFD:
            continue  # compare like with like against the CHAINS entry
        w = (W_PSR[node["PSR"]] * W_PFD[node["PFD"]])
        rows = [r for r in csv.DictReader(open(outdir / b["curve_file"]))
                if int(r["site_id"]) == site_id]
        y = np.zeros(len(d0_ref))
        for r in rows:
            y[int(np.searchsorted(d0_ref, float(r["D0"])))] = float(
                r["annual_rate_distributed"])
        tot_w, tot_y = acc.get(cell, (0.0, np.zeros(len(d0_ref))))
        acc[cell] = (tot_w + w, tot_y + w * y)
    return {c: y / w for c, (w, y) in acc.items()}


def assign_to_traces(tlon, tlat):
    """Assign on-trace sites to their nearest declared fault trace and
    return, per trace, the site indices ordered by along-strike arc
    length together with that arc length in km.

    The Norcia traces overlap in latitude, so latitude is not a usable
    ordering variable: sorting by it interleaves the two faults. This is
    the same nearest-trace/arc-length assignment used to colour the
    traces on the map, factored out so both can share it.
    """
    traces = fault_traces()
    if not len(tlon) or not traces:
        return []
    lat0 = float(np.mean(tlat))
    kx, ky = 111.195 * np.cos(np.radians(lat0)), 111.195
    P = np.column_stack((tlon * kx, tlat * ky))

    polys, best = [], np.full(len(P), np.inf)
    owner = np.full(len(P), -1)
    for k, (_name, xy) in enumerate(traces):
        q = np.column_stack((xy[:, 0] * kx, xy[:, 1] * ky))
        polys.append(q)
        seg = np.diff(q, axis=0)
        L2 = np.einsum("ij,ij->i", seg, seg)
        L2[L2 == 0] = 1e-12
        w = P[:, None, :] - q[None, :-1, :]
        t = np.clip(np.einsum("ijk,jk->ij", w, seg) / L2, 0.0, 1.0)
        proj = q[None, :-1, :] + t[..., None] * seg[None, :, :]
        d = np.linalg.norm(P[:, None, :] - proj, axis=2)
        j = np.argmin(d, axis=1)
        dmin = d[np.arange(len(P)), j]
        take = dmin < best
        best[take], owner[take] = dmin[take], k

    out = []
    for k, q in enumerate(polys):
        sel = np.flatnonzero(owner == k)
        if sel.size < 2:
            continue
        seg = np.diff(q, axis=0)
        L2 = np.einsum("ij,ij->i", seg, seg)
        L2[L2 == 0] = 1e-12
        cum = np.concatenate(([0.0], np.cumsum(np.sqrt(L2))))
        w = P[sel][:, None, :] - q[None, :-1, :]
        t = np.clip(np.einsum("ijk,jk->ij", w, seg) / L2, 0.0, 1.0)
        proj = q[None, :-1, :] + t[..., None] * seg[None, :, :]
        d = np.linalg.norm(P[sel][:, None, :] - proj, axis=2)
        j = np.argmin(d, axis=1)
        sarc = cum[j] + t[np.arange(sel.size), j] * np.sqrt(L2[j])
        o = np.argsort(sarc)
        out.append((fault_traces()[k][0], sel[o], sarc[o]))
    return out


def _draw_colored_traces(ax, tlon, tlat, tval, cmap, norm):
    """Draw each fault trace as a dashed line coloured by the on-trace
    displacement, replacing the flat-colour trace and the separate
    on-trace markers.

    The on-trace map sites are an unordered flat list, so each is
    assigned to its nearest declared trace and ordered by along-strike
    arc length before the segments are built.
    """
    from matplotlib.collections import LineCollection

    traces = fault_traces()
    if not len(tlon) or not traces:
        return
    lat0 = float(np.mean(tlat))
    kx, ky = 111.195 * np.cos(np.radians(lat0)), 111.195
    P = np.column_stack((tlon * kx, tlat * ky))

    polys, best = [], np.full(len(P), np.inf)
    owner = np.full(len(P), -1)
    for k, (_name, xy) in enumerate(traces):
        q = np.column_stack((xy[:, 0] * kx, xy[:, 1] * ky))
        polys.append(q)
        seg = np.diff(q, axis=0)
        L2 = np.einsum("ij,ij->i", seg, seg)
        L2[L2 == 0] = 1e-12
        # perpendicular distance from every site to every segment
        w = P[:, None, :] - q[None, :-1, :]
        t = np.clip(np.einsum("ijk,jk->ij", w, seg) / L2, 0.0, 1.0)
        proj = q[None, :-1, :] + t[..., None] * seg[None, :, :]
        d = np.linalg.norm(P[:, None, :] - proj, axis=2)
        j = np.argmin(d, axis=1)
        dmin = d[np.arange(len(P)), j]
        take = dmin < best
        best[take], owner[take] = dmin[take], k

    for k, q in enumerate(polys):
        sel = owner == k
        if sel.sum() < 2:
            continue
        seg = np.diff(q, axis=0)
        L2 = np.einsum("ij,ij->i", seg, seg)
        L2[L2 == 0] = 1e-12
        cum = np.concatenate(([0.0], np.cumsum(np.sqrt(L2))))
        w = P[sel][:, None, :] - q[None, :-1, :]
        t = np.clip(np.einsum("ijk,jk->ij", w, seg) / L2, 0.0, 1.0)
        proj = q[None, :-1, :] + t[..., None] * seg[None, :, :]
        d = np.linalg.norm(P[sel][:, None, :] - proj, axis=2)
        j = np.argmin(d, axis=1)
        sarc = cum[j] + t[np.arange(sel.sum()), j] * np.sqrt(L2[j])

        order = np.argsort(sarc)
        pts = np.column_stack((tlon[sel][order], tlat[sel][order]))
        val = tval[sel][order]
        segs = np.stack([pts[:-1], pts[1:]], axis=1)
        lc = LineCollection(
            segs, cmap=cmap, norm=norm, linewidths=1.8, zorder=4,
            linestyles=[(0, (4.0, 2.6))],
            path_effects=[pe.Stroke(linewidth=3.2, foreground="0.2"),
                          pe.Normal()])
        lc.set_array(0.5 * (val[:-1] + val[1:]))
        ax.add_collection(lc)


# ------------------------------------------------------------ figure 1
def fig_ms_tunnel(d0, branches):
    """Distributed hazard at the two antithetic-fault sites, MS and SL,
    against what was actually observed there in 2016.

    Both sites sit on pre-existing antithetic structures of the Mount
    Vettore system that ruptured at the surface in 2016, which is why
    both are classified case1 (all three Visini combinations active).
    The regional map lives in fig_principal_fractiles; this figure is
    the site comparison alone.
    """
    cells = {}
    cells3 = {}
    for sid in (MS_SITE, SL_SITE):
        cells[sid] = load_cellsize_curves(d0, "out_cellsize", sid)
        cells3[sid] = load_cellsize_curves(d0, "out_cellsize_case3", sid)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True)
    for n, (sid, ax) in enumerate(zip((MS_SITE, SL_SITE), axes)):
        curves = {}
        for key, _lbl, _c in CHAINS:
            sel = [b for b in branches
                   if (b["node"]["SSR"], b["node"]["SFD"]) == key]
            wb = np.array([b["w"] for b in sel])
            arr = np.array([b["distributed"][sid] for b in sel])
            curves[key] = (wb[:, None] * arr).sum(0) / wb.sum()

        for key, lbl, colour in CHAINS:
            ax.plot(d0, np.maximum(curves[key], 1e-16), "-", color=colour,
                    lw=1.5, label=lbl)
        # Visini under both decision-tree classifications, labelled by
        # what they physically assume rather than by the paper's
        # combination letters:
        #   "pre-existing fault"     = case1 = combinations A + B + C
        #   "no pre-existing fault"  = case3 = combination A alone
        # Solid = 500 m cell, dashed = 100 m. With a pre-existing fault
        # the two cells overplot; without one they separate, because
        # only A carries the cell-dependent along-strike Monte Carlo
        # term (C sets the along-strike probability to 1 by design,
        # Visini et al. 2025 p. 12).
        if cells[sid]:
            ax.plot(d0, np.maximum(cells[sid][CELL_ALT_M], 1e-16), "--",
                    color=C_V25, lw=1.4,
                    label=f"Visini, {CELL_ALT_M} m cell (pre-existing fault)")
        if cells3[sid]:
            ax.plot(d0, np.maximum(cells3[sid][CELL_TREE_M], 1e-16), "-",
                    color=C_V25_TB, lw=1.5,
                    label=f"Visini, {CELL_TREE_M} m cell (no pre-existing fault)")
            ax.plot(d0, np.maximum(cells3[sid][CELL_ALT_M], 1e-16), "--",
                    color=C_V25_TB, lw=1.4,
                    label=f"Visini, {CELL_ALT_M} m cell (no pre-existing fault)")

        obs, note = OBSERVED[sid]
        ax.axvline(obs, color="0.15", lw=1.0, ls="--", zorder=3)
        ax.text(0.03, 0.03, note, transform=ax.transAxes,
                ha="left", va="bottom", linespacing=1.35)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1e-3, 3)
        ax.set_ylim(1e-11, 1e-3)
        ax.set_xlabel("Fault displacement (m)")
        ax.set_title(f"({'ab'[n]}) Mean hazard curves at {SITES[sid]}",
                     loc="left")
        ax.grid(True, which="both")
        if n == 0:
            ax.set_ylabel("Annual frequency of exceedance (1/yr)")

        print(f"\n=== {SITES[sid]}: rate/RP of exceeding the observed "
              f"{obs} m ===")
        for key, lbl, _c in CHAINS:
            r = rate_at(d0, curves[key][None, :], obs)[0]
            print(f"  {lbl:42s} {r:.3e} /a  RP {1/r:,.0f} a")
        if cells3[sid]:
            for cell in (CELL_TREE_M, CELL_ALT_M):
                r3 = rate_at(d0, cells3[sid][cell][None, :], obs)[0]
                print(f"  {'Visini ' + str(cell) + ' m, no pre-existing fault':42s} "
                      f"{r3:.3e} /a  RP {1/r3:,.0f} a")
        if cells[sid]:
            ra = rate_at(d0, cells[sid][CELL_ALT_M][None, :], obs)[0]
            print(f"  {'Visini, 100 m cell, pre-existing fault':42s} {ra:.3e} /a  "
                  f"RP {1/ra:,.0f} a")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.06), handlelength=2.4,
               columnspacing=1.6)
    save(fig, "fig_ms_tunnel")


# ------------------------------------------------------------ figure 2
def fig_three_site_fans(d0, branches):
    w = np.array([b["w"] for b in branches])
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.7), sharey=True)
    for i, ax in enumerate(axes):
        allr = np.array([b["rates"][i] for b in branches])
        for y in allr:
            ax.plot(d0, np.maximum(y, 1e-16), "-", color="0.6", lw=0.4,
                    alpha=0.35)
        mean = (w[:, None] * allr).sum(0) / w.sum()
        # weighted fractiles per displacement level
        q05 = np.zeros_like(mean)
        q95 = np.zeros_like(mean)
        for j in range(len(d0)):
            order = np.argsort(allr[:, j])
            cw = np.cumsum(w[order]) / w.sum()
            q05[j] = allr[order][np.searchsorted(cw, 0.05), j]
            q95[j] = allr[order][np.searchsorted(cw, 0.95), j]
        ax.fill_between(d0, np.maximum(q05, 1e-16),
                        np.maximum(q95, 1e-16), color=C_Y03, alpha=0.18,
                        lw=0)
        ax.plot(d0, np.maximum(mean, 1e-16), color=C_MEAN, lw=1.8)
        ax.axhline(1.0 / RP_TARGET, color="0.3", lw=0.7, ls=":")
        if i == 2:
            ax.annotate("$10^{5}$ a", xy=(2.5, 1.1e-5), fontsize=7,
                        color="0.3")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1e-3, 10)
        ax.set_ylim(1e-11, 1e-2)
        ax.set_xlabel("displacement (m)")
        ax.set_title(SITE_TITLES[i], loc="left", fontsize=8)
        ax.grid(True, which="both")
        # spread annotation at the 1e5 a level
        dmin = displ_at_rate(d0, q05, 1.0 / RP_TARGET)
        dmax = displ_at_rate(d0, q95, 1.0 / RP_TARGET)
        print(f"site {SITES[i]}: D(1e5a) q05={dmin:.3g} q95={dmax:.3g} "
              f"mean={displ_at_rate(d0, mean, 1.0/RP_TARGET):.3g} m")
    axes[0].set_ylabel("annual rate of exceedance (a$^{-1}$)")
    handles = [
        Line2D([], [], color="0.6", lw=0.6,
               label=f"{len(branches)} end branches"),
        Line2D([], [], color=C_MEAN, lw=1.8, label="weighted mean"),
        plt.Rectangle((0, 0), 1, 1, fc=C_Y03, alpha=0.18,
                      label="5-95% weighted band"),
    ]
    axes[0].legend(handles=handles, loc="lower left", frameon=False)
    fig.tight_layout()
    save(fig, "fig_three_site_fans")


# ------------------------------------------------------------ figure 3
def fig_profile():
    """Cross-trace hazard profile under the Petersen Gaussian W_p."""
    outdir = HERE / "out_transect"
    manifest = json.loads((outdir / "manifest.json").read_text())
    offs = np.loadtxt(HERE / "transect_offsets_km.txt") * 1000.0  # m
    b = manifest["branches"][0]
    rows = list(csv.DictReader(open(outdir / b["curve_file"])))
    d0 = np.array(sorted({float(r["D0"]) for r in rows}))
    tot = np.zeros((len(offs), len(d0)))
    pri = np.zeros_like(tot)
    dis = np.zeros_like(tot)
    for r in rows:
        i = int(r["site_id"])
        j = int(np.searchsorted(d0, float(r["D0"])))
        tot[i, j] = float(r["annual_rate"])
        pri[i, j] = float(r["annual_rate_principal"])
        dis[i, j] = float(r["annual_rate_distributed"])

    dRP = np.array([displ_at_rate(d0, tot[i], 1.0 / RP_TARGET)
                    for i in range(len(offs))])
    dRP_p = np.array([displ_at_rate(d0, pri[i], 1.0 / RP_TARGET)
                      for i in range(len(offs))])
    dRP_d = np.array([displ_at_rate(d0, dis[i], 1.0 / RP_TARGET)
                      for i in range(len(offs))])
    r02 = rate_at(d0, tot, D_OBS_V)

    def vis(a):
        """Mask absent components so the line stops instead of
        flatlining on the axis floor."""
        return np.where(a > 0, a, np.nan)

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8))
    axes[0].plot(offs, vis(dRP), color=C_MEAN, lw=2.6, label="total",
                 zorder=1)
    axes[0].plot(offs, vis(dRP_p), color=C_Y03, lw=1.4,
                 label="principal", zorder=3)
    axes[0].plot(offs, vis(dRP_d), color=C_V25, lw=1.4,
                 label="distributed", zorder=3)
    axes[0].set_ylim(1e-2, 8)
    axes[1].plot(offs, vis(r02), color=C_MEAN, lw=1.8)

    for ax in axes:
        for x in (-2 * 26.89, 2 * 26.89):
            ax.axvline(x, color="0.4", lw=0.8, ls=":")
        ax.set_yscale("log")
        ax.set_xlim(-600, 600)
        ax.set_xlabel("distance from MVF trace (m)\n"
                      "(negative = hanging wall)")
        ax.grid(True, which="both")
    from matplotlib.transforms import blended_transform_factory
    for ax in axes:
        tr = blended_transform_factory(ax.transData, ax.transAxes)
        ax.text(60, 0.93, "$\\pm 2\\sigma$", transform=tr, fontsize=7,
                color="0.35", ha="left")
    axes[0].set_ylabel("displacement at $10^{5}$ a (m)")
    axes[0].set_title("(a) design displacement", loc="left")
    axes[0].legend(loc="lower left", frameon=False, fontsize=7)
    axes[1].set_ylabel("rate of exceeding 0.20 m (a$^{-1}$)")
    axes[1].set_title("(b) rate at the tunnel-offset level", loc="left")
    fig.tight_layout()
    save(fig, "fig_profile")

    above1 = offs[dRP >= 1.0]
    print(f"profile: peak {dRP.max():.2f} m; D>=1 m from {above1.min():.0f} "
          f"to {above1.max():.0f} m ({above1.max() - above1.min():.0f} m wide)")
    for x in (-300, -100, 100, 300):
        i = int(np.argmin(np.abs(offs - x)))
        print(f"  D(1e5a) at {x:+4d} m: {dRP[i]:.3f} m")


def _load_map(outdir: Path, name: str):
    """(lon, lat, value) from a displacement-map CSV."""
    rows = list(csv.DictReader(
        (line for line in open(outdir / name) if not line.startswith("#"))))
    lon = np.array([float(r["lon"]) for r in rows])
    lat = np.array([float(r["lat"]) for r in rows])
    col = [c for c in rows[0] if c.startswith("displ")][0]
    return lon, lat, np.array([float(r[col]) for r in rows])


def _raster(lon, lat, step):
    """Bin scattered engine-grid samples onto a regular raster.

    The engine's region grid is not an exact lon/lat lattice (row-wise
    longitudes drift), so bin by index rather than np.unique.
    """
    ix = np.round((lon - lon.min()) / step).astype(int)
    iy = np.round((lat - lat.min()) / step).astype(int)
    ulon = lon.min() + step * np.arange(ix.max() + 1)
    ulat = lat.min() + step * np.arange(iy.max() + 1)

    def grid(v):
        g = np.full((len(ulat), len(ulon)), np.nan)
        g[iy, ix] = v
        return g

    return ulon, ulat, grid


# ------------------------------------------------------------ figure 4
def principal_fractile_maps(quantiles=(0.84,), with_mean=True):
    """Weighted EPISTEMIC fractiles of the principal displacement at the
    map return period, rebuilt from the per-branch curves.

    The aggregate ``rates_fractiles.h5`` stores TOTAL rates only, and the
    exported quantile CSVs likewise carry no principal column, so a
    principal-only fractile cannot be read off them: it has to be
    recomputed from ``rates_principal`` in each branch file, weighted by
    the branch weight.

    The mean is the displacement read from the weighted-mean RATE curve
    (the standard mean hazard curve), not the mean of the per-branch
    displacements - the two differ, and the former is what the equal-
    hazard design value means.

    The fractiles are taken ACROSS LOGIC-TREE BRANCHES - epistemic. They are
    a different quantity from the aleatory percentiles of the
    displacement distribution given rupture (fig_design_basis), and the
    two must not be read as versions of one another.
    """
    import h5py
    files = sorted((HERE / "out_map" / "source_model_branches").glob(
        "*/branches/branch_*.h5"))
    if not files:
        return None
    curves, wts = [], []
    for fp in files:
        with h5py.File(fp, "r") as f:
            d0m = f["d0"][:]
            curves.append(f["rates_principal"][:])
            wts.append(float(f.attrs["weight"]))
            lons, lats = f["site_lons"][:], f["site_lats"][:]
    C = np.stack(curves)                       # (branch, site, displ)
    w = np.asarray(wts, dtype=float)
    w /= w.sum()

    # displacement at the target rate, per branch and site
    target = 1e-5
    nb, ns, _ = C.shape
    D = np.zeros((nb, ns))
    for b in range(nb):
        for i in range(ns):
            r = C[b, i]
            m = r > 0
            if m.sum() < 2 or target > r[m][0]:
                continue                        # left clamp: D = 0
            if target < r[m][-1]:
                D[b, i] = d0m[m][-1]
                continue
            D[b, i] = np.exp(np.interp(np.log(target), np.log(r[m])[::-1],
                                       np.log(d0m[m])[::-1]))
    out = {}
    if with_mean:
        # displacement off the weighted-mean rate curve
        Cm = np.tensordot(w, C, axes=(0, 0))
        dm = np.zeros(ns)
        for i in range(ns):
            r = Cm[i]
            m = r > 0
            if m.sum() < 2 or target > r[m][0]:
                continue
            dm[i] = (d0m[m][-1] if target < r[m][-1] else
                     np.exp(np.interp(np.log(target), np.log(r[m])[::-1],
                                      np.log(d0m[m])[::-1])))
        out["mean"] = dm

    # weighted quantile per site, over branches
    order = np.argsort(D, axis=0)
    Ds = np.take_along_axis(D, order, axis=0)
    Ws = w[order]
    cw = np.cumsum(Ws, axis=0)
    cw /= cw[-1]
    for q in quantiles:
        idx = (cw >= q).argmax(axis=0)
        out[q] = Ds[idx, np.arange(ns)]
    return lons, lats, out


def fig_design_basis():
    """Design-displacement basis for the on-fault principal hazard,
    following Abrahamson and Liou (2025, BSSA 115, 2845-2856).

    That paper distinguishes the equal-hazard displacement (the value at
    the design return period read off the MEAN hazard curve) from three
    alternatives that become necessary when the hazard curve is FLAT at
    the design return period, in which case the mean-hazard value is a
    poor representation of the epistemic centre. Their scope statement
    is explicit (p. 2853): the conditional approach "is intended to be
    used for cases in which the hazard curve is flat at the selected
    return period, and not to reduce the displacement for cases in which
    the hazard curve is steep."

    So the slope is the deciding diagnostic, and it is computed and
    plotted here rather than assumed. Panel (b) shows the CCDF of
    equation (9), P(D > z | D > 0) = Haz(D > z)/Haz(D > 0), from which
    any percentile given surface rupture can be read.
    """
    fmap = HERE / "out_map" / "aggregate" / "rates_mean.h5"
    if not fmap.exists():
        return
    import h5py
    with h5py.File(fmap, "r") as f:
        d0m = f["d0"][:]
        R = f["rates_mean_principal"][:]
    mp = HERE / "out_map" / "aggregate" / "displacement_map_mean.csv"
    rows = [r for r in csv.DictReader(
        l for l in open(mp) if not l.startswith("#"))]
    mtr = np.array([r.get("is_trace", "0") == "1" for r in rows])
    idx = np.flatnonzero(mtr)
    Rt = R[idx]
    tlon = np.array([float(r["lon"]) for r in rows])[idx]
    tlat = np.array([float(r["lat"]) for r in rows])[idx]
    per_fault = assign_to_traces(tlon, tlat)

    def interp_d(curve, target):
        m = curve > 0
        if m.sum() < 2 or target > curve[m][0] or target < curve[m][-1]:
            return np.nan
        return float(np.exp(np.interp(np.log(target),
                                      np.log(curve[m])[::-1],
                                      np.log(d0m[m])[::-1])))

    eq = np.array([interp_d(r, 1e-5) for r in Rt])
    ccdf = Rt / Rt[:, :1]                      # equation (9)
    p50 = np.array([interp_d(c, 0.50) for c in ccdf])
    p85 = np.array([interp_d(c, 0.15) for c in ccdf])
    rup = Rt[:, 0]                             # Haz(D > 0)

    slope = []
    for r in Rt:
        z = interp_d(r, 1e-5)
        if not np.isfinite(z):
            slope.append(np.nan)
            continue
        zl, zh = z * 0.9, z * 1.1
        rl, rh = np.interp([np.log(zl), np.log(zh)], np.log(d0m),
                           np.log(np.maximum(r, 1e-300)))
        slope.append((rh - rl) / (np.log(zh) - np.log(zl)))
    slope = np.array(slope)

    fig = plt.figure(figsize=(7.5, 3.1))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.42)

    ax = fig.add_subplot(gs[0, 0])
    for n, (name, sel, sarc) in enumerate(per_fault):
        first = n == 0
        ax.plot(sarc, eq[sel], "-", color=C_Y03, lw=1.7,
                label="equal hazard, mean curve, $10^{5}$ a" if first else None)
        ax.plot(sarc, p85[sel], "-", color=C_V25, lw=1.5,
                label="85th percentile given rupture" if first else None)
        ax.plot(sarc, p50[sel], "--", color=C_V25, lw=1.4,
                label="median given rupture (eq. 9)" if first else None)
        ax.annotate(name, (sarc[len(sarc) // 2], eq[sel][len(sel) // 2]),
                    fontsize=6.2, ha="center", va="bottom",
                    xytext=(0, 4), textcoords="offset points")
    ax.set_yscale("log")
    ax.set_xlabel("distance along strike (km)")
    ax.set_ylabel("principal displacement (m)")
    ax.set_title("(a) design basis along the fault", loc="left")
    ax.legend(fontsize=6.2, frameon=False, loc="lower right")
    ax.grid(True, which="both", color="0.9", lw=0.4)

    ax2 = fig.add_subplot(gs[0, 1])
    for c in ccdf[::7]:
        ax2.plot(d0m, c, "-", color=C_V25, lw=0.7, alpha=0.45)
    ax2.axhline(0.15, color="0.15", lw=1.0, ls="--")
    ax2.axhline(0.50, color="0.15", lw=1.0, ls=":")
    ax2.text(1.5e-4, 0.17, "85th percentile", fontsize=6.2)
    ax2.text(1.5e-4, 0.52, "median", fontsize=6.2)
    ax2.set_xscale("log")
    ax2.set_xlim(1e-4, 30)
    ax2.set_ylim(0, 1)
    ax2.set_xlabel("displacement (m)")
    ax2.set_ylabel("$P(D>z\\,|\\,D>0)$")
    ax2.set_title("(b) given surface rupture", loc="left")
    ax2.grid(True, which="both", color="0.9", lw=0.4)

    save(fig, "fig_design_basis")

    print("=== design basis, on-fault principal (Abrahamson & Liou 2025) ===")
    print(f"  rate of surface rupture at the site  "
          f"{np.median(rup):.3e} /a   RP {1/np.median(rup):,.0f} a")
    print("  log-log hazard SLOPE at 1e5 a        "
          f"{np.nanmedian(slope):+.2f}  -> "
          + ("STEEP: the paper's alternatives do NOT apply"
             if np.nanmedian(slope) < -1 else
             "FLAT: the paper's alternatives apply"))
    print(f"  equal hazard, mean curve             "
          f"{np.nanmedian(eq):.2f} m")
    print(f"  85th percentile given rupture        "
          f"{np.nanmedian(p85):.2f} m")
    print(f"  median given rupture                 "
          f"{np.nanmedian(p50):.2f} m")


# ------------------------------------------------------------ figure 5
def fig_principal_fractiles(d0=None, branches=None):
    """Epistemic fractile maps of the on-fault principal displacement,
    with the branch-by-branch hazard curve at the on-fault site PF.

    Panel (c) shows where the map fractiles come from: every logic-tree
    branch as a grey curve, the MEAN hazard in black, and the 16th-84th
    band (+-1 sigma) shaded. The black line is the weighted mean of the
    branch RATES, which is the standard mean hazard curve and the same
    quantity the equal-hazard map value is read from - not the median of
    panel (a), which is a fractile. All three panels are PRINCIPAL
    only and all three fractiles are taken ACROSS BRANCHES, so the
    spread shown is epistemic; the aleatory variability is already
    integrated inside each individual grey curve.
    """
    res = principal_fractile_maps((0.84,))
    if res is None:
        return
    lons, lats, qmaps = res
    rows = [r for r in csv.DictReader(
        l for l in open(HERE / "out_map" / "aggregate"
                        / "displacement_map_mean.csv")
        if not l.startswith("#"))]
    mtr = np.array([r.get("is_trace", "0") == "1" for r in rows])

    panels = [("mean", qmaps["mean"]), (0.84, qmaps[0.84])]
    allv = np.concatenate([v[mtr][v[mtr] > 0] for _, v in panels])
    lo = np.floor(allv.min() * 10) / 10
    hi = np.ceil(allv.max() * 10) / 10
    norm = Normalize(vmin=lo, vmax=hi)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#f2f2f2")

    fig = plt.figure(figsize=(6.4, 6.4))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 1.0, 0.055],
                          height_ratios=[1.0, 0.82], wspace=0.16,
                          hspace=0.34)
    axes = [fig.add_subplot(gs[0, 0])]
    axes.append(fig.add_subplot(gs[0, 1], sharey=axes[0]))
    for n, (ax, (q, v)) in enumerate(zip(axes, panels)):
        _draw_colored_traces(ax, lons[mtr], lats[mtr], v[mtr], cmap, norm)
        thalo = [pe.withStroke(linewidth=1.6, foreground="white")]
        for name, xy in fault_traces():
            tx, ty = ((13.262, 42.905) if name == "MVFS"
                      else (13.06, 42.925))
            ax.text(tx, ty, name, color="0.15", fontsize=6.5,
                    fontweight="bold", zorder=6, path_effects=thalo)
        for i, (x, y) in SITE_XY.items():
            ax.scatter([x], [y], marker="^", s=30, zorder=5,
                       facecolor=SITE_COLORS[i], edgecolor="0.15",
                       linewidth=0.5)
            dy = 0.008 if SITES[i] == "SL" else -0.020
            ax.text(x + 0.005, y + dy, SITES[i], fontsize=6.5,
                    fontweight="bold", color="0.15", zorder=6,
                    path_effects=thalo)
        lbl = ("mean hazard" if q == "mean"
               else f"epistemic {int(q * 100)}th percentile")
        ax.set_title(f"({'ab'[n]}) {lbl}", loc="left")
        ax.set_xlabel("longitude (deg)")
        ax.set_xlim(13.05, 13.35)
        ax.set_ylim(42.65, 42.95)
        # sparse ticks: the default set collides across the panel seam
        ax.set_xticks([13.1, 13.2, 13.3])
    axes[0].set_ylabel("latitude (deg)")
    plt.setp(axes[1].get_yticklabels(), visible=False)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, cax=fig.add_subplot(gs[0, 2]))
    cb.set_label("principal $D$ at $10^{5}$ a (m)", fontsize=7,
                 labelpad=3)
    cb.ax.tick_params(labelsize=7)

    # (c) branch-by-branch principal hazard curve at PF
    if branches is not None:
        axc = fig.add_subplot(gs[1, :2])
        w = np.array([b["w"] for b in branches], dtype=float)
        allr = np.array([b["principal"][0] for b in branches])
        for y in allr:
            axc.plot(d0, np.maximum(y, 1e-16), "-", color="0.72", lw=0.45,
                     alpha=0.55, zorder=1)
        mean = (w[:, None] * allr).sum(0) / w.sum()
        q16 = np.zeros(len(d0)); q84 = np.zeros(len(d0))
        for j in range(len(d0)):
            o = np.argsort(allr[:, j])
            cw = np.cumsum(w[o]) / w.sum()
            col = allr[o, j]
            q16[j] = col[min(np.searchsorted(cw, 0.16), len(col) - 1)]
            q84[j] = col[min(np.searchsorted(cw, 0.84), len(col) - 1)]
        axc.fill_between(d0, np.maximum(q16, 1e-16), np.maximum(q84, 1e-16),
                         color="#0072B2", alpha=0.25, lw=0, zorder=2,
                         label="16th-84th (\u00b11$\\sigma$)")
        axc.plot(d0, np.maximum(mean, 1e-16), "-", color="black", lw=1.8,
                 zorder=3, label="mean hazard")
        axc.plot([], [], "-", color="0.72", lw=0.9,
                 label=f"{len(branches)} branches")
        axc.axhline(1e-5, color="0.3", lw=0.8, ls=":", zorder=1)
        axc.text(1.15e-2, 1.18e-5, "$10^{5}$ a", fontsize=6.5,
                 color="0.3")
        axc.set_xscale("log"); axc.set_yscale("log")
        axc.set_xlim(1e-2, 30); axc.set_ylim(1e-7, 1e-3)
        axc.set_xlabel("principal displacement (m)")
        axc.set_ylabel("annual rate of exceedance (a$^{-1}$)")
        axc.set_title("(c) epistemic spread at PF", loc="left")
        axc.legend(fontsize=6.2, frameon=False, loc="lower left")
        axc.grid(True, which="both", color="0.92", lw=0.4)
        d_at = lambda c: float(np.exp(np.interp(
            np.log(1e-5), np.log(np.maximum(c, 1e-300))[::-1],
            np.log(d0)[::-1])))
        print(f"  PF principal at 1e5 a: mean {d_at(mean):.2f} m, "
              f"16th {d_at(q16):.2f} m, 84th {d_at(q84):.2f} m")

    save(fig, "fig_principal_fractiles")
    for q, v in panels:
        t = v[mtr]
        print(f"  principal on-fault {str(q):<5}: median of sites "
              f"{np.median(t):.2f} m  [{t.min():.2f} - {t.max():.2f}]")


# ------------------------------------------------------------ figure 6
#: display names for the node ids that parse_branch_id returns
NODE_LABEL = {
    "PIZZA23": "Pizza et al. (2023)\nItalian normal",
    "Y03EC": "Youngs et al. (2003)\nExtensional Cordillera",
    "Y03NBR": "Youngs et al. (2003)\nN. Basin & Range",
    "Y03GB": "Youngs et al. (2003)\nGreat Basin",
    "Y03_AD": "Youngs et al. (2003)\naverage displacement",
    "Y03_MD": "Youngs et al. (2003)\nmaximum displacement",
    "LAV23": "Lavrentiadis &\nAbrahamson (2023)",
    "Y03": "Youngs et al. (2003)\n500 m cell",
    "V25": "Visini et al. (2025)\n500 m cell",
    "Y03_P85": "85th-percentile fit",
    "Y03_P95": "95th-percentile fit",
    "V25_WC94": "Wells & Coppersmith\n(1994) scaling",
    "V25_TB17": "Thingbaijam et al.\n(2017) scaling",
    "V25_LEO10": "Leonard (2010)\nscaling",
    "ACCURATE": "$\\sigma_r$ = 26.89 m\n(accurately located)",
}
#: which family colour each node belongs to
NODE_HUE = {"PIZZA23": "#009E73", "Y03EC": C_Y03, "Y03NBR": C_Y03,
            "Y03GB": C_Y03, "Y03_AD": C_Y03, "Y03_MD": C_Y03,
            "LAV23": "#CC79A7", "Y03": C_Y03, "V25": C_V25,
            "Y03_P85": C_Y03, "Y03_P95": C_Y03_95, "V25_WC94": C_V25,
            "V25_TB17": C_V25_TB, "V25_LEO10": "#E69F00",
            "ACCURATE": "0.35"}


def fig_logic_tree(branches):
    """The FDHA logic tree, drawn from the branches that were actually
    run so the figure cannot drift from the calculation.

    Node membership, level-to-level connectivity and weights are all
    recovered from the executed end branches. Testing that connectivity
    shows only ONE transition is a genuine pairing: distributed-rupture
    model to distributed-displacement model, where each model carries
    its own published displacement variants. Every other transition is a
    full cross product, so those are drawn as a product operator rather
    than as a bundle of lines that would carry no information.

    The tree is 4 x 3 x (2 + 3) x 1 = 60 end branches; weights are
    normalised and multiply along each path. Two points belong in the
    figure caption rather than on the page: the S_sr level is subsumed
    in the pairing term, so it does not appear separately in that count;
    and the 0.25 / 0.25 / 0.5 weighting of the principal-displacement
    level is what expresses the two displacement definitions (Youngs
    principal, LA23 sum-of-principal) being held in separate
    definition-consistent branch sets in the input file.
    """
    LEVELS = [("PSR", "Surface-rupture\nprobability", W_PSR),
              ("PFD", "Principal\ndisplacement", W_PFD),
              ("SSR", "Distributed-rupture\nprobability", W_SSR),
              ("SFD", "Distributed\ndisplacement", W_SFD),
              ("RS", "Rupture-location\nuncertainty", W_RS)]
    keys = [k for k, _t, _w in LEVELS]

    order = {}
    for k in keys:
        seen = []
        for b in branches:
            if b["node"][k] not in seen:
                seen.append(b["node"][k])
        order[k] = seen
    # the one real pairing: which S_fd nodes hang off which S_sr node
    pairing = {}
    for b in branches:
        pairing.setdefault(b["node"]["SSR"], [])
        if b["node"]["SFD"] not in pairing[b["node"]["SSR"]]:
            pairing[b["node"]["SSR"]].append(b["node"]["SFD"])

    BH, BG = 0.46, 0.46          # box height, vertical gap
    def stack(n):                # centred vertical slots for n boxes
        span = (n - 1) * (BH + BG)
        return [span / 2 - i * (BH + BG) for i in range(n)]

    fig, ax = plt.subplots(figsize=(8.4, 4.4))
    XW, GAP = 1.72, 0.86
    x = 0.0
    xpos = {}
    for k in keys:
        xpos[k] = x
        x += XW + (GAP if k != "SSR" else 0.52)

    ytop = max(stack(len(order["SFD"]))) + 0.62
    for li, (k, title, wtab) in enumerate(LEVELS):
        xs = xpos[k]
        ax.text(xs + XW / 2, ytop + 0.10, title, ha="center", va="bottom",
                fontsize=7.6, fontweight="bold", color="0.2",
                linespacing=1.3)
        slots = stack(len(order[k]))
        for n, y in zip(order[k], slots):
            hue = NODE_HUE.get(n, "0.4")
            ax.text(xs + XW / 2, y + 0.115, NODE_LABEL.get(n, n),
                    ha="center", va="center", fontsize=6.4,
                    linespacing=1.3, zorder=3, color="0.12")
            ax.text(xs + XW / 2, y - 0.235, f"w = {wtab[n]:g}",
                    ha="center", va="center", fontsize=6.4,
                    fontweight="bold", zorder=3,
                    color=hue if hue.startswith("#") else "0.3")
        # connector to the next level
        if li == len(LEVELS) - 1:
            continue
        nk = keys[li + 1]
        if k == "SSR":                       # the one real pairing
            pslots = dict(zip(order[k], slots))
            cslots = dict(zip(order[nk], stack(len(order[nk]))))
            for par, kids in pairing.items():
                for kid in kids:
                    xm = xs + XW + 0.26
                    ax.plot([xs + XW + 0.06, xm, xm, xpos[nk] + 0.12],
                            [pslots[par]] * 2 + [cslots[kid]] * 2, "-",
                            color="0.6", lw=0.8, zorder=1,
                            solid_joinstyle="round")
        else:
            ax.text(xs + XW + GAP / 2, 0.0, "$\\times$", ha="center",
                    va="center", fontsize=13, color="0.45", zorder=3)

    ax.set_xlim(-0.18, xpos["RS"] + XW + 0.18)
    ax.set_ylim(min(stack(len(order["SFD"]))) - 0.55, ytop + 0.72)
    ax.axis("off")
    save(fig, "fig_logic_tree")
    print(f"  logic tree: {len(branches)} end branches, "
          f"nodes per level {[len(order[k]) for k in keys]}")


def main():
    d0, branches = load_curves(HERE / "out_curve")
    total_w = sum(b["w"] for b in branches)
    assert abs(total_w - 1.0) < 1e-9, total_w
    print(f"loaded {len(branches)} branches, weight sum {total_w:.12f}")
    fig_logic_tree(branches)
    fig_ms_tunnel(d0, branches)
    fig_three_site_fans(d0, branches)
    fig_profile()
    fig_design_basis()
    fig_principal_fractiles(d0, branches)
    print("figures written to", FIGDIR)


if __name__ == "__main__":
    main()
