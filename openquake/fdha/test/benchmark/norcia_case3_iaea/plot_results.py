#!/usr/bin/env python
"""Generate SRL manuscript figures for the Norcia Case 3 IAEA benchmark."""

from __future__ import annotations

import csv
import json
import os
import xml.etree.ElementTree as ET
from pathlib import Path

os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

from openquake.fdha.calc.calculators import (
    BaseFaultRuptureCalculator,
    FaultRuptureProbabilityCalculator,
)
from openquake.fdha.calc.config_loader import load_config
from openquake.fdha.calc.hazard import calculate_fdha_hazard
from openquake.fdha.calc.utils.interpolation import get_map_from_curves
from openquake.fdha.logic_tree.driver import FdhaLogicTree
from openquake.fdha.logic_tree.site_builder import build_hazard_map_sites


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
FIGURES = ROOT / "Figures"
SOURCE_MODEL = HERE / "source_model_norcia_case3.xml"
CURVE_INI = HERE / "job_norcia_case3_iaea_curve.ini"
MAP_INI = HERE / "job_norcia_case3_iaea_map.ini"
OUT_CURVE = HERE / "out_curve"
OUT_MAP = HERE / "out_map"

APPLICATION_SITES = {
    "Mount Serra": (13.188, 42.749),
    "San Lorenzo": (13.212, 42.853),
    "Principal": (13.278, 42.767),
}
SITE_LABELS = {
    "Mount Serra": "MS",
    "San Lorenzo": "SL",
    "Principal": "PF",
}
SITE_COLORS = {
    "Mount Serra": "#0072B2",
    "San Lorenzo": "#D55E00",
    "Principal": "#111111",
}
CURVE_SITE_ORDER = ("Principal", "Mount Serra", "San Lorenzo")
SOURCE_ID = "MVFS,NFS"


def _curve_output_matches_sites() -> bool:
    path = OUT_CURVE / "aggregate_hazard.csv"
    if not path.is_file():
        return False
    try:
        manifest = json.loads((OUT_CURVE / "manifest.json").read_text())
        if any(
            branch.get("fdha_source_id") != SOURCE_ID
            for branch in manifest.get("branches", [])
        ):
            return False
        with path.open(newline="") as f:
            reader = csv.DictReader(f)
            if not reader.fieldnames or "site_id" not in reader.fieldnames:
                return False
            seen: dict[int, tuple[float, float]] = {}
            for row in reader:
                sid = int(row["site_id"])
                if sid not in seen:
                    seen[sid] = (float(row["lon"]), float(row["lat"]))
        if len(seen) != len(CURVE_SITE_ORDER):
            return False
        for sid, label in enumerate(CURVE_SITE_ORDER):
            expected_lon, expected_lat = APPLICATION_SITES[label]
            got_lon, got_lat = seen[sid]
            if not (
                np.isclose(got_lon, expected_lon)
                and np.isclose(got_lat, expected_lat)
            ):
                return False
    except (KeyError, OSError, TypeError, ValueError):
        return False
    return True


def ensure_outputs(run_map: bool = True) -> None:
    if not _curve_output_matches_sites():
        FdhaLogicTree.from_ini(CURVE_INI).run(outdir=OUT_CURVE)
    if run_map and not (OUT_MAP / "aggregate" / "rates_mean.h5").is_file():
        FdhaLogicTree.from_ini(MAP_INI).run(outdir=OUT_MAP)


def read_aggregate_curves(outdir: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    rows_by_sid: dict[int, list[tuple[float, float]]] = {
        sid: [] for sid in range(len(CURVE_SITE_ORDER))
    }
    with (outdir / "aggregate_hazard.csv").open(newline="") as f:
        for row in csv.DictReader(f):
            sid = int(row["site_id"])
            if sid in rows_by_sid:
                rows_by_sid[sid].append((float(row["D0"]), float(row["mean"])))
    curves = {}
    for sid, label in enumerate(CURVE_SITE_ORDER):
        arr = np.asarray(rows_by_sid[sid], dtype=float)
        if arr.size == 0:
            raise ValueError(f"Missing curve rows for site {sid} ({label})")
        curves[label] = (arr[:, 0], arr[:, 1])
    return curves


def branch_config_paths(outdir: Path) -> list[Path]:
    return sorted((outdir / "source_model_branches").glob("*/branch_configs/branch_*.ini"))


def branch_weights(outdir: Path) -> np.ndarray:
    manifest = json.loads((outdir / "manifest.json").read_text())
    weights = np.array([float(b["combined_branch_weight"]) for b in manifest["branches"]])
    return weights / weights.sum()


def branch_labels(outdir: Path) -> list[str]:
    manifest = json.loads((outdir / "manifest.json").read_text())
    labels = []
    for branch in manifest["branches"]:
        secondary = branch["fdha_models"].get("secondary_surf_rup", "")
        labels.append("V24" if secondary.startswith("Visini2025") else "Y03")
    return labels


def curve_components() -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    configs = branch_config_paths(OUT_CURVE)
    weights = branch_weights(OUT_CURVE)
    labels = branch_labels(OUT_CURVE)
    total = principal = distributed = None
    branch_curves: dict[str, np.ndarray] = {}
    d0 = None
    for cfg_path, weight, label in zip(configs, weights, labels):
        calc = FaultRuptureProbabilityCalculator(
            str(cfg_path),
            [str(SOURCE_MODEL)],
            rupture_mesh_spacing=0.5,
            width_of_mfd_bin=0.1,
        )
        res = calc.run()
        d0 = np.asarray(res["imls"], dtype=float)
        p = np.asarray(res["rate_principal"], dtype=float)[0]
        s = np.asarray(res["rate_distributed"], dtype=float)[0]
        t = p + s
        branch_curves[label] = t
        if total is None:
            total = np.zeros_like(t)
            principal = np.zeros_like(p)
            distributed = np.zeros_like(s)
        total += weight * t
        principal += weight * p
        distributed += weight * s
    return d0, principal, distributed, branch_curves


def load_fault_traces() -> dict[str, np.ndarray]:
    ns = {
        "nrml": "http://openquake.org/xmlns/nrml/0.5",
        "gml": "http://www.opengis.net/gml",
    }
    root = ET.parse(SOURCE_MODEL).getroot()
    traces = {}
    for src in root.findall(".//*[@id]"):
        name = src.get("name") or src.get("id")
        pos = src.find(".//gml:posList", ns)
        if pos is None or not pos.text:
            continue
        vals = [float(v) for v in pos.text.split()]
        traces[name] = np.asarray(list(zip(vals[::2], vals[1::2])), dtype=float)
    return traces


def make_curve_figure() -> Path:
    ensure_outputs(run_map=False)
    curves = read_aggregate_curves(OUT_CURVE)

    FIGURES.mkdir(parents=True, exist_ok=True)
    out = FIGURES / "norcia_case3_curve_vs_iaea.pdf"

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    pf_d0, pf_mean = curves["Principal"]
    ms_d0, ms_mean = curves["Mount Serra"]
    sl_d0, sl_mean = curves["San Lorenzo"]
    ax.loglog(pf_d0, pf_mean, color=SITE_COLORS["Principal"], lw=2.4, label="PF principal-site")
    ax.loglog(ms_d0, ms_mean, color=SITE_COLORS["Mount Serra"], lw=2.2, label="MS distributed-site")
    ax.loglog(sl_d0, sl_mean, color=SITE_COLORS["San Lorenzo"], lw=2.2, label="SL distributed-site")

    ax.set_xlabel("Fault displacement (m)")
    ax.set_ylabel("Annual frequency of exceedance (1/yr)")
    ax.set_xlim(8e-4, 10.5)
    ax.set_ylim(1e-10, 6e-4)
    ax.grid(True, which="both", ls=":", color="#bbbbbb", alpha=0.55)
    ax.legend(loc="lower left", fontsize=8.8, frameon=True, framealpha=0.92)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def _map_sites_and_config():
    cfg = load_config(MAP_INI)
    geom = cfg["geometry"]
    erf = cfg["erf"]
    sites = build_hazard_map_sites(
        region=geom["region"],
        spacing=float(geom["region_grid_spacing"]),
        max_distance_km=float(geom.get("max_distance_km", 50.0)),
        vs30=cfg.get("site_location", {}).get("vs30"),
        source_model_paths=[str(SOURCE_MODEL)],
        hdf5path="",
        rupture_mesh_spacing=float(erf.get("rupture_mesh_spacing", 0.5)),
        width_of_mfd_bin=float(erf.get("width_of_mfd_bin", 0.1)),
    )
    return cfg, sites


def map_components() -> dict[str, np.ndarray]:
    ensure_outputs(run_map=True)
    cache = OUT_MAP / "component_maps.npz"
    if cache.is_file():
        data = np.load(cache)
        source_id = str(data["source_id"][0]) if "source_id" in data.files else ""
        if source_id == SOURCE_ID:
            return {k: data[k] for k in data.files if k != "source_id"}

    cfg, sites = _map_sites_and_config()
    target_displ = np.asarray(cfg["parameters"]["target_displacement"], dtype=float)
    return_period = float(cfg["calculation"].get("return_period", 100000.0))
    target_rate = 1.0 / return_period

    site_lons = np.array([s.location.longitude for s in sites.combined_sitecol], dtype=float)
    site_lats = np.array([s.location.latitude for s in sites.combined_sitecol], dtype=float)
    is_trace = np.zeros(len(site_lons), dtype=bool)
    is_trace[sites.n_active_grid:] = True

    weights = branch_weights(OUT_MAP)
    principal_rates = np.zeros((len(site_lons), len(target_displ)), dtype=float)
    distributed_rates = np.zeros_like(principal_rates)
    for cfg_path, weight in zip(branch_config_paths(OUT_MAP), weights):
        calc = BaseFaultRuptureCalculator(
            str(cfg_path),
            [str(SOURCE_MODEL)],
            hdf5path="",
            rupture_mesh_spacing=0.5,
            width_of_mfd_bin=0.1,
        )
        res = calculate_fdha_hazard(calc, sites.combined_sitecol, show_progress=False)
        principal_rates += weight * np.asarray(res["rate_principal"], dtype=float)
        distributed_rates += weight * np.asarray(res["rate_distributed"], dtype=float)

    total_rates = principal_rates + distributed_rates
    arrays = {
        "lon": site_lons,
        "lat": site_lats,
        "is_trace": is_trace,
        "total": get_map_from_curves(target_displ, total_rates, target_rate),
        "principal": get_map_from_curves(target_displ, principal_rates, target_rate),
        "distributed": get_map_from_curves(target_displ, distributed_rates, target_rate),
        "source_id": np.array([SOURCE_ID]),
    }
    np.savez(cache, **arrays)
    return arrays


def _plot_panel(ax, arrays: dict[str, np.ndarray], key: str, title: str | None, norm, cmap) -> None:
    lon = arrays["lon"]
    lat = arrays["lat"]
    values = arrays[key]
    is_trace = arrays["is_trace"].astype(bool)
    grid = ~is_trace
    plot_values = np.ma.masked_less_equal(values, 0.0)
    ax.scatter(lon[grid], lat[grid], c=plot_values[grid], s=8, marker="s", cmap=cmap, norm=norm, linewidths=0)
    ax.scatter(lon[is_trace], lat[is_trace], c=plot_values[is_trace], s=12, marker="o", cmap=cmap, norm=norm, linewidths=0)
    for name, trace in load_fault_traces().items():
        color = "#D62728" if name == "MVFS" else "#7F3C8D"
        ax.plot(trace[:, 0], trace[:, 1], color=color, lw=0.65, ls="--", solid_capstyle="round")
        if name == "MVFS":
            ax.text(13.245, 42.895, "MVFS", color=color, fontsize=7.0, weight="bold")
        else:
            ax.text(13.075, 42.842, "NFS", color=color, fontsize=7.0, weight="bold")
    for label, (x, y) in APPLICATION_SITES.items():
        ax.scatter(
            [x],
            [y],
            marker="^",
            s=48,
            facecolor=SITE_COLORS[label],
            edgecolor="#222222",
            linewidth=0.6,
            zorder=5,
        )
        ax.text(
            x + 0.004,
            y + (0.008 if label == "San Lorenzo" else -0.014),
            SITE_LABELS[label],
            fontsize=7.5,
            weight="bold",
            color="#222222",
            zorder=6,
        )
    if title:
        ax.set_title(title, fontsize=10)
    ax.set_xlim(13.05, 13.35)
    ax.set_ylim(42.65, 42.95)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Longitude")
    ax.grid(True, color="#d0d0d0", lw=0.4, alpha=0.5)


def make_map_figure() -> Path:
    arrays = map_components()
    FIGURES.mkdir(parents=True, exist_ok=True)
    out = FIGURES / "norcia_case3_map.pdf"

    positive = arrays["total"][arrays["total"] > 0]
    vmax = max(float(np.nanmax(positive)), 1e-2) if positive.size else 1.0
    norm = mcolors.LogNorm(vmin=1e-4, vmax=vmax)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#f3f3f3")

    fig = plt.figure(figsize=(5.4, 4.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.055], wspace=0.08)
    ax = fig.add_subplot(gs[0, 0])
    _plot_panel(ax, arrays, "total", None, norm, cmap)
    ax.set_ylabel("Latitude")

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cax = fig.add_subplot(gs[0, 1])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label("Displacement at 1e5 yr return period (m)")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> int:
    curve = make_curve_figure()
    hazard_map = make_map_figure()
    print(f"Saved {curve}")
    print(f"Saved {hazard_map}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
