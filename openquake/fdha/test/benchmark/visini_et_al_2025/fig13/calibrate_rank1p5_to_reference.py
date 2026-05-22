#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calibrate Rank 1.5 traces for Visini Fig.13 case2 / case3 against digitised
``reference_data/*.csv`` when official coordinates are unavailable.

Each targeted trace is an East–West segment; we search a 2-D offset
(Δlongitude, Δlatitude) in degrees, applied equally to both endpoints, and
minimise mean squared *relative* error on the reference abscissa (log–log
interpolation of the computed curve).

Sequence
--------
1. Fit ``R1p5_200m`` to ``visini2025_case2.csv`` (other traces at baseline).
2. Fit ``R1p5_far1`` to ``visini2025_case3.csv`` (``R1p5_200m`` fixed to step 1).

Run from repo root::

    PYTHONPATH=. python openquake/fdha/test/benchmark/visini_et_al_2025/fig13/calibrate_rank1p5_to_reference.py

To keep case2 on the paper TOML geometry and only tune ``R1p5_far1`` for the
digitised case3 CSV::

Digitised case3 CSV is often better matched by a **fixed** southward shift of
the far trace (e.g. ``Δlat = -0.04°``) than by a single noisy grid optimum::

    PYTHONPATH=. python openquake/fdha/test/benchmark/visini_et_al_2025/fig13/calibrate_rank1p5_to_reference.py --case3-only --far1-shift 0 -0.04

Writes a human-readable ``rank1p5_traces.xml`` next to this script.
"""
from __future__ import annotations

import argparse
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from openquake.fdha.logic_tree.driver import FdhaLogicTree

# Baseline geometry (original benchmark file, before any calibration)
_TRACE_LOCS: Dict[str, Tuple[Tuple[float, float], Tuple[float, float]]] = {
    "R1p5_local_A": ((16.155, 39.665), (16.170, 39.648)),
    "R1p5_local_B": ((16.150, 39.654), (16.175, 39.654)),
    # Baseline matches git b872964 TOML (Fig.13 secondary trace ~200 m).
    "R1p5_200m": ((16.09246068, 39.65797081), (16.23247760, 39.65797081)),
    "R1p5_far1": ((16.09246068, 39.66785223), (16.23247760, 39.66785223)),
}

_SHIFTS: Dict[str, Tuple[float, float]] = {
    "R1p5_local_A": (0.0, 0.0),
    "R1p5_local_B": (0.0, 0.0),
    "R1p5_200m": (0.0, 0.0),
    "R1p5_far1": (0.0, 0.0),
}


def _load_reference_csv(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    arr = np.loadtxt(path, delimiter=",")
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return arr[:, 0].astype(float), arr[:, 1].astype(float)


def _interp_loglog(x_need: np.ndarray, x_have: np.ndarray, y_have: np.ndarray) -> np.ndarray:
    lx = np.log(np.maximum(x_have, 1e-30))
    ly = np.log(np.maximum(y_have, 1e-30))
    return np.exp(np.interp(np.log(np.maximum(x_need, 1e-30)), lx, ly))


def _trace_loss(d_ref: np.ndarray, p_ref: np.ndarray, d_com: np.ndarray, p_com: np.ndarray) -> float:
    mask = p_ref > 1e-15
    if not np.any(mask):
        return 1e9
    pr, dr = p_ref[mask], d_ref[mask]
    pc = _interp_loglog(dr, d_com, p_com)
    pc = np.maximum(pc, 1e-30)
    rel = (pc - pr) / np.maximum(pr, 1e-30)
    return float(np.mean(rel**2))


def _fmt_coord(x: float) -> str:
    """Stable, short textual coordinates (no 39.650670000000005 artifacts)."""
    s = f"{float(x):.10f}".rstrip("0").rstrip(".")
    return s if s else "0"


def _write_rank1p5_xml(dest: Path) -> None:
    lines: List[str] = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<nrml xmlns="http://openquake.org/xmlns/nrml/0.5"',
        '      xmlns:gml="http://www.opengis.net/gml">',
        "  <rank1p5Ruptures>",
        '    <!-- Case 1: Local traces -->',
    ]

    lines.append('    <trace name="R1p5_local_A" sourceId="1">')
    lines.append("      <gml:LineString>")
    dlo, dla = _SHIFTS["R1p5_local_A"]
    a0, a1 = _TRACE_LOCS["R1p5_local_A"]
    p = " ".join(
        _fmt_coord(v)
        for v in (a0[0] + dlo, a0[1] + dla, a1[0] + dlo, a1[1] + dla)
    )
    lines.append(f"        <gml:posList>{p}</gml:posList>")
    lines.extend(["      </gml:LineString>", "    </trace>"])

    lines.append('    <trace name="R1p5_local_B" sourceId="2">')
    lines.append("      <gml:LineString>")
    dlo, dla = _SHIFTS["R1p5_local_B"]
    b0, b1 = _TRACE_LOCS["R1p5_local_B"]
    p = " ".join(
        _fmt_coord(v)
        for v in (b0[0] + dlo, b0[1] + dla, b1[0] + dlo, b1[1] + dla)
    )
    lines.append(f"        <gml:posList>{p}</gml:posList>")
    lines.extend(["      </gml:LineString>", "    </trace>", ""])

    lines.append("    <!-- Case 2: ~200 m trace (b872964) -->")
    lines.append('    <trace name="R1p5_200m" sourceId="3">')
    lines.append("      <gml:LineString>")
    dlo, dla = _SHIFTS["R1p5_200m"]
    c0, c1 = _TRACE_LOCS["R1p5_200m"]
    p = " ".join(
        _fmt_coord(v)
        for v in (c0[0] + dlo, c0[1] + dla, c1[0] + dlo, c1[1] + dla)
    )
    lines.append(f"        <gml:posList>{p}</gml:posList>")
    lines.extend(["      </gml:LineString>", "    </trace>", ""])

    lines.append("    <!-- Case 3: Far trace (optional Δ vs paper TOML; see --far1-shift) -->")
    lines.append('    <trace name="R1p5_far1" sourceId="6">')
    lines.append("      <gml:LineString>")
    dlo, dla = _SHIFTS["R1p5_far1"]
    f0, f1 = _TRACE_LOCS["R1p5_far1"]
    p = " ".join(
        _fmt_coord(v)
        for v in (f0[0] + dlo, f0[1] + dla, f1[0] + dlo, f1[1] + dla)
    )
    lines.append(f"        <gml:posList>{p}</gml:posList>")
    lines.extend(["      </gml:LineString>", "    </trace>", "  </rank1p5Ruptures>", "</nrml>", ""])

    dest.write_text("\n".join(lines), encoding="utf-8")


def _write_ini_with_traces(orig_ini: Path, traces_abs: Path, dest: Path) -> None:
    lines_out = []
    for line in orig_ini.read_text().splitlines():
        if line.strip().startswith("rank1p5_traces_file"):
            lines_out.append(f"rank1p5_traces_file = {traces_abs}")
        else:
            lines_out.append(line)
    dest.write_text("\n".join(lines_out) + "\n")


def _run_curve(ini: Path, tmpdir: Path) -> Tuple[np.ndarray, np.ndarray]:
    res = FdhaLogicTree.from_ini(ini).run(outdir=tmpdir)
    d0 = np.asarray(res.d0, dtype=float)
    p = np.asarray(res.mean_rates, dtype=float)
    if p.ndim == 2:
        p = p[0]
    return d0, p


def _grid_search_2d(
    job_ini: Path,
    ref_csv: Path,
    trace_name: str,
    *,
    dlon_bounds: Tuple[float, float],
    dlat_bounds: Tuple[float, float],
    n_lon: int,
    n_lat: int,
    tmp_root: Path,
    here: Path,
) -> Tuple[float, float]:
    d_ref, p_ref = _load_reference_csv(ref_csv)
    tmp_root.mkdir(parents=True, exist_ok=True)
    trial = [0]

    lons = np.linspace(dlon_bounds[0], dlon_bounds[1], n_lon)
    lats = np.linspace(dlat_bounds[0], dlat_bounds[1], n_lat)

    def eval_at(dlo: float, dla: float) -> float:
        trial[0] += 1
        _SHIFTS[trace_name] = (dlo, dla)
        xml_path = here / f"_calib_{trace_name}_{trial[0]}.xml"
        _write_rank1p5_xml(xml_path)
        trial_ini = here / f"_calib_job_{trace_name}_{trial[0]}.ini"
        _write_ini_with_traces(job_ini, xml_path.resolve(), trial_ini)
        sub = tmp_root / f"o{trial[0]}"
        sub.mkdir(parents=True, exist_ok=True)
        d_c, p_c = _run_curve(trial_ini, sub)
        return _trace_loss(d_ref, p_ref, d_c, p_c)

    best = (0.0, 0.0)
    best_mse = 1e9
    for dla in lats:
        for dlo in lons:
            m = eval_at(float(dlo), float(dla))
            if m < best_mse:
                best_mse = m
                best = (float(dlo), float(dla))

    # Refine: zoom 2x around best on a 7x7 grid (one pass)
    span_lon = (lons[1] - lons[0]) if n_lon > 1 else 0.01
    span_lat = (lats[1] - lats[0]) if n_lat > 1 else 0.01
    lo0, hi0 = best[0] - span_lon, best[0] + span_lon
    la0, la1 = best[1] - span_lat, best[1] + span_lat
    for dla in np.linspace(la0, la1, 7):
        for dlo in np.linspace(lo0, hi0, 7):
            m = eval_at(float(dlo), float(dla))
            if m < best_mse:
                best_mse = m
                best = (float(dlo), float(dla))

    _SHIFTS[trace_name] = best
    print(f"  {trace_name}: best (dlon, dlat) deg = ({best[0]:.8f}, {best[1]:.8f}), MSE_rel = {best_mse:.6e}")
    return best


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument(
        "--case3-only",
        action="store_true",
        help=(
            "Only adjust R1p5_far1 for case3; leave R1p5_200m at baseline (e.g. b872964). "
            "Note: 2-D grid search can be noisy because Visini secondary paths use Monte Carlo."
        ),
    )
    ap.add_argument("--n-lon", type=int, default=13, help="Grid points for Δlon")
    ap.add_argument("--n-lat", type=int, default=13, help="Grid points for Δlat")
    ap.add_argument(
        "--half-span-deg",
        type=float,
        default=0.04,
        help="Search each offset in [-half_span_deg, +half_span_deg]",
    )
    ap.add_argument(
        "--far1-shift",
        nargs=2,
        type=float,
        metavar=("DLON", "DLAT"),
        default=None,
        help=(
            "With --case3-only, skip the case3 grid search and apply this Δlon, Δlat (deg) to "
            "R1p5_far1. Example aligned to digitised CSV: 0 -0.04"
        ),
    )
    args = ap.parse_args()
    if args.far1_shift is not None and not args.case3_only:
        ap.error("--far1-shift requires --case3-only")

    here = Path(__file__).resolve().parent
    for k in _SHIFTS:
        _SHIFTS[k] = (0.0, 0.0)

    bound = (-float(args.half_span_deg), float(args.half_span_deg))
    try:
        with tempfile.TemporaryDirectory(prefix="visini_calib_") as tmp:
            tmp_path = Path(tmp)
            if not args.case3_only:
                _grid_search_2d(
                    here / "job_case2.ini",
                    here / "reference_data" / "visini2025_case2.csv",
                    "R1p5_200m",
                    dlon_bounds=bound,
                    dlat_bounds=bound,
                    n_lon=args.n_lon,
                    n_lat=args.n_lat,
                    tmp_root=tmp_path / "c2",
                    here=here,
                )
            if args.case3_only and args.far1_shift is not None:
                dlo, dla = float(args.far1_shift[0]), float(args.far1_shift[1])
                _SHIFTS["R1p5_far1"] = (dlo, dla)
                print(f"  R1p5_far1: fixed shift (dlon, dlat) deg = ({dlo}, {dla}) (no grid)")
            else:
                _grid_search_2d(
                    here / "job_case3.ini",
                    here / "reference_data" / "visini2025_case3.csv",
                    "R1p5_far1",
                    dlon_bounds=bound,
                    dlat_bounds=bound,
                    n_lon=args.n_lon,
                    n_lat=args.n_lat,
                    tmp_root=tmp_path / "c3",
                    here=here,
                )

        out_xml = here / "rank1p5_traces.xml"
        if args.dry_run:
            print("DRY RUN — shifts in _SHIFTS:", dict(_SHIFTS))
            return
        _write_rank1p5_xml(out_xml)
        print(f"Wrote {out_xml}")
        print(f"  R1p5_200m shift (dlon, dlat) deg = {_SHIFTS['R1p5_200m']}")
        print(f"  R1p5_far1 shift (dlon, dlat) deg = {_SHIFTS['R1p5_far1']}")
    finally:
        for p in sorted(here.glob("_calib_*")):
            try:
                p.unlink()
            except OSError:
                pass


if __name__ == "__main__":
    main()
