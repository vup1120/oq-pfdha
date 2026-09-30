#!/usr/bin/env python
"""Build the fault-parallel strip site set for the unrolled corridor map.

The principal corridor is only ~108 m wide (the +-2 sigma truncation of
the Petersen rupture-location weight), so a regional grid fine enough to
resolve it is unaffordable. Instead this samples a narrow band that
follows each trace and parameterises every site by its position ALONG
the fault (s) and ACROSS it (x), which is the natural coordinate system
for the corridor and lets the result be plotted unrolled.

Why not simply set a small [geometry].max_distance_km on a map job:
that parameter is overloaded. `build_hazard_map_sites` uses it to keep
grid sites near the trace, but `calc.hazard` ALSO passes it to
FDHAContextMaker as the rupture integration distance, so shrinking it to
a few hundred metres would silently truncate the integration. Supplying
the band explicitly via [geometry].sites_csv keeps max_distance_km at
its proper value.

Writes:
  strip_sites.csv  lon,lat            -> consumed by job_strip.ini
  strip_meta.csv   fault,s_km,x_m,... -> consumed by make_figures.py
"""

from __future__ import annotations

import csv
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

ALONG_STEP_M = 200.0    # along-strike sampling
ACROSS_STEP_M = 10.0    # across-strike sampling (corridor edge is at 54 m)
ACROSS_MAX_M = 150.0    # half-width of the band

EARTH_KM_PER_DEG = 111.195


def traces():
    """(name, Nx2 lon/lat) for each fault, deduplicated by name."""
    gml = "{http://www.opengis.net/gml}"
    root = ET.parse(HERE / "source_model_norcia_case3.xml").getroot()
    out = {}
    for src in root.iter():
        # only real source elements: parent containers (sourceModel, the
        # nrml root) also carry a name and would match .//posList
        if not src.tag.endswith("Source"):
            continue
        name = src.get("name") or src.get("id")
        pls = src.findall(f".//{gml}posList")
        if not name or not pls or name in out:
            continue
        xy = np.array(pls[0].text.split(), dtype=float).reshape(-1, 2)
        if len(xy) >= 2:
            out[name] = xy
    return out


def resample(xy, lat0, step_m):
    """Resample a lon/lat polyline at fixed arc length; return the
    resampled points, their cumulative distance (m) and unit tangents,
    all in the local east/north km frame."""
    kx = EARTH_KM_PER_DEG * np.cos(np.radians(lat0))
    ky = EARTH_KM_PER_DEG
    p = np.column_stack((xy[:, 0] * kx, xy[:, 1] * ky))  # km
    seg = np.diff(p, axis=0)
    seglen = np.hypot(seg[:, 0], seg[:, 1])
    cum = np.concatenate(([0.0], np.cumsum(seglen)))
    total_m = cum[-1] * 1000.0
    s_m = np.arange(0.0, total_m + step_m, step_m)
    s_km = s_m / 1000.0
    x = np.interp(s_km, cum, p[:, 0])
    y = np.interp(s_km, cum, p[:, 1])
    pts = np.column_stack((x, y))
    # tangent by central difference on the resampled line
    tan = np.gradient(pts, axis=0)
    tan /= np.linalg.norm(tan, axis=1)[:, None]
    return pts, s_m, tan, (kx, ky)


def main():
    sites, meta = [], []
    for name, xy in traces().items():
        lat0 = float(xy[:, 1].mean())
        pts, s_m, tan, (kx, ky) = resample(xy, lat0, ALONG_STEP_M)
        # unit normal; orient so +x points EAST (footwall side - both
        # Norcia faults dip SW, so the hanging wall is the -x side,
        # matching the sign convention of the fig_profile transect)
        nrm = np.column_stack((-tan[:, 1], tan[:, 0]))
        flip = nrm[:, 0] < 0
        nrm[flip] *= -1.0

        offs = np.arange(-ACROSS_MAX_M, ACROSS_MAX_M + ACROSS_STEP_M,
                         ACROSS_STEP_M)
        for i in range(len(pts)):
            for xo in offs:
                q = pts[i] + nrm[i] * (xo / 1000.0)   # km frame
                lon, lat = q[0] / kx, q[1] / ky
                sites.append((lon, lat))
                meta.append((name, s_m[i] / 1000.0, xo, lon, lat))

    with open(HERE / "strip_sites.csv", "w", newline="") as fh:
        fh.write("lon,lat\n")
        for lon, lat in sites:
            fh.write(f"{lon:.6f},{lat:.6f}\n")

    with open(HERE / "strip_meta.csv", "w", newline="") as fh:
        fh.write("site_id,fault,s_km,x_m,lon,lat\n")
        for i, (name, s, xo, lon, lat) in enumerate(meta):
            fh.write(f"{i},{name},{s:.4f},{xo:.1f},{lon:.6f},{lat:.6f}\n")

    per_fault = {}
    for name, s, xo, _, _ in meta:
        per_fault.setdefault(name, set()).add(round(s, 4))
    print(f"{len(sites)} strip sites "
          f"(along {ALONG_STEP_M:.0f} m, across {ACROSS_STEP_M:.0f} m "
          f"over +-{ACROSS_MAX_M:.0f} m)")
    for name, ss in per_fault.items():
        print(f"  {name}: {len(ss)} along-strike stations, "
              f"{max(ss):.1f} km long")
    print(f"  projected 45-branch run time: "
          f"{len(sites) * 6.9 * 45 / 1000 / 3600:.2f} h")


if __name__ == "__main__":
    main()
