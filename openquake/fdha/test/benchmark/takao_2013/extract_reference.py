#!/usr/bin/env python
"""Extract the Fig. 10 / Fig. 11 curves of Takao et al. (2013) from the paper PDF.

The paper's figures are vector graphics (Excel charts): solid curves are
stroked polyline segments, dashed/dotted curves are sequences of small
*filled* rectangles. This script converts the figure page to SVG
(``pdftocairo -svg``), parses every path, transforms it to page
coordinates, calibrates each chart's plot-area rectangle against its known
axis ranges, and writes one CSV of extracted data points per figure into
``reference/``:

    fig10_points.csv   x = AD, D/AD or D;      y = probability density
    fig11a_points.csv  x = displacement (m);   y = annual rate (case (a))
    fig11b_points.csv  x = displacement (m);   y = annual rate (case (b))

Each row carries the ink color (black/gray, which in Fig. 11(a)
distinguishes the R = 3,000 yr curves (gray) from R = 30,000 yr (black))
and the source element type. Curve identity is NOT resolved here — the
points are matched to computed curves by the reproduce_* scripts.

Axis ranges (hard-coded from the figure axis labels):

    Fig 10   x log10: 1e-3 .. 1e3     y linear: 0 .. 1.2e-2
    Fig 11a  x log10: 1e-2 .. 1e1     y log10: 1e-10 .. 1e-3
    Fig 11b  x log10: 1e-2 .. 1e1     y log10: 1e-10 .. 1e-3

Usage::

    python extract_reference.py /path/to/Takao_et_al_2013.pdf

(The PDF is not committed; the extracted CSVs in ``reference/`` are.)
"""

from __future__ import annotations

import csv
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PAGE = 13  # 図10 / 図11

# Known axis ranges: (x0, x1, xlog, y0, y1, ylog)
AXES = {
    "fig10": (1e-3, 1e3, True, 0.0, 1.2e-2, False),
    "fig11a": (1e-2, 1e1, True, 1e-10, 1e-3, True),
    "fig11b": (1e-2, 1e1, True, 1e-10, 1e-3, True),
}


def parse_matrix(tr: str):
    m = re.search(r"matrix\(([^)]*)\)", tr or "")
    if not m:
        return (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)
    return tuple(float(v) for v in m.group(1).split(","))


def apply_matrix(mat, x, y):
    a, b, c, d, e, f = mat
    return a * x + c * y + e, b * x + d * y + f


def path_points(d: str, mat):
    """All vertex coordinates of an M/L/Z path, transformed to page space."""
    pts = []
    for cmd, coords in re.findall(r"([MLZ])\s*((?:-?[\d.e]+\s+-?[\d.e]+\s*)*)", d):
        vals = [float(v) for v in coords.split()]
        for i in range(0, len(vals) - 1, 2):
            pts.append(apply_matrix(mat, vals[i], vals[i + 1]))
    return pts


def extract(pdf: Path):
    with tempfile.TemporaryDirectory() as tmp:
        svg_path = Path(tmp) / "page.svg"
        subprocess.run(
            ["pdftocairo", "-svg", "-f", str(PAGE), "-l", str(PAGE),
             str(pdf), str(svg_path)],
            check=True,
        )
        body = svg_path.read_text().split("</defs>")[1]

    raw_paths = re.findall(r'<path ([^>]*?)d="([^"]*)"([^>]*?)/>', body)
    elements = []  # (kind, color, points)
    plot_rects, legend_rects = [], []
    for pre, d, post in raw_paths:
        attrs = dict(re.findall(r'([\w-]+)="([^"]*)"', pre + post))
        mat = parse_matrix(attrs.get("transform", ""))
        pts = path_points(d, mat)
        if not pts:
            continue
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        w, h = max(xs) - min(xs), max(ys) - min(ys)
        fill = attrs.get("fill", "none")
        stroke = attrs.get("stroke", "none")
        sw = attrs.get("stroke-width", "")

        if fill.startswith("rgb(100%"):
            # White rectangles: chart borders, plot areas, legend boxes.
            # Plot areas are the untransformed inner rectangles (> 100 pt
            # wide); legend boxes are the smaller ones.
            if w > 100 and h > 100:
                plot_rects.append((min(xs), min(ys), max(xs), max(ys)))
            elif w > 20 and h > 10:
                legend_rects.append((min(xs), min(ys), max(xs), max(ys)))
        elif fill.startswith(("rgb(0%", "rgb(50%")) and "stroke" not in pre + post:
            # Filled mini-rectangles: dashes/dots of the dashed curves.
            # The element size distinguishes dot patterns (< ~1 pt) from
            # dash patterns (several pt).
            if w < 12 and h < 12:
                color = "black" if fill.startswith("rgb(0%") else "gray"
                cx = (min(xs) + max(xs)) / 2
                cy = (min(ys) + max(ys)) / 2
                elements.append(("dash", color, [(cx, cy)], max(w, h)))
        elif sw in ("0.777", "0.775") and stroke.startswith(("rgb(0%", "rgb(50%")):
            # Stroked segments: solid curves (and legend swatch lines)
            color = "black" if stroke.startswith("rgb(0%") else "gray"
            elements.append(("stroke", color, pts, max(w, h)))

    # The chart borders (transformed rects) are also >100x100; keep only
    # rects whose aspect/size matches the inner plot areas: choose, per
    # vertical band, the *smallest* qualifying rect (plot area lies inside
    # the chart border).
    plot_rects.sort(key=lambda r: (r[1], (r[2] - r[0]) * (r[3] - r[1])))
    bands = []
    for r in plot_rects:
        for band in bands:
            if abs(r[1] - band[0][1]) < 60:
                band.append(r)
                break
        else:
            bands.append([r])
    areas = [min(b, key=lambda r: (r[2] - r[0]) * (r[3] - r[1])) for b in bands]
    areas.sort(key=lambda r: r[1])  # top to bottom: fig10, fig11a, fig11b
    if len(areas) != 3:
        raise RuntimeError(f"expected 3 plot areas, found {len(areas)}: {areas}")
    figures = dict(zip(["fig10", "fig11a", "fig11b"], areas))

    def in_rect(x, y, rect, pad=0.0):
        return (rect[0] - pad <= x <= rect[2] + pad
                and rect[1] - pad <= y <= rect[3] + pad)

    out = {name: [] for name in figures}
    for kind, color, pts, size in elements:
        for x, y in pts:
            for name, rect in figures.items():
                if not in_rect(x, y, rect, pad=0.5):
                    continue
                if any(in_rect(x, y, lb, pad=2.0) for lb in legend_rects):
                    continue
                x0, x1, xlog, y0, y1, ylog = AXES[name]
                fx = (x - rect[0]) / (rect[2] - rect[0])
                fy = (rect[3] - y) / (rect[3] - rect[1])  # page y is down
                if xlog:
                    xd = 10 ** (np.log10(x0) + fx * (np.log10(x1) - np.log10(x0)))
                else:
                    xd = x0 + fx * (x1 - x0)
                if ylog:
                    yd = 10 ** (np.log10(y0) + fy * (np.log10(y1) - np.log10(y0)))
                else:
                    yd = y0 + fy * (y1 - y0)
                out[name].append((xd, yd, color, kind, size))

    refdir = HERE / "reference"
    refdir.mkdir(exist_ok=True)
    for name, rows in out.items():
        rows.sort()
        with (refdir / f"{name}_points.csv").open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["x", "y", "color", "source", "size_pt"])
            for xd, yd, color, kind, size in rows:
                w.writerow([f"{xd:.6g}", f"{yd:.6g}", color, kind, f"{size:.3g}"])
        print(f"{name}: {len(rows)} points -> reference/{name}_points.csv")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    extract(Path(sys.argv[1]))
