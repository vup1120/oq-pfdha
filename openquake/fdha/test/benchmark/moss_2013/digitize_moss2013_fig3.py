#!/usr/bin/env python
"""Digitise Figure 3 of Moss et al. (2013), the source publication.

The figure is one embedded greyscale JPEG (515 x 703 px, 150 ppi) on PDF
page 4 of the open-access copy, with two panels:

    (a) Reverse       x: Mw 5..8.5   y: 0..100 %
    (b) Strike slip   x: Mw 4..8.5   y: 0..100 %

each with a dashed mid-grey "Stiff" curve and a solid dark "Soft" curve;
the legend of each panel prints the four logistic regressions. Both axes
are calibrated by least squares on the light-grey major gridlines. In every
column the dark ink is split into runs; a run whose darkest pixel is below
``SOLID_MAX`` is the solid soft-soil curve, a run whose darkest pixel lies
in ``DASHED_RANGE`` is the dashed stiff curve, and runs longer than
``MAX_RUN`` (where the two curves touch) and the columns on the vertical
gridlines are skipped; points farther than
``OUTLIER_PX`` from the rolling median of their own curve are dropped. The
regression text and the legend are masked out. Output: ``reference/moss2013_fig3.csv``
(columns: panel, curve, mw, p).

Usage::

    python digitize_moss2013_fig3.py /path/to/Mossetal.2013-OpenAccess.pdf
"""
import csv
import os
import subprocess
import sys
import tempfile

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "reference", "moss2013_fig3.csv")
PDF_PAGE = 4

INK_MAX = 140              # grey level below which a pixel is curve ink
SOLID_MAX = 75             # darkest pixel of the solid (soft) curve
DASHED_RANGE = (80, 135)   # darkest pixel of the dashed (stiff) curve
MAX_RUN = 12               # longer runs are two touching curves: skipped
OUTLIER_PX = 2.5           # max offset from the curve's own rolling median
STEP_PX = 2

# per panel: row band, gridline values (top to bottom, P) and magnitudes of
# the vertical gridlines, plus rectangles (x0, x1, y0, y1) to mask
PANELS = {
    "reverse": dict(rows=(28, 300), p_grid=np.linspace(1.0, 0.0, 11),
                    mw_grid=[5, 6, 7, 8], x_max=508,
                    masks=[(88, 245, 38, 96), (340, 510, 222, 265)]),
    # panel (b): the 100 % line sits on the frame; first gridline is 90 %
    "strike-slip": dict(rows=(418, 662), p_grid=np.linspace(0.9, 0.0, 10),
                        mw_grid=[4, 5, 6, 7, 8], x_max=508,
                        masks=[(88, 250, 402, 456), (335, 510, 585, 625)]),
}


def extract_image(pdf):
    tmp = tempfile.mkdtemp()
    subprocess.run(["pdfimages", "-f", str(PDF_PAGE), "-l", str(PDF_PAGE),
                    "-j", pdf, os.path.join(tmp, "fig")], check=True)
    return np.asarray(Image.open(os.path.join(tmp, "fig-000.jpg"))
                      .convert("L")).astype(int)


def _runs(idx):
    return [r for r in np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1)
            if r.size]


def calibrate(img, spec):
    r0, r1 = spec["rows"]
    sub = img[r0:r1]
    grid = (sub > 100) & (sub < 200)
    ys = [r0 + r.mean() for r in _runs(np.flatnonzero(grid.sum(1) > 250))]
    xs = [r.mean() for r in _runs(np.flatnonzero(
        grid.sum(0) > 0.5 * (r1 - r0)))]
    ys = ys[:len(spec["p_grid"])]
    xs = [x for x in xs if 70 < x < spec["x_max"]][:len(spec["mw_grid"])]
    if len(ys) != len(spec["p_grid"]) or len(xs) != len(spec["mw_grid"]):
        raise RuntimeError(f"gridline detection failed: {ys} {xs}")
    p_of_y = np.polyfit(ys, spec["p_grid"], 1)
    mw_of_x = np.polyfit(xs, spec["mw_grid"], 1)
    ry = np.abs(np.polyval(p_of_y, ys) - spec["p_grid"]).max()
    rx = np.abs(np.polyval(mw_of_x, xs) - spec["mw_grid"]).max()
    return mw_of_x, p_of_y, xs, ry, rx


def trace(img, spec, grid_xs):
    r0, r1 = spec["rows"]
    ink = img < INK_MAX
    ink[:r0] = False
    ink[r1:] = False
    for mx0, mx1, my0, my1 in spec["masks"]:
        ink[my0:my1, mx0:mx1] = False
    pts = {"stiff": [], "soft": []}
    for x in range(int(grid_xs[0]) + 3, spec["x_max"]):
        if min(abs(x - g) for g in grid_xs) < 2:   # gridline + axis tick
            continue
        for run in _runs(np.flatnonzero(ink[:, x])):
            if run.size < 2 or run.size > MAX_RUN:
                continue
            darkest = img[run, x].min()
            w = 255.0 - img[run, x]            # ink-weighted centre
            y = float((run * w).sum() / w.sum())
            if darkest < SOLID_MAX:
                pts["soft"].append((x, y))
            elif DASHED_RANGE[0] <= darkest <= DASHED_RANGE[1]:
                pts["stiff"].append((x, y))
    return {k: _drop_outliers(v)[::STEP_PX] for k, v in pts.items()}


def _drop_outliers(pts, half=6):
    """Remove isolated misclassified runs: keep points within OUTLIER_PX of
    the rolling median of the same curve's neighbouring points (no model
    information is used)."""
    if not pts:
        return pts
    arr = np.array(sorted(pts))
    keep = []
    for i, (x, y) in enumerate(arr):
        near = arr[max(0, i - half):i + half + 1, 1]
        if abs(y - np.median(near)) <= OUTLIER_PX:
            keep.append((x, y))
    return keep


def main(pdf):
    img = extract_image(pdf)
    out = []
    for panel, spec in PANELS.items():
        mw_of_x, p_of_y, grid_xs, ry, rx = calibrate(img, spec)
        print(f"{panel}: calibration residual {ry * 100:.2f} % (y), "
              f"{rx:.3f} Mw (x)")
        for curve, pts in trace(img, spec, grid_xs).items():
            for x, y in pts:
                out.append((panel, curve,
                            round(float(np.polyval(mw_of_x, x)), 4),
                            round(float(np.polyval(p_of_y, y)), 4)))
            print(f"  {curve}: {len(pts)} points")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["panel", "curve", "mw", "p"])
        w.writerows(out)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main(sys.argv[1])
