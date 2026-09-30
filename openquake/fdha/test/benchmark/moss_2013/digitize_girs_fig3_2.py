#!/usr/bin/env python
"""Digitise the four Moss et al. (2013) P(sr | Mw) curves of GIRS-2022-05 Fig. 3.2.

Moss et al. (2013) is not available as a released table or code; the
curves are reproduced in Figure 3.2 of the GIRS-2022-05 report (Moss et
al., 2022, revised 1/17/2024, p. 7 of the report = PDF page 20), captioned
"after Moss et al., 2013". The figure is one embedded raster image
(1050 x 1509 px, 222 ppi) with two panels:

    top     (a) Reverse faulting       x: Mw 5..9    y: 0..100 %
    bottom  (a) Strike Slip faulting   x: Mw 4..9    y: 0..100 %

each with a purple solid "Stiff Soil/Rock" curve and a green dashed "Soft
Soil" curve. This script extracts the image with ``pdfimages``, calibrates
both axes of each panel by least squares on the gray major gridlines (every
10 % and every magnitude unit), masks the two curve colours, and writes the
column-wise centre line of each curve to ``reference/girs_fig3_2.csv``
(columns: panel, curve, mw, p). The legend swatches are masked out.

Usage::

    python digitize_girs_fig3_2.py /path/to/Moss_REV_1.17.2024.pdf
"""
import csv
import os
import subprocess
import sys
import tempfile

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "reference", "girs_fig3_2.csv")
PDF_PAGE = 20

# ink colours of the two curves (RGB), from the image palette
COLOURS = {"stiff": (122, 37, 166), "soft": (110, 149, 41)}
COLOUR_TOL = 60.0          # max RGB distance counted as curve ink
GRID_GRAY = (100, 200)     # gray level of gridlines and frame

# per panel: pixel row band, magnitude range, gridline magnitudes, and the
# legend rectangle (x0, x1, y0, y1) to ignore
PANELS = {
    "reverse": dict(rows=(10, 630), mw_grid=[5, 6, 7, 8, 9],
                    legend=(630, 730, 515, 605)),
    "strike-slip": dict(rows=(783, 1400), mw_grid=[4, 5, 6, 7, 8, 9],
                        legend=(630, 730, 1287, 1377)),
}
STEP_PX = 3                # keep one digitised point every STEP_PX columns


def extract_image(pdf):
    tmp = tempfile.mkdtemp()
    subprocess.run(["pdfimages", "-f", str(PDF_PAGE), "-l", str(PDF_PAGE),
                    "-png", pdf, os.path.join(tmp, "fig")], check=True)
    return np.asarray(Image.open(os.path.join(tmp, "fig-000.png"))
                      .convert("RGB")).astype(float)


def _line_centres(profile, threshold):
    """Centres of runs of ``profile > threshold`` (pixel indices)."""
    idx = np.flatnonzero(profile > threshold)
    runs = np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1)
    return [r.mean() for r in runs if r.size]


def calibrate(img, rows, mw_grid):
    """Linear pixel -> (Mw, P) maps from the panel's gridlines and frame."""
    r0, r1 = rows
    sub = img[r0:r1]
    lvl = sub.mean(axis=2)
    neutral = np.ptp(sub, axis=2) < 12
    ink = neutral & (lvl < GRID_GRAY[1])  # gridlines, frame and axis text
    # horizontal gridlines: rows where much of the plot width is gray ink
    # (the 10 % line is partly hidden by the legend box)
    ys = [r0 + y for y in _line_centres(ink[:, 180:960].sum(axis=1), 350)]
    # vertical gridlines + frame: columns inked over most of the panel height
    xs = _line_centres(ink.sum(axis=0), 0.6 * (r1 - r0))
    xs = [x for x in xs if x > 150]  # drop the y-axis label column
    if len(ys) != 11 or len(xs) != len(mw_grid):
        raise RuntimeError(f"gridline detection failed: {len(ys)} rows "
                           f"{len(xs)} cols")
    p_of_y = np.polyfit(ys, np.linspace(1.0, 0.0, 11), 1)
    mw_of_x = np.polyfit(xs, mw_grid, 1)
    resid_y = np.polyval(p_of_y, ys) - np.linspace(1.0, 0.0, 11)
    resid_x = np.polyval(mw_of_x, xs) - np.asarray(mw_grid, float)
    return mw_of_x, p_of_y, xs, float(np.abs(resid_y).max()), \
        float(np.abs(resid_x).max())


def trace(img, colour, rows, x_range, legend):
    r0, r1 = rows
    dist = np.sqrt(((img - np.asarray(colour)) ** 2).sum(axis=2))
    mask = dist < COLOUR_TOL
    mask[:r0] = False
    mask[r1:] = False
    lx0, lx1, ly0, ly1 = legend
    mask[ly0:ly1, lx0:lx1] = False
    pts = []
    for x in range(int(x_range[0]) + 2, int(x_range[1]) - 1):
        ys = np.flatnonzero(mask[:, x])
        if ys.size >= 3 and np.ptp(ys) < 25:  # one clean stroke only
            pts.append((x, ys.mean()))
    return pts[::STEP_PX]


def main(pdf):
    img = extract_image(pdf)
    rows_out = []
    for panel, spec in PANELS.items():
        mw_of_x, p_of_y, xs, ry, rx = calibrate(img, spec["rows"],
                                                spec["mw_grid"])
        print(f"{panel}: calibration residual {ry * 100:.2f} % (y), "
              f"{rx:.3f} Mw (x)")
        for curve, colour in COLOURS.items():
            pts = trace(img, colour, spec["rows"], (xs[0], xs[-1]),
                        spec["legend"])
            for x, y in pts:
                rows_out.append((panel, curve,
                                 round(float(np.polyval(mw_of_x, x)), 4),
                                 round(float(np.polyval(p_of_y, y)), 4)))
            print(f"  {curve}: {len(pts)} points")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["panel", "curve", "mw", "p"])
        w.writerows(rows_out)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main(sys.argv[1])
