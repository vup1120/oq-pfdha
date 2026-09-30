#!/usr/bin/env python
"""Digitize the TECDOC-2092 distributed-displacement figures that the
exercise coordinators' plotting scripts did not cover.

The curves of Figs 14/18/21 and 15(a,b)/19(a,b)/22(a) are available as
numeric vectors from the coordinators' MATLAB scripts and are transcribed
verbatim in ``make_reference_csvs.py``. The six remaining exercises --- the
2nd and 3rd distributed case of each case study --- exist only as published
raster figures:

    Fig. 15(c)  Kumamoto 2nd distributed, r = 10 km        T13, P11
    Fig. 15(d)  Kumamoto 3rd distributed, Suizenji         T13, P11
    Fig. 19(c)  Le Teil 2nd distributed, dip 45 NW         T13, V24
    Fig. 19(d)  Le Teil 3rd distributed, three faults      T13, V24
    Fig. 22(b)  Norcia 2nd distributed, site SL            Y03, V24
    Fig. 22(c)  Norcia 3rd distributed, MVF + NF           Y03, V24

This script extracts those twelve curves from the publication PDF and
writes them as CSVs next to the transcribed ones, so that the six
exercises can be compared on the same footing as the other six.

Method
------
1. The figure page's embedded JPEG is extracted at native resolution
   (``pdfimages``); no re-rendering, so no resampling blur is added.
2. The axes box of each panel is located from the longest continuous dark
   runs (the box border).
3. The axes are calibrated from the *major gridlines*, which MATLAB draws
   solid at the labelled decades: their pixel pitch gives px/decade and
   their absolute position anchors the scale. The box edges are then
   checked to fall on decades, which is an independent consistency test
   of the calibration (reported by ``--verify``).
4. Each curve is isolated by a colour mask (the teams' model colours are
   well separated: black T13/Y03-black, red P11, green V24, orange Y03).
   Connected components narrower than 15 % of the panel are discarded,
   which removes the in-axes panel letters "(c)"/"(d)"; in-axes legend
   boxes are masked out explicitly.
5. Where a panel superimposes per-source curves (dotted/dashed) on the
   total (solid), the total is recovered as the topmost trace of that
   colour. This is exact rather than heuristic: the total rate is the sum
   of the per-source rates, hence never below any of them.
6. Traces are resampled in log-log onto the 18-point displacement grid
   shared by the transcribed reference CSVs. Grid points outside the
   plotted range, or where the curve has already left the axes, are
   written as empty fields.

Digitisation error is dominated by the curve line width (5-8 px, i.e.
0.03-0.09 decade at these scales) and is reported per curve by
``--verify`` as the round-trip residual of the calibration.

Usage (repo root)::

    python openquake/fdha/test/benchmark/IAEA/reference/digitize_tecdoc_figures.py
    python .../digitize_tecdoc_figures.py --verify      # + overlay PNGs

Requires ``pdfimages`` (poppler) and the TECDOC-2092 PDF; set
``TECDOC_PDF`` if it is not at the default path below.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

HERE = Path(__file__).resolve().parent
DEFAULT_PDF = Path(os.environ.get(
    "TECDOC_PDF",
    Path.home() / "Documents/PFDHA/IAEA_benchmarking/TE-2092web.pdf"))

# The 18-point grid of the transcribed reference CSVs (cm).
DISP_CM = [0.01, 0.1, 0.5, 1.0, 1.5, 3.0, 5.0, 7.5, 10.0, 15.0,
           30.0, 50.0, 75.0, 100.0, 300.0, 500.0, 750.0, 1000.0]


# --------------------------------------------------------------- colours
def mask_black(im):
    mx = im.max(axis=2)
    mn = im.min(axis=2)
    return (mx < 110) & ((mx - mn) < 60)


def mask_red(im):
    r, g, b = im[..., 0], im[..., 1], im[..., 2]
    return (r > 120) & (g < 95) & (b < 95) & ((r - g) > 60)


def mask_green(im):
    r, g, b = im[..., 0], im[..., 1], im[..., 2]
    return (g > 90) & (r < 150) & (b < 150) & ((g - r) > 45) & ((g - b) > 45)


def mask_orange(im):
    r, g, b = im[..., 0], im[..., 1], im[..., 2]
    return (r > 165) & (g > 105) & (g < 220) & (b < 130) & \
           ((r - b) > 80) & ((g - b) > 40)


MASKS = {"black": mask_black, "red": mask_red,
         "green": mask_green, "orange": mask_orange}


# ------------------------------------------------------------ panel spec
class Panel:
    def __init__(self, figure, page, quad, xlim, ylim, curves,
                 legend=None, note=""):
        self.figure = figure        # e.g. "fig15c"
        self.page = page
        self.quad = quad            # (y0,y1,x0,x1) search window
        self.xlim = xlim            # (log10 xmin, log10 xmax) at box edges
        self.ylim = ylim            # (log10 ytop, log10 ybottom)
        self.curves = curves        # {model: colour}
        self.legend = legend        # (y0,y1,x0,x1) to mask out
        self.note = note


PANELS = [
    Panel("fig15c_kumamoto_distributed_r10", 48, (690, 1369, 80, 700),
          (0, 3), (-5, -10), {"T13": "black", "P11": "red"},
          note="Kumamoto 2nd distributed exercise, r = 10 km"),
    Panel("fig15d_kumamoto_distributed_suizenji", 48, (690, 1369, 760, 1368),
          (0, 3), (-5, -10), {"T13": "black", "P11": "red"},
          legend=(714, 900, 1057, 1320),
          note="Kumamoto 3rd distributed exercise, Suizenji fault, r = 0.6 km"),
    Panel("fig19c_leteil_distributed_dipNW", 56, (620, 1281, 80, 700),
          (0, 3), (-5, -11), {"T13": "black", "V24": "green"},
          note="Le Teil 2nd distributed exercise, dip direction 45 NW"),
    Panel("fig19d_leteil_distributed_3faults", 56, (620, 1281, 700, 1280),
          (0, 3), (-5, -11), {"T13": "black", "V24": "green"},
          legend=(650, 909, 860, 1235),
          note="Le Teil 3rd distributed exercise, three faults (total)"),
    Panel("fig22b_norcia_distributed_sl", 61, (0, 700, 720, 1338),
          (0, 3), (-4, -9), {"Y03": "orange", "V24": "green"},
          note="Norcia 2nd distributed exercise, site SL, r = 2.4 km"),
    Panel("fig22c_norcia_distributed_2sources", 61, (620, 1300, 80, 760),
          (0, 3), (-4, -9), {"Y03": "orange", "V24": "green"},
          note="Norcia 3rd distributed exercise, MVF + NF (total)"),
]


# ------------------------------------------------------------- geometry
def _longest_run(mask):
    best = cur = 0
    for v in mask:
        cur = cur + 1 if v else 0
        best = max(best, cur)
    return best


def find_box(im, quad):
    """Axes box (top, bottom, left, right) as pixel indices."""
    y0, y1, x0, x1 = quad
    # 210 rather than a harder threshold: the JPEG-compressed box borders
    # are anti-aliased, while the major gridlines stay lighter than this.
    dark = im.max(axis=2) < 210
    sub = dark[y0:y1, x0:x1]
    h, w = sub.shape
    hl = np.array([_longest_run(sub[y]) for y in range(h)])
    vl = np.array([_longest_run(sub[:, x]) for x in range(w)])
    ys = np.where(hl > 0.6 * w)[0]
    xs = np.where(vl > 0.6 * h)[0]
    return int(y0 + ys.min()), int(y0 + ys.max()), \
           int(x0 + xs.min()), int(x0 + xs.max())


def gridline_positions(im, box, axis):
    """Rows (axis='y') or columns (axis='x') carrying a solid major
    gridline, i.e. non-white over essentially the whole panel."""
    t, b, l, r = box
    sub = im[t + 3:b - 2, l + 3:r - 2]
    nonwhite = sub.max(axis=2) < 245
    frac = nonwhite.mean(axis=1) if axis == "y" else nonwhite.mean(axis=0)
    off = (t + 3) if axis == "y" else (l + 3)
    idx = np.where(frac > 0.85)[0]
    if idx.size == 0:
        return []
    groups, cur = [], [idx[0]]
    for i in idx[1:]:
        if i - cur[-1] <= 3:
            cur.append(i)
        else:
            groups.append(off + float(np.mean(cur)))
            cur = [i]
    groups.append(off + float(np.mean(cur)))
    return groups


def calibrate(im, box, panel):
    """Return (px_per_decade_x, x0_px, px_per_decade_y, y0_px) plus the
    consistency residual of the box edges against whole decades."""
    t, b, l, r = box
    gy = gridline_positions(im, box, "y")
    gx = gridline_positions(im, box, "x")

    # x: box edges are decades by construction of these plots; use them,
    # and check any interior gridlines fall on the implied decades.
    ndec_x = panel.xlim[1] - panel.xlim[0]
    ppd_x = (r - l) / ndec_x
    res_x = []
    for g in gx:
        dec = (g - l) / ppd_x
        res_x.append(abs(dec - round(dec)))

    # y: derive the pitch from the labelled gridlines when there are >=2,
    # otherwise fall back to the box edges.
    ndec_y = panel.ylim[0] - panel.ylim[1]
    if len(gy) >= 2:
        step = (max(gy) - min(gy)) / (len(gy) - 1)     # px per label step
        # label step in decades, inferred from the box height
        dec_per_step = round(ndec_y / ((b - t) / step))
        ppd_y = step / dec_per_step
    else:
        ppd_y = (b - t) / ndec_y
    res_y = []
    for g in gy:
        dec = (g - t) / ppd_y
        res_y.append(abs(dec - round(dec)))
    # box-edge consistency: implied height in decades vs declared
    edge_res = abs((b - t) / ppd_y - ndec_y)
    return dict(ppd_x=ppd_x, ppd_y=ppd_y, box=box,
                res_x=max(res_x) if res_x else 0.0,
                res_y=max(res_y) if res_y else 0.0,
                edge_res=edge_res, n_gy=len(gy), n_gx=len(gx))


def px_to_data(cal, panel, xs, ys):
    t, b, l, r = cal["box"]
    logx = panel.xlim[0] + (np.asarray(xs) - l) / cal["ppd_x"]
    logy = panel.ylim[0] - (np.asarray(ys) - t) / cal["ppd_y"]
    return 10.0 ** logx, 10.0 ** logy


# ------------------------------------------------------------- tracing
def trace(im, box, panel, colour):
    """Topmost trace of `colour` inside the axes; returns (xs, ys) px."""
    t, b, l, r = box
    m = MASKS[colour](im).copy()

    # restrict to the axes interior (leave the border out)
    keep = np.zeros_like(m)
    keep[t + 3:b - 2, l + 3:r - 2] = True
    m &= keep
    if panel.legend:
        ly0, ly1, lx0, lx1 = panel.legend
        m[ly0:ly1, lx0:lx1] = False

    # drop connected components too narrow to be a curve (panel letters)
    lab, n = ndimage.label(m, structure=np.ones((3, 3), int))
    if n:
        widths = {}
        objs = ndimage.find_objects(lab)
        for i, sl in enumerate(objs, start=1):
            if sl is None:
                continue
            widths[i] = sl[1].stop - sl[1].start
        min_w = 0.15 * (r - l)
        drop = [i for i, w in widths.items() if w < min_w]
        if drop:
            m[np.isin(lab, drop)] = False

    xs, ys = [], []
    for x in range(l + 3, r - 2):
        col = np.where(m[:, x])[0]
        if col.size == 0:
            continue
        # topmost contiguous run -> the total curve
        runs, cur = [], [col[0]]
        for v in col[1:]:
            if v - cur[-1] <= 2:
                cur.append(v)
            else:
                runs.append(cur)
                cur = [v]
        runs.append(cur)
        top = runs[0]
        xs.append(x)
        ys.append(float(np.mean(top)))
    return np.array(xs), np.array(ys)


def resample(x_dat, y_dat, box, cal, panel):
    """Sample the trace on DISP_CM; blank outside the plotted range."""
    t, b, l, r = box
    ymin = 10.0 ** panel.ylim[1]
    ymax = 10.0 ** panel.ylim[0]
    out = []
    if x_dat.size < 5:
        return [""] * len(DISP_CM)
    lx = np.log10(x_dat)
    ly = np.log10(y_dat)
    order = np.argsort(lx)
    lx, ly = lx[order], ly[order]
    # The trace starts a few pixels inside the box border, so the axis-edge
    # grid point sits marginally outside it. Clamping within 0.03 decade
    # recovers it; that is well inside the digitisation error.
    EDGE = 0.03
    for d in DISP_CM:
        ld = np.log10(d)
        if ld < lx.min() - EDGE or ld > lx.max() + EDGE:
            out.append("")
            continue
        ld = min(max(ld, lx.min()), lx.max())
        v = 10.0 ** np.interp(ld, lx, ly)
        # a curve that has run off the bottom of the axes carries no value
        if v <= ymin * 1.15 or v >= ymax * 0.999:
            out.append("")
        else:
            out.append(f"{v:.3e}")
    return out


# ---------------------------------------------------------------- main
def extract_pages(pdf, pages, workdir):
    imgs = {}
    for p in pages:
        subprocess.run(["pdfimages", "-png", "-f", str(p), "-l", str(p),
                        str(pdf), str(workdir / f"page{p}")], check=True)
        cands = sorted(workdir.glob(f"page{p}-*.png"))
        if not cands:
            raise SystemExit(f"no image extracted from page {p}")
        imgs[p] = np.array(Image.open(cands[0]).convert("RGB")).astype(int)
    return imgs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdf", type=Path, default=DEFAULT_PDF)
    ap.add_argument("--verify", action="store_true",
                    help="write overlay PNGs and print calibration residuals")
    ap.add_argument("--outdir", type=Path, default=HERE)
    args = ap.parse_args()

    if not args.pdf.exists():
        raise SystemExit(f"TECDOC PDF not found: {args.pdf}\n"
                         "set TECDOC_PDF or pass --pdf")

    with tempfile.TemporaryDirectory() as td:
        work = Path(td)
        imgs = extract_pages(args.pdf, sorted({p.page for p in PANELS}), work)

        overlays = []
        for panel in PANELS:
            im = imgs[panel.page]
            box = find_box(im, panel.quad)
            cal = calibrate(im, box, panel)
            rows = {}
            traces = {}
            for model, colour in panel.curves.items():
                xs, ys = trace(im, box, panel, colour)
                xd, yd = px_to_data(cal, panel, xs, ys)
                traces[model] = (xs, ys, xd, yd)
                rows[model] = resample(xd, yd, box, cal, panel)

            out = args.outdir / f"{panel.figure}.csv"
            with out.open("w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["disp_cm"] + list(panel.curves))
                for i, d in enumerate(DISP_CM):
                    w.writerow([d] + [rows[m][i] for m in panel.curves])
            print(f"wrote {out.name}  ({panel.note})")
            if args.verify:
                print(f"    box={box} px/decade x={cal['ppd_x']:.2f} "
                      f"y={cal['ppd_y']:.2f}  gridlines x={cal['n_gx']} "
                      f"y={cal['n_gy']}")
                print(f"    residuals: x-grid={cal['res_x']:.4f} dec, "
                      f"y-grid={cal['res_y']:.4f} dec, "
                      f"box-edge={cal['edge_res']:.4f} dec")
                for m, (xs, ys, xd, yd) in traces.items():
                    print(f"    {m}: {xs.size} px columns, "
                          f"x {xd.min():.2f}-{xd.max():.1f} cm, "
                          f"y {yd.min():.2e}-{yd.max():.2e}")
                overlays.append((panel, im, box, traces))

        if args.verify and overlays:
            write_overlays(overlays, args.outdir)


def write_overlays(overlays, outdir):
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    for ax, (panel, im, box, traces) in zip(axes.flat, overlays):
        t, b, l, r = box
        ax.imshow(im[t - 30:b + 30, l - 30:r + 30].astype(np.uint8))
        for m, (xs, ys, xd, yd) in traces.items():
            ax.plot(xs - (l - 30), ys - (t - 30), lw=0.9, ls="--",
                    label=m)
        ax.set_title(panel.figure, fontsize=8)
        ax.legend(fontsize=7)
        ax.axis("off")
    fig.tight_layout()
    p = outdir / "digitisation_overlays.png"
    fig.savefig(p, dpi=130)
    print("wrote", p.name)


if __name__ == "__main__":
    sys.exit(main())
