# -*- coding: utf-8 -*-
"""Shared helpers: match extracted figure points to computed curves.

The extraction (``extract_reference.py``) yields anonymous data points
(x, y, color); curve identity is resolved here by nearest-curve assignment
in log-log space. Because neighbouring curves in the paper's figures are
separated by factors of ~3-10 while the expected reproduction error is a
few percent, the assignment is unambiguous away from curve crossings; the
reported per-curve statistics are unaffected by how crossing-point ties
break (both candidate curves pass within the tolerance there).
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def load_points(name, color=None, source=None):
    """Load extracted points for figure `name` ('fig10', 'fig11a', 'fig11b')."""
    rows = list(csv.DictReader((HERE / "reference" / f"{name}_points.csv").open()))
    pts = [(float(r["x"]), float(r["y"]), r["color"], r["source"]) for r in rows]
    if color is not None:
        pts = [p for p in pts if p[2] == color]
    if source is not None:
        pts = [p for p in pts if p[3] == source]
    return np.array([(p[0], p[1]) for p in pts])


def match_stats(points, curves, x_range=None, y_floor=None):
    """Assign each extracted point to the nearest computed curve (log-log).

    :param points: (N, 2) array of extracted (x, y)
    :param curves: dict name -> (x_grid, y_values); log-log interpolated
    :param x_range: optional (xmin, xmax) to restrict the comparison
    :param y_floor: drop points below this y (e.g. at the plot frame floor)
    :returns: dict name -> {n, median_ratio, max_abs_dlog10, ...}
    """
    interp = {}
    for name, (xg, yg) in curves.items():
        pos = np.asarray(yg) > 0
        interp[name] = (np.log10(np.asarray(xg)[pos]), np.log10(np.asarray(yg)[pos]))

    assigned = {name: [] for name in curves}
    for x, y in points:
        if x_range and not (x_range[0] <= x <= x_range[1]):
            continue
        if y_floor is not None and y <= y_floor:
            continue
        lx, ly = np.log10(x), np.log10(y)
        best, best_d = None, None
        for name, (cx, cy) in interp.items():
            if lx < cx[0] or lx > cx[-1]:
                continue
            d = ly - np.interp(lx, cx, cy)
            if best is None or abs(d) < abs(best_d):
                best, best_d = name, d
        if best is not None:
            assigned[best].append(best_d)

    stats = {}
    for name, dlogs in assigned.items():
        if not dlogs:
            stats[name] = {"n": 0}
            continue
        dlogs = np.array(dlogs)
        ratios = 10 ** dlogs
        stats[name] = {
            "n": int(dlogs.size),
            "median_ratio": float(np.median(ratios)),
            "median_abs_relerr": float(np.median(np.abs(ratios - 1))),
            "p90_abs_relerr": float(np.percentile(np.abs(ratios - 1), 90)),
            "max_abs_relerr": float(np.max(np.abs(ratios - 1))),
        }
    return stats
