"""Weighted aggregation of per-branch rate curves: mean and fractiles."""
from __future__ import annotations

import numpy as np


def weighted_mean(rates: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """
    rates: (n_branches, ...), weights: (n_branches,)
    """
    rates = np.asarray(rates, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if rates.shape[0] != weights.shape[0]:
        raise ValueError("rates and weights mismatch on branch axis")
    return np.tensordot(weights, rates, axes=(0, 0))


def weighted_fractiles(
    rates: np.ndarray,
    weights: np.ndarray,
    qs=(0.05, 0.16, 0.5, 0.84, 0.95),
) -> dict[float, np.ndarray]:
    """
    Weighted empirical fractiles along axis 0 (branch axis).

    For each element position in the remaining dimensions, sort branch values,
    compute cumulative weights, and linearly interpolate within the CDF.
    """
    rates = np.asarray(rates, dtype=float)
    w = np.asarray(weights, dtype=float)
    if rates.shape[0] != w.shape[0]:
        raise ValueError("rates and weights mismatch on branch axis")
    if np.any(w < 0) or not np.isfinite(w).all():
        raise ValueError("invalid weights")
    if w.sum() == 0:
        raise ValueError("zero total weight")
    w = w / w.sum()

    flat = rates.reshape((rates.shape[0], -1))
    out: dict[float, np.ndarray] = {}
    for q in qs:
        vals = np.empty(flat.shape[1], dtype=float)
        for i in range(flat.shape[1]):
            x = flat[:, i]
            idx = np.argsort(x)
            xs = x[idx]
            ws = w[idx]
            cdf = np.cumsum(ws)
            # leftmost
            if q <= cdf[0]:
                vals[i] = xs[0]
                continue
            # rightmost
            if q >= cdf[-1]:
                vals[i] = xs[-1]
                continue
            j = int(np.searchsorted(cdf, q, side="left"))
            x0, x1 = xs[j - 1], xs[j]
            c0, c1 = cdf[j - 1], cdf[j]
            if c1 == c0:
                vals[i] = x1
            else:
                t = (q - c0) / (c1 - c0)
                vals[i] = x0 + t * (x1 - x0)
        out[float(q)] = vals.reshape(rates.shape[1:])
    return out

