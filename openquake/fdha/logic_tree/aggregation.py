"""Weighted aggregation of per-branch rate curves: mean and fractiles.

Relation to openquake.hazardlib.stats
-------------------------------------
``weighted_fractiles`` DELEGATES to :func:`openquake.hazardlib.stats.quantile_curve`.
The two are provably numerically equivalent: ``quantile_curve`` sorts each
element's branch values, forms the cumulative-weight CDF and evaluates
``numpy.interp(q, cdf, values)``. ``numpy.interp`` clamps ``q`` outside
``[cdf[0], cdf[-1]]`` to the endpoint values, exactly reproducing this
module's historical left-clamp/right-clamp convention, and interpolates
linearly between adjacent CDF points otherwise (identical to the old
argsort/cumsum/searchsorted loop). A probe over equal weights, weights
summing to 1±ulps (3x1/3, 7x1/7), a zero-weight branch, tied branch values,
q below the first / above the last cumulative weight, q exactly on a CDF
node, and single/two/random-cube inputs shows a maximum relative difference
of 1.7e-16 (ulp level; most cases bit-identical). The fractile values this
function returns are therefore what the engine writes for its ``quantiles``
outputs: classical PSHA binds ``functools.partial(stats.quantile_curve, q)``
per ``oqparam.quantiles`` entry (``commonlib.oqvalidation.OqParam
.hazard_stats``) and applies it to realization weights renormalized to sum
to 1 (``hazardlib.logictree.FullLogicTree.get_realizations``).
``quantile_curve`` does not renormalize its weights, so the renormalization
the callers rely on is still performed here before delegating.

``weighted_mean`` DELEGATES to
:func:`openquake.hazardlib.stats.mean_curve` (``numpy.average(values,
axis=0, weights=weights)``), so the aggregate ``mean`` curve is identical to
the engine's ``mean`` output. ``mean_curve`` renormalizes by ``sum(weights)``,
which is harmless here: the driver already normalizes the branch weights to
sum to 1 before calling (``logic_tree/driver.py``), so the double
normalization is a no-op on production inputs. (It also makes the function
robust to an un-normalized caller instead of silently returning a scaled
result.) The former local ``tensordot`` differed from this only at ULP level
for normalized weights; the engine's weighted mean is the coherent choice.
"""
from __future__ import annotations

import numpy as np

from openquake.hazardlib.stats import mean_curve, quantile_curve


def weighted_mean(rates: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """
    rates: (n_branches, ...), weights: (n_branches,)

    Weighted mean over the branch axis, delegated to
    :func:`openquake.hazardlib.stats.mean_curve` (``numpy.average``) so the
    aggregate ``mean`` curve is byte-for-byte what the OpenQuake engine
    writes for its own ``mean`` output. See the module docstring.
    """
    rates = np.asarray(rates, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if rates.shape[0] != weights.shape[0]:
        raise ValueError("rates and weights mismatch on branch axis")
    return mean_curve(rates, weights)


def weighted_fractiles(
    rates: np.ndarray,
    weights: np.ndarray,
    qs=(0.05, 0.16, 0.5, 0.84, 0.95),
) -> dict[float, np.ndarray]:
    """
    Weighted empirical fractiles along axis 0 (branch axis).

    For each element position in the remaining dimensions, sort branch values,
    compute cumulative weights, and linearly interpolate within the CDF,
    clamping ``q`` outside ``[cdf[0], cdf[-1]]`` to the extreme branch values.

    This is a thin wrapper over
    :func:`openquake.hazardlib.stats.quantile_curve`; the two are numerically
    equivalent (see the module docstring), so the values returned here are
    exactly what the OpenQuake engine writes for its ``quantiles`` outputs.
    Weights are renormalized to sum to 1 before delegating, matching the
    historical convention and independent of ``quantile_curve`` (which does
    not renormalize).
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

    out: dict[float, np.ndarray] = {}
    for q in qs:
        out[float(q)] = quantile_curve(q, rates, w)
    return out
