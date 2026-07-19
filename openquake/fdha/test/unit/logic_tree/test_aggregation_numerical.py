"""Numerical correctness tests for the logic-tree aggregation statistics.

Modelled on the OpenQuake engine's logic-tree tests
(``openquake/commonlib/tests/logictree_test.py`` and
``openquake/hazardlib/tests/lt_test.py``): every expected value is either
worked out by hand in the test body or produced by an *independent*
brute-force implementation of the weighted-percentile definition, never by
the code under test itself.

The weighted-fractile convention verified here is the one documented in
:mod:`openquake.fdha.logic_tree.aggregation`: sort the branch values, build
the cumulative-weight CDF, clamp ``q`` outside ``[cdf[0], cdf[-1]]`` to the
extreme values, and interpolate linearly between adjacent CDF points
otherwise.
"""
from __future__ import annotations

import bisect

import numpy as np
import pytest

from openquake.fdha.logic_tree.aggregation import weighted_fractiles, weighted_mean


# ------------------------------------------------------------------------
# Independent reference implementation of the weighted percentile.
# Written from the definition (piecewise-linear inverse of the weighted
# empirical CDF), using bisect + plain Python so it shares no code with
# aggregation.weighted_fractiles.
# ------------------------------------------------------------------------

def _reference_weighted_percentile(values, weights, q):
    pairs = sorted(zip(values, weights))
    xs = [p[0] for p in pairs]
    ws = [p[1] for p in pairs]
    total = sum(ws)
    cdf = []
    acc = 0.0
    for w in ws:
        acc += w / total
        cdf.append(acc)
    if q <= cdf[0]:
        return xs[0]
    if q >= cdf[-1]:
        return xs[-1]
    j = bisect.bisect_left(cdf, q)
    c0, c1 = cdf[j - 1], cdf[j]
    x0, x1 = xs[j - 1], xs[j]
    if c1 == c0:
        return x1
    return x0 + (q - c0) / (c1 - c0) * (x1 - x0)


def _reference_fractile_cube(rates, weights, q):
    rates = np.asarray(rates, dtype=float)
    flat = rates.reshape(rates.shape[0], -1)
    out = np.array(
        [
            _reference_weighted_percentile(flat[:, i].tolist(), list(weights), q)
            for i in range(flat.shape[1])
        ]
    )
    return out.reshape(rates.shape[1:])


# ----------------------------------------------------------- weighted_mean


def test_weighted_mean_hand_computed_cube():
    # (3 branches, 2 sites, 2 D0 levels) with non-uniform weights.
    rates = np.array(
        [
            [[1.0, 2.0], [3.0, 4.0]],
            [[5.0, 6.0], [7.0, 8.0]],
            [[9.0, 10.0], [11.0, 12.0]],
        ]
    )
    w = np.array([0.5, 0.3, 0.2])
    mean = weighted_mean(rates, w)
    expected = np.array(
        [
            [0.5 * 1 + 0.3 * 5 + 0.2 * 9, 0.5 * 2 + 0.3 * 6 + 0.2 * 10],
            [0.5 * 3 + 0.3 * 7 + 0.2 * 11, 0.5 * 4 + 0.3 * 8 + 0.2 * 12],
        ]
    )  # = [[3.8, 4.8], [5.8, 6.8]]
    assert mean.shape == (2, 2)
    # rtol=1e-12 (not exact): weighted_mean delegates to hazardlib
    # stats.mean_curve (numpy.average), whose multiply+sum reduction differs
    # from a hand-computed sum only at ULP level. The weighted mean is a
    # physical quantity, not a specific summation order.
    np.testing.assert_allclose(mean, expected, rtol=1e-12)


def test_weighted_mean_matches_engine_mean_curve():
    # Coherence contract: weighted_mean IS the engine's weighted mean, so the
    # aggregate `mean` curve equals what OpenQuake writes for its own `mean`
    # output. In production the driver normalizes weights to sum to 1 before
    # calling; mean_curve renormalizes too, so on normalized weights the
    # result is the plain weighted average.
    from openquake.hazardlib.stats import mean_curve

    rates = np.array([[[2.0]], [[4.0]]])
    w = np.array([0.5, 0.5])  # already normalized, as the driver supplies
    got = weighted_mean(rates, w)
    np.testing.assert_allclose(got, [[3.0]])
    np.testing.assert_allclose(got, mean_curve(rates, w), rtol=0, atol=0)


def test_weighted_mean_branch_axis_mismatch_raises():
    with pytest.raises(ValueError, match="mismatch"):
        weighted_mean(np.zeros((3, 2)), np.array([0.5, 0.5]))


# ------------------------------------------------------ weighted_fractiles


def test_weighted_fractiles_hand_computed_interpolation():
    # Branch values 10, 2, 6 with weights .5, .2, .3.
    # Sorted: xs = [2, 6, 10], ws = [.2, .3, .5], cdf = [.2, .5, 1.0].
    #   q=0.05 -> q <= cdf[0]           -> 2
    #   q=0.20 -> q <= cdf[0] (boundary)-> 2
    #   q=0.35 -> between (.2,2) (.5,6): t=(.35-.2)/.3=0.5 -> 4
    #   q=0.50 -> between (.2,2) (.5,6): t=1               -> 6
    #   q=0.75 -> between (.5,6) (1,10): t=0.5             -> 8
    #   q=1.00 -> q >= cdf[-1]          -> 10
    rates = np.array([[[10.0]], [[2.0]], [[6.0]]])
    w = np.array([0.5, 0.2, 0.3])
    qs = (0.05, 0.20, 0.35, 0.50, 0.75, 1.00)
    fr = weighted_fractiles(rates, w, qs=qs)
    expected = {0.05: 2.0, 0.20: 2.0, 0.35: 4.0, 0.50: 6.0, 0.75: 8.0, 1.00: 10.0}
    for q, exp in expected.items():
        # rtol accounts for one ulp introduced by the cumulative-weight sum.
        np.testing.assert_allclose(fr[q], [[exp]], rtol=1e-12)


def test_weighted_fractiles_tied_values():
    # xs = [1, 5, 5], ws = [.5, .25, .25], cdf = [.5, .75, 1.0].
    #   q=0.60 -> between (.5,1) (.75,5): t=0.4 -> 1 + 0.4*4 = 2.6
    #   q=0.80 -> between (.75,5) (1,5)         -> 5
    rates = np.array([[[5.0]], [[5.0]], [[1.0]]])
    w = np.array([0.25, 0.25, 0.5])
    fr = weighted_fractiles(rates, w, qs=(0.60, 0.80))
    np.testing.assert_allclose(fr[0.60], [[2.6]])
    np.testing.assert_allclose(fr[0.80], [[5.0]])


def test_weighted_fractiles_unnormalized_weights_equivalent():
    rates = np.array([[[10.0]], [[2.0]], [[6.0]]])
    qs = (0.05, 0.35, 0.75)
    fr_norm = weighted_fractiles(rates, np.array([0.5, 0.2, 0.3]), qs=qs)
    fr_scaled = weighted_fractiles(rates, np.array([5.0, 2.0, 3.0]), qs=qs)
    for q in qs:
        np.testing.assert_allclose(fr_scaled[q], fr_norm[q], rtol=0, atol=0)


def test_weighted_fractiles_single_branch_returns_that_branch():
    rates = np.array([[[3.0, 7.0], [1.0, 0.0]]])
    fr = weighted_fractiles(rates, np.array([1.0]), qs=(0.05, 0.5, 0.95))
    for q in (0.05, 0.5, 0.95):
        np.testing.assert_allclose(fr[q], rates[0], rtol=0, atol=0)


def test_weighted_fractiles_matches_bruteforce_reference():
    rng = np.random.default_rng(42)
    rates = rng.lognormal(mean=-8.0, sigma=1.5, size=(9, 4, 5))
    weights = rng.uniform(0.05, 1.0, size=9)
    qs = (0.01, 0.05, 0.16, 0.33, 0.5, 0.84, 0.95, 0.99)
    fr = weighted_fractiles(rates, weights, qs=qs)
    for q in qs:
        expected = _reference_fractile_cube(rates, weights, q)
        np.testing.assert_allclose(fr[q], expected, rtol=1e-12)


def test_weighted_fractiles_monotone_in_q_and_bounded():
    rng = np.random.default_rng(7)
    rates = rng.exponential(scale=1e-4, size=(7, 3, 4))
    weights = rng.uniform(0.1, 1.0, size=7)
    qs = np.linspace(0.01, 0.99, 25)
    fr = weighted_fractiles(rates, weights, qs=qs)
    stacked = np.stack([fr[float(q)] for q in qs], axis=0)
    assert np.all(np.diff(stacked, axis=0) >= 0), "fractiles must not decrease in q"
    assert np.all(stacked >= rates.min(axis=0))
    assert np.all(stacked <= rates.max(axis=0))


def test_weighted_fractiles_branch_permutation_invariant():
    rng = np.random.default_rng(11)
    rates = rng.random(size=(6, 2, 3))
    weights = rng.uniform(0.1, 1.0, size=6)
    perm = rng.permutation(6)
    qs = (0.05, 0.5, 0.95)
    fr = weighted_fractiles(rates, weights, qs=qs)
    fr_perm = weighted_fractiles(rates[perm], weights[perm], qs=qs)
    for q in qs:
        # The cumulative sum runs in a different order after the permutation,
        # so allow float-epsilon differences.
        np.testing.assert_allclose(fr_perm[q], fr[q], rtol=1e-12)


def test_weighted_fractiles_shape_preserved():
    rates = np.zeros((4, 5, 6))
    fr = weighted_fractiles(rates, np.ones(4), qs=(0.5,))
    assert fr[0.5].shape == (5, 6)


def test_weighted_fractiles_invalid_weights_raise():
    rates = np.zeros((2, 1, 1))
    with pytest.raises(ValueError, match="invalid weights"):
        weighted_fractiles(rates, np.array([0.5, -0.1]))
    with pytest.raises(ValueError, match="invalid weights"):
        weighted_fractiles(rates, np.array([0.5, np.nan]))
    with pytest.raises(ValueError, match="zero total weight"):
        weighted_fractiles(rates, np.array([0.0, 0.0]))
    with pytest.raises(ValueError, match="mismatch"):
        weighted_fractiles(rates, np.array([1.0]))
