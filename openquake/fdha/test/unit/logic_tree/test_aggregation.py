import numpy as np

from openquake.fdha.logic_tree.aggregation import weighted_mean, weighted_fractiles


def test_weighted_mean_two_branch():
    rates = np.array(
        [
            [[1.0, 2.0]],
            [[3.0, 4.0]],
        ]
    )  # (2,1,2)
    w = np.array([0.3, 0.7])
    mean = weighted_mean(rates, w)
    assert mean.shape == (1, 2)
    assert np.allclose(mean, [[0.3 * 1.0 + 0.7 * 3.0, 0.3 * 2.0 + 0.7 * 4.0]])


def test_weighted_fractiles_monotone():
    rates = np.array(
        [
            [[1.0]],
            [[2.0]],
            [[3.0]],
        ]
    )
    w = np.array([0.2, 0.3, 0.5])
    fr = weighted_fractiles(rates, w, qs=(0.5,))
    assert np.allclose(fr[0.5], [[2.0]])

