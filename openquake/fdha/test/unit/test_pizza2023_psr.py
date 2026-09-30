# -*- coding: utf-8 -*-
"""
Pizza2023PrimarySR against its published regression.

Pizza et al. (2023), BSSA 113(5), Eq. 1:
P = exp(a + b M) / (1 + exp(a + b M)), with the coefficients of Table 2
("Data Fitting and Goodness-of-Fit Parameters for the Equations Derived in
the Current Study") for the four data subsets, each checked over its own
Mw range from the same table.
"""
import numpy as np
import pytest

from openquake.fdha.primary_surf_rup import Pizza2023PrimarySR

pytestmark = pytest.mark.unit

# Table 2: subset -> (a, b, Mw min, Mw max)
TABLE_2 = {
    "all": (-14.47, 2.177, 5.5, 7.9),
    "normal": (-13.50, 2.159, 5.5, 7.1),
    "reverse": (-10.75, 1.427, 5.5, 7.9),
    "strike-slip": (-28.56, 4.436, 5.5, 7.8),
}


def _eq1(a, b, mag):
    fx = a + b * mag
    return np.exp(fx) / (1.0 + np.exp(fx))


@pytest.mark.parametrize("style", sorted(TABLE_2))
def test_matches_table_2(style):
    a, b, m0, m1 = TABLE_2[style]
    mags = np.linspace(m0, m1, 25)
    np.testing.assert_allclose(
        Pizza2023PrimarySR().get_prob(mags, style=style), _eq1(a, b, mags),
        rtol=1e-12)


@pytest.mark.parametrize("style", sorted(TABLE_2))
def test_pinned_style_equals_call_time_style(style):
    mags = np.array([5.5, 6.5, 7.5])
    np.testing.assert_array_equal(
        Pizza2023PrimarySR(style=style).get_prob(mags),
        Pizza2023PrimarySR().get_prob(mags, style=style))


def test_default_subset_is_all_types():
    mags = np.array([5.5, 6.5, 7.5])
    a, b, _, _ = TABLE_2["all"]
    np.testing.assert_allclose(
        Pizza2023PrimarySR().get_prob(mags), _eq1(a, b, mags), rtol=1e-12)


def test_unknown_style_raises():
    with pytest.raises(ValueError):
        Pizza2023PrimarySR().get_prob(6.5, style="oblique")
