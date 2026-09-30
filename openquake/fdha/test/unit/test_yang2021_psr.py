# -*- coding: utf-8 -*-
"""
Yang2021PrimarySR against its published regression.

Yang et al. (2021), GSA Bulletin 133, Fig. 11A legend: the "SCR Oz"
logistic P(Slip | m) = 1 / (1 + exp(a + b Mw)) with a = 24.59 and
b = -4.00, valid for 4.0 <= Mw <= 6.6 (Fig. 11 caption). Pizza et al.
(2023, BSSA 113, Table 1) restate the same regression in their Eq. 1 form
exp(A + B M) / (1 + exp(A + B M)) as A = -24.59, B = 4.00; both
transcriptions are checked, since a sign slip between the two forms is
the likely implementation error.
"""
import numpy as np
import pytest

from openquake.fdha.primary_surf_rup import Yang2021PrimarySR

pytestmark = pytest.mark.unit

MAGS = np.linspace(4.0, 6.6, 27)


def test_matches_fig_11a_legend():
    a, b = 24.59, -4.00
    expected = 1.0 / (1.0 + np.exp(a + b * MAGS))
    np.testing.assert_allclose(
        Yang2021PrimarySR().get_prob(MAGS), expected, rtol=1e-12)


def test_matches_pizza2023_table_1_restatement():
    a, b = -24.59, 4.00
    fx = a + b * MAGS
    expected = np.exp(fx) / (1.0 + np.exp(fx))
    np.testing.assert_allclose(
        Yang2021PrimarySR().get_prob(MAGS), expected, rtol=1e-12)


def test_half_probability_at_published_midpoint():
    """P = 0.5 where a + b Mw = 0, i.e. Mw = 24.59 / 4 = 6.1475; the
    plotted SCR Oz curve crosses 50 % there."""
    assert Yang2021PrimarySR().get_prob(6.1475) == pytest.approx(
        0.5, abs=1e-12)


def test_scalar_in_scalar_out():
    assert isinstance(Yang2021PrimarySR().get_prob(6.0), float)
