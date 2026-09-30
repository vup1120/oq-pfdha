# -*- coding: utf-8 -*-
"""
Rodriguez2023SecondarySR against its published regression.

Rodriguez Padilla and Oskin (2023), BSSA 113(6), Eq. 2:
nu(x) = nu_o ((x + x_fr) / x_fr) ** -gamma, the probability of a rupture
per 1 m^2, with x and x_fr in metres, and the "General Model" column of
Table 1: nu_o = 0.13, x_fr = 6.7 m, gamma = 1.19. The model takes the
distance in km.
"""
import numpy as np
import pytest

from openquake.fdha.secondary_surf_rup import Rodriguez2023SecondarySR

pytestmark = pytest.mark.unit

NU_O, X_FR, GAMMA = 0.13, 6.7, 1.19   # Table 1, General Model
X_M = np.array([0.0, 1.0, 6.7, 10.0, 50.0, 100.0, 500.0, 1000.0, 3000.0])


def test_matches_eq_2_general_model():
    got = Rodriguez2023SecondarySR().get_prob(X_M / 1000.0)
    expected = NU_O * ((X_M + X_FR) / X_FR) ** -GAMMA
    np.testing.assert_allclose(got, expected, rtol=1e-12)


def test_anchor_values():
    """nu(0) = nu_o and nu(x_fr) = nu_o 2^-gamma, independent of units."""
    m = Rodriguez2023SecondarySR()
    assert m.get_prob(0.0) == pytest.approx(NU_O, rel=1e-12)
    assert m.get_prob(X_FR / 1000.0) == pytest.approx(
        NU_O * 2.0 ** -GAMMA, rel=1e-12)


def test_only_the_1_m_pixel_is_calibrated():
    with pytest.raises(ValueError):
        Rodriguez2023SecondarySR().get_prob(0.1, pixel_size=100)
