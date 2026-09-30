# -*- coding: utf-8 -*-
"""
Takao2014SecondarySR against its published regression.

Takao et al. (2014), J. JAEE 14(2), 16-36 (in Japanese), Eq. (1), p. 19:
P2d = exp(z) / (1 + exp(z)), z = C1 + C2 ln(r + C3), with r the distance
from the principal fault in km and, per unit cell size,

    500 m x 500 m:  C1 = -3.859, C2 = -1.499, C3 = 0.2
    250 m x 250 m:  C1 = -4.903, C2 = -1.459, C3 = 0.2
    100 m x 100 m:  C1 = -6.135, C2 = -1.427, C3 = 0.2
     50 m x  50 m:  C1 = -6.988, C2 = -1.410, C3 = 0.2

The text below the equation reads the fitted curves (their Fig. 2) at 5 km
from the principal fault: 1.8e-3, 6.7e-4, 2.1e-4 and 9.0e-5 for the four
cell sizes. The model is independent of magnitude.
"""
import numpy as np
import pytest

from openquake.fdha.secondary_surf_rup import Takao2014SecondarySR

pytestmark = pytest.mark.unit

EQ_1 = {
    500: (-3.859, -1.499, 0.2),
    250: (-4.903, -1.459, 0.2),
    100: (-6.135, -1.427, 0.2),
    50: (-6.988, -1.410, 0.2),
}
AT_5_KM = {500: 1.8e-3, 250: 6.7e-4, 100: 2.1e-4, 50: 9.0e-5}
R_KM = np.array([0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0])


@pytest.mark.parametrize("cell", sorted(EQ_1))
def test_matches_eq_1(cell):
    c1, c2, c3 = EQ_1[cell]
    z = c1 + c2 * np.log(R_KM + c3)
    np.testing.assert_allclose(
        Takao2014SecondarySR().get_prob(R_KM, pixel_size=cell),
        np.exp(z) / (1.0 + np.exp(z)), rtol=1e-12)


@pytest.mark.parametrize("cell", sorted(AT_5_KM))
def test_matches_values_quoted_at_5_km(cell):
    """Independent of the transcribed coefficients: fixes their signs and
    the km unit (quoted to two significant figures)."""
    assert float(Takao2014SecondarySR().get_prob(5.0, pixel_size=cell)) == \
        pytest.approx(AT_5_KM[cell], rel=0.03)


def test_pinned_cell_size_equals_call_time_cell_size():
    np.testing.assert_array_equal(
        Takao2014SecondarySR(pixel_size=250).get_prob(R_KM),
        Takao2014SecondarySR().get_prob(R_KM, pixel_size=250))


def test_uncalibrated_cell_size_raises():
    with pytest.raises(ValueError):
        Takao2014SecondarySR().get_prob(1.0, pixel_size=200)
