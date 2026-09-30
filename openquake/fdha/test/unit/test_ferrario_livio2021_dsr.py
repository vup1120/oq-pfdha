# -*- coding: utf-8 -*-
"""
FerrarioLivio2021SecondarySR against its published regression.

Ferrario and Livio (2021), Solid Earth 12, Eq. 2:
P(x) = exp(a + b ln(x + c)) / (1 + exp(a + b ln(x + c))), x = distance
from the principal fault in km, 500 m cells, with the Table 2 coefficients
for the regular and conservative scenarios on the hanging wall (HW) and
footwall (FW). The model's rx > 0 selects the hanging wall.
"""
import numpy as np
import pytest

from openquake.fdha.secondary_surf_rup import FerrarioLivio2021SecondarySR

pytestmark = pytest.mark.unit

# Table 2: (scenario, side) -> (a, b, c)
TABLE_2 = {
    ("regular", "HW"): (-2.254, -1.175, 1e-5),
    ("regular", "FW"): (-3.459, -1.903, 1.008e-5),
    ("conservative", "HW"): (-1.888, -0.8802, 1.009e-5),
    ("conservative", "FW"): (-2.505, -1.181, 1.006e-5),
}
R_KM = np.array([0.0, 0.05, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0])


def _eq2(a, b, c, x):
    f = a + b * np.log(x + c)
    return np.exp(f) / (1.0 + np.exp(f))


@pytest.mark.parametrize("version,side", sorted(TABLE_2))
def test_matches_table_2(version, side):
    a, b, c = TABLE_2[(version, side)]
    rx = R_KM if side == "HW" else -R_KM
    rx = np.where(rx == 0.0, 1e-9 if side == "HW" else 0.0, rx)
    got = FerrarioLivio2021SecondarySR().get_prob(
        R_KM, rx, version=version)
    np.testing.assert_allclose(got, _eq2(a, b, c, R_KM), rtol=1e-12)


def test_hanging_wall_exceeds_footwall():
    """The paper's key contrast: a slower decay on the hanging wall."""
    m = FerrarioLivio2021SecondarySR()
    r = np.array([0.5, 1.0, 5.0])
    assert np.all(m.get_prob(r, r) > m.get_prob(r, -r))


def test_conservative_exceeds_regular_on_hanging_wall():
    m = FerrarioLivio2021SecondarySR()
    r = np.array([0.5, 1.0, 5.0, 10.0])
    assert np.all(m.get_prob(r, r, version="conservative")
                  > m.get_prob(r, r, version="regular"))


def test_pinned_version_equals_call_time_version():
    r = np.array([0.5, 2.0])
    np.testing.assert_array_equal(
        FerrarioLivio2021SecondarySR(version="conservative").get_prob(r, r),
        FerrarioLivio2021SecondarySR().get_prob(r, r,
                                                version="conservative"))
