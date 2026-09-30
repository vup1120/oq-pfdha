# -*- coding: utf-8 -*-
"""
Moss2013PrimarySR against its published regressions.

Moss et al. (2013), SRL 84(3), Eq. 1: P = f(z) = 1 / (1 + exp(-z)), with
z = beta + beta_1 Mw per faulting style and Vs30 bin (stiff > 600 m/s), as
printed in the legend of their Figure 3:

    reverse      z_stiff = -13.9745 + 2.1395 Mw   z_soft =  -6.2548 + 0.8308 Mw
    strike-slip  z_stiff = -11.4071 + 1.8465 Mw   z_soft = -12.2908 + 1.9520 Mw

The reverse pair is restated as Eqs 3.4-3.5 of GIRS-2022-05. That report's
redrawn Fig. 3.2 has a wrong reverse soft-soil curve; see benchmark/moss_2013.
"""
import numpy as np
import pytest

from openquake.fdha.primary_surf_rup import Moss2013PrimarySR

pytestmark = pytest.mark.unit

MAGS = np.linspace(4.2, 8.7, 46)   # data range stated in the paper

# Figure 3 legend: (style, Vs30 bin) -> (beta, beta_1)
FIG_3 = {
    ("reverse", "stiff"): (-13.9745, 2.1395),
    ("reverse", "soft"): (-6.2548, 0.8308),
    ("strike-slip", "stiff"): (-11.4071, 1.8465),
    ("strike-slip", "soft"): (-12.2908, 1.9520),
}
VS30 = {"stiff": (760.0, 601.0), "soft": (200.0, 599.0)}


def _eq1(z):
    return 1.0 / (1.0 + np.exp(-z))


@pytest.mark.parametrize("style,site", sorted(FIG_3))
def test_matches_fig_3_legend(style, site):
    a, b = FIG_3[(style, site)]
    for vs30 in VS30[site]:
        got = Moss2013PrimarySR().get_prob(MAGS, style=style, vs30=vs30)
        np.testing.assert_allclose(got, _eq1(a + b * MAGS), rtol=1e-12)


def test_stiff_soft_boundary_is_600_m_per_s():
    """600 m/s itself falls on the soft side (the code tests vs30 > 600)."""
    a, b = FIG_3[("reverse", "soft")]
    assert Moss2013PrimarySR().get_prob(6.5, "reverse", 600.0) == \
        pytest.approx(float(_eq1(a + b * 6.5)), rel=1e-12)


def test_pinned_parameters_equal_call_time_parameters():
    np.testing.assert_array_equal(
        Moss2013PrimarySR(style="strike-slip", vs30=300.0).get_prob(MAGS),
        Moss2013PrimarySR().get_prob(MAGS, style="strike-slip", vs30=300.0))


def test_normal_style_and_missing_vs30_raise():
    """The paper has too few normal-faulting events to regress."""
    with pytest.raises(ValueError):
        Moss2013PrimarySR().get_prob(6.5, style="normal", vs30=760.0)
    with pytest.raises(ValueError):
        Moss2013PrimarySR().get_prob(6.5, style="reverse")
