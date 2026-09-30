# -*- coding: utf-8 -*-
"""
Youngs2003SecondarySR and Youngs2003SecondaryFD against the published model.

Youngs et al. (2003), Earthquake Spectra 19(1), Appendix (checked against
the accepted manuscript):

* occurrence, 0.5 km x 0.5 km cells, P = exp(f) / (1 + exp(f)) with
  Eq. 7:  f = 2.06 + (-4.62 + 0.118 m + 0.682 h) ln(r + 3.32)
  Eq. 8:  f = 3.27 + (-8.28 + 0.577 m + 0.629 h) ln(r + 4.14) + 0.611 z_i,
          z_i = 0 for the average event;
  r in km, h = 1 on the hanging wall and 0 on the footwall;
* displacement: the 85th-95th percentile of D_distributed / MD_principal
  (Fig. 11) is 0.35 exp(-0.091 r) on the hanging wall and
  0.16 exp(-0.137 r) on the footwall; the ratio is gamma distributed with
  shape a = 2.5, whose 95th and 85th percentiles lie at x/b = 5.535 and
  4.058; it is convolved with the Wells and Coppersmith (1994, Table 2B,
  normal faulting) MD distribution log10 MD = -5.90 + 0.89 M, sigma 0.38.
  Truncating that distribution at +/-3 sigma is the library's choice.

The library's version "3" (the 50/50 average of Eqs 7 and 8) is a
convenience, not a published model.
"""
import numpy as np
import pytest
from scipy import integrate, stats

from openquake.fdha.secondary_surf_displ import Youngs2003SecondaryFD
from openquake.fdha.secondary_surf_rup import Youngs2003SecondarySR

pytestmark = pytest.mark.unit

R_KM = np.array([0.0, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0])
MAGS = (5.5, 6.5, 7.4)
GAMMA_SHAPE = 2.5
PERCENTILE_X_OVER_B = {"95": 5.535, "85": 4.058}


def _logistic(f):
    return np.exp(f) / (1.0 + np.exp(f))


def _eq7(m, r, h):
    return _logistic(2.06 + (-4.62 + 0.118 * m + 0.682 * h) * np.log(r + 3.32))


def _eq8(m, r, h):
    return _logistic(3.27 + (-8.28 + 0.577 * m + 0.629 * h) * np.log(r + 4.14))


def _fig11(r, hanging_wall):
    return (0.35 * np.exp(-0.091 * r) if hanging_wall
            else 0.16 * np.exp(-0.137 * r))


# ---------------------------------------------------------------------------
# Distributed occurrence
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("mag", MAGS)
@pytest.mark.parametrize("version,eq", [("1", _eq7), ("2", _eq8)])
def test_occurrence_matches_appendix_equations(mag, version, eq):
    model = Youngs2003SecondarySR()
    hw = model.get_prob(mag, R_KM + 1e-9, R_KM, version=version)
    fw = model.get_prob(mag, -R_KM, R_KM, version=version)
    np.testing.assert_allclose(hw, eq(mag, R_KM, 1), rtol=1e-12)
    np.testing.assert_allclose(fw, eq(mag, R_KM, 0), rtol=1e-12)


@pytest.mark.parametrize("mag", MAGS)
def test_version_3_is_the_equal_weight_average(mag):
    model = Youngs2003SecondarySR()
    rx = np.concatenate([R_KM + 1e-9, -R_KM])
    r = np.concatenate([R_KM, R_KM])
    avg = 0.5 * model.get_prob(mag, rx, r, version="1") \
        + 0.5 * model.get_prob(mag, rx, r, version="2")
    np.testing.assert_allclose(
        model.get_prob(mag, rx, r, version="3"), avg, rtol=1e-12)


def test_hanging_wall_exceeds_footwall():
    """The trend the regression was built to capture (paper, Fig. 9)."""
    model = Youngs2003SecondarySR()
    r = np.array([0.5, 2.0, 10.0])
    for version in ("1", "2"):
        assert np.all(model.get_prob(6.5, r, r, version=version)
                      > model.get_prob(6.5, -r, r, version=version))


# ---------------------------------------------------------------------------
# Distributed displacement
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("pct,q", [("95", 0.95), ("85", 0.85)])
def test_published_percentile_constants(pct, q):
    """x/b = 5.535 and 4.058 are the 95th and 85th percentiles of a
    gamma distribution with shape 2.5 (to the four digits printed)."""
    assert stats.gamma.ppf(q, GAMMA_SHAPE) == pytest.approx(
        PERCENTILE_X_OVER_B[pct], abs=5e-4)
    assert Youngs2003SecondaryFD._PERCENTILE_SCALING[pct] == \
        PERCENTILE_X_OVER_B[pct]


@pytest.mark.parametrize("pct", ["85", "95"])
def test_ratio_distribution_matches_fig_11_anchoring(pct):
    """get_prob_D_MD is the gamma survival function whose scale is the
    Fig. 11 curve divided by the percentile constant."""
    model = Youngs2003SecondaryFD()
    ratio = np.array([0.01, 0.05, 0.1, 0.2, 0.4, 0.8])
    for hanging_wall in (True, False):
        rx = R_KM + 1e-9 if hanging_wall else -R_KM
        got = model.get_prob_D_MD(ratio, rx, R_KM, percentile=pct)
        b = _fig11(R_KM, hanging_wall) / PERCENTILE_X_OVER_B[pct]
        expected = stats.gamma.sf(ratio[None, :], GAMMA_SHAPE,
                                  scale=b[:, None])
        np.testing.assert_allclose(got, expected, rtol=1e-12)


def test_md_distribution_is_wc94_normal_faulting():
    assert (Youngs2003SecondaryFD._WC94_MD_INTERCEPT,
            Youngs2003SecondaryFD._WC94_MD_SLOPE,
            Youngs2003SecondaryFD._WC94_MD_SIGMA) == (-5.90, 0.89, 0.38)


def _convolution(d0, mag, r, hanging_wall, pct, trunc=3.0):
    """P(d > d0 | m, r) by adaptive quadrature over log10 MD."""
    mu, sigma = -5.90 + 0.89 * mag, 0.38
    b = _fig11(r, hanging_wall) / PERCENTILE_X_OVER_B[pct]
    mass = stats.norm.cdf(trunc) - stats.norm.cdf(-trunc)

    def integrand(y):
        return (stats.gamma.sf(d0 / 10.0 ** y, GAMMA_SHAPE, scale=b)
                * stats.norm.pdf(y, mu, sigma) / mass)
    return integrate.quad(integrand, mu - trunc * sigma, mu + trunc * sigma,
                          epsabs=1e-16, epsrel=1e-11, limit=400)[0]


@pytest.fixture(scope="module")
def convolution_table():
    """(library, quadrature) pairs over magnitude, distance, wall side,
    percentile and displacement, computed once for both tests below."""
    model = Youngs2003SecondaryFD()
    d = np.logspace(-3, 1, 25)
    got, ref = [], []
    for mag in MAGS:
        for r in (0.2, 1.0, 5.0, 15.0):
            for hanging_wall in (True, False):
                for pct in ("85", "95"):
                    rx = r if hanging_wall else -r
                    got.append(model.get_prob(d, mag, np.array([rx]),
                                              np.array([r]),
                                              percentile=pct)[0])
                    ref.append([_convolution(di, mag, r, hanging_wall, pct)
                                for di in d])
    return np.concatenate(got), np.concatenate(ref)


def _worst(table, p_min, p_max):
    got, ref = table
    sel = (ref >= p_min) & (ref < p_max)
    return np.max(np.abs(got[sel] - ref[sel]) / ref[sel])


def test_convolution_matches_quadrature_above_1_percent(convolution_table):
    """Observed: <= 0.8 % for conditional probabilities >= 0.01."""
    assert _worst(convolution_table, 1e-2, 1.1) < 0.01


@pytest.mark.xfail(strict=True, reason=(
    "the library sums the MD distribution on 100 points with full weight "
    "at both ends of the +/-3 sigma range; the tail, controlled by the upper "
    "end, is over-predicted (about +3 % at P ~ 1e-3, +18 % at P ~ 1e-5). "
    "Trapezoid end weights bring it within 0.1 %"))
def test_convolution_tail_matches_quadrature(convolution_table):
    assert _worst(convolution_table, 1e-6, 1e-2) < 0.01
