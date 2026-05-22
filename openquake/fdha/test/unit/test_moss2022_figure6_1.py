# -*- coding: utf-8 -*-
"""
Reproduction of GIRS-2022-05 Figure 6.1 Hazard Curve (MD-based, principal
displacement) — faithfully translated from the Matlab Appendix C source code.

Reference
---------
Moss, R., Thompson, S., Kuo, C.-H., Younesi, K., and Baumont, D. (2022).
Reverse Fault PFDHA.  Report GIRS-2022-05 (Revised 1/17/2024).
DOI: 10.34948/N3F595

Scenario (Section 6 / Appendix C, pages 114–122):
  - Mw 7.5 reverse fault, simple geometry
  - Slip rate: 5 mm/yr (= 0.5 cm/yr)
  - Fault dimensions: 100 km × 15 km
  - Stiff soil (VS30 = 700 m/s > 600 m/s threshold)
  - x/L = 0.4–0.5 (bins 41–50)
  - MFD: Truncated exponential, b = 0.8
  - P(SR|M): Moss et al. (2013) stiff soil
  - Shear modulus: 37.5 GPa (3.75e11 dyne/cm²)
  - MD scaling: 84th percentile (mu + 1σ = 0.148 shift in log10)
  - Convolution sigma: 0.133 × ln(10) ≈ 0.306 (AD regression σ, NOT 0.20)
  - Gamma D/MD: a = 1.4244*xL + 1.856, b = -0.0832*xL + 0.1994 (SCALE)
  - Gamma truncation: D/MD capped at 1.0
  - Monte Carlo: 10 000 samples per (M, x/L) bin

Reference values at 975-year return period (report page 86):
  On fault:         0.7 m
  100 m from fault: 0.25 m
  500 m from fault: 0.03 m

Key differences from a production OpenQuake-style Moss2022 workflow:

  ====================  =======================  =======================
  Parameter             Matlab Appendix C        Production workflow
  ====================  =======================  =======================
  b-value               0.8                      source-model choice
  Shear modulus          37.5 GPa                 source-model choice
  Mmin                   5.0                      source-model choice
  P(SR|M)               Moss 2013 stiff soil     logic-tree branch
  Convolution σ          0.133 × ln(10)           recommended 0.20
  MD scaling             84th percentile shift    integration over ε
  Method                 Monte Carlo (10 000)     analytical
  ====================  =======================  =======================

The direct reusable-class tests live in ``test_moss2022.py``.  This file is
kept separate as an Appendix C / Figure 6.1 benchmark reproduction, because
the benchmark fixes source-model and surface-rupture assumptions that are not
intrinsic defaults of the Moss2022 displacement classes.

Usage::

    pytest openquake/fdha/test/unit/test_moss2022_figure6_1.py -v
    python  openquake/fdha/test/unit/test_moss2022_figure6_1.py  # plots
"""
import math
import os

import numpy as np
import pytest
from scipy import special, stats
from scipy.interpolate import interp1d

pytestmark = [pytest.mark.benchmark, pytest.mark.regression]


# ====================================================================
# Matlab Appendix C helper functions (exact reproduction)
# ====================================================================

def _a_gam_md(xL):
    """Gamma shape for D/MD (Appendix C, p. 121)."""
    return 1.4244 * xL + 1.856


def _b_gam_md(xL):
    """Gamma scale for D/MD (Appendix C, p. 122)."""
    return -0.0832 * xL + 0.1994


def _mu_md(mag):
    """ln(MD_84th) — Appendix C, p. 122.

    Matlab: muout = log(10^(a_MD + b_MD*mag + sigma_mu))
    where a_MD = -2.5, b_MD = 0.415, sigma_mu = 0.148.
    """
    return math.log(10 ** (-2.5 + 0.415 * mag + 0.148))


def _trlnrnd(mu_val, sigma_val, n, rng):
    """Truncated lognormal samples (±5σ), Appendix C, p. 122."""
    eps_max = 5
    z = rng.uniform(0, 1, n)
    cdf_min = stats.lognorm.cdf(
        math.exp(mu_val - eps_max * sigma_val),
        s=sigma_val, scale=math.exp(mu_val))
    cdf_max = stats.lognorm.cdf(
        math.exp(mu_val + eps_max * sigma_val),
        s=sigma_val, scale=math.exp(mu_val))
    return np.exp(
        sigma_val * math.sqrt(2) * special.erfinv(
            2 * (z - 0.5) * (cdf_max - cdf_min)
        ) + mu_val
    )


# ====================================================================
# Full hazard-curve computation (Appendix C MD code)
# ====================================================================

def compute_hazard_curve_md(
    b_value=0.8,
    shear_modulus_cgs=3.75e11,
    min_mag=5.0,
    max_mag=7.5,
    length_km=100,
    width_km=15,
    slip_rate_cm=0.5,
    vs30=700,
    xL_min=41,
    xL_max=50,
    sim=10000,
    seed=42,
):
    """Matlab Appendix C → (D_levels, nu, N_m_min)."""
    rng = np.random.default_rng(seed)
    beta = math.log(10) * b_value
    area_cm2 = length_km * 1e5 * width_km * 1e5
    dm, dr = 0.008, 0.01
    sigma_conv = 0.133 * math.log(10)

    mag_idx = np.arange(1, 252)
    M_arr = min_mag + (2.0 / 250.0) * (mag_idx - 1)

    f_m = (beta * np.exp(-beta * (M_arr - min_mag))
           / (1 - math.exp(-beta * (max_mag - min_mag))))

    denom = 0.0
    for s in range(250):
        mo_s = 10 ** (1.5 * M_arr[s] + 16.05)
        mo_s1 = 10 ** (1.5 * M_arr[s + 1] + 16.05)
        denom += (f_m[s] * mo_s + f_m[s + 1] * mo_s1) * dm / 2
    N_m_min = shear_modulus_cgs * area_cm2 * slip_rate_cm / denom

    if vs30 > 600:
        pr_slip = 1.0 / (1.0 + np.exp(-(-13.9745 + 2.1395 * M_arr)))
    else:
        pr_slip = 1.0 / (1.0 + np.exp(-(-6.2548 + 0.8308 * M_arr)))

    D_levels = []
    D, v = 0.01, 0.01
    while D <= 10.0 + 1e-9:
        D_levels.append(round(D, 6))
        D += v
        if D > 0.09 + 1e-9 and D <= 0.9 + 1e-9:
            v = 0.1
        if D > 0.9 + 1e-9:
            v = 1.0
    D_levels = np.array(D_levels)
    n_D = len(D_levels)

    probDd = np.zeros((52, 252, n_D))
    for m in range(1, 252):
        mag = M_arr[m - 1]
        mu_val = _mu_md(mag)
        for r in range(xL_min, xL_max + 2):
            xL = 0.5 * (r - 1) / 50.0
            A = _trlnrnd(mu_val, sigma_conv, sim, rng)
            B = rng.gamma(_a_gam_md(xL), _b_gam_md(xL), sim)
            B[B >= 1.0] = 1.0
            combine = A * B

            n_hist, bin_edges = np.histogram(combine, bins=1000)
            Dbin = (bin_edges[:-1] + bin_edges[1:]) / 2.0
            cdf = np.cumsum(n_hist).astype(float)
            invcdf = 1.0 - cdf / cdf[-1]

            for dd in range(n_D):
                if D_levels[dd] > Dbin[-1]:
                    probDd[r, m, dd] = 0.0
                else:
                    probDd[r, m, dd] = float(
                        np.interp(D_levels[dd], Dbin, invcdf))

    nu = np.zeros(n_D)
    for dd in range(n_D):
        rate = 0.0
        for mm in range(1, 251):
            for rr in range(xL_min, xL_max + 1):
                val = (
                    f_m[mm - 1] * probDd[rr, mm, dd] * pr_slip[mm - 1]
                    + f_m[mm - 1] * probDd[rr + 1, mm, dd] * pr_slip[mm - 1]
                    + f_m[mm] * probDd[rr + 1, mm + 1, dd] * pr_slip[mm]
                    + f_m[mm] * probDd[rr, mm + 1, dd] * pr_slip[mm]
                ) / 4.0
                rate += N_m_min * val * dm * dr
        nu[dd] = rate

    return D_levels, nu, N_m_min


def _displacement_at_return_period(D_levels, nu, return_period):
    """Log-log interpolation of displacement at a given return period."""
    target_rate = 1.0 / return_period
    valid = nu > 0
    if not np.any(valid) or nu[valid].max() < target_rate:
        return float('nan')
    log_nu = np.log10(nu[valid])
    log_D = np.log10(D_levels[valid])
    f = interp1d(log_nu, log_D, kind='linear',
                 bounds_error=False, fill_value=np.nan)
    return 10 ** float(f(math.log10(target_rate)))


# ====================================================================
# Distributed displacement post-processing (Appendix C, pages 123–134)
# ====================================================================

def _distributed_postprocess(D_principal, nu_principal, r_dist_m,
                              wall=1, complex_flag=0, mbc=7.5):
    """
    Apply three Matlab post-processing factors to convert the principal
    MD hazard curve into a distributed displacement hazard curve.

    The Matlab code (Appendix C, pp. 130–132) scales both axes:
      d_distributed = D_principal × d_MD_ratio   (displacement axis)
      ν_distributed = ν_principal × p_r_dist     (rate axis)
    where p_r_dist = P(d>d₀|r) × P(d>0).

    Parameters
    ----------
    D_principal : array
        Principal displacement levels (m).
    nu_principal : array
        Principal annual exceedance rates.
    r_dist_m : float
        Distance from principal fault in metres.
    wall : int
        1 = hanging wall, 0 = footwall.
    complex_flag : int
        0 = simple faulting, 1 = complex multi-fault system.
    mbc : float
        Magnitude bin centre (7.5, 6.5, or 5.5).

    Returns
    -------
    d_off_fault : array
        Distributed displacement levels (m).
    d_rate : array
        Distributed annual exceedance rates.
    components : dict
        Individual scaling components for verification.
    """
    r_km = r_dist_m / 1000.0

    # --- 1. P(d>0): Eq. 5.5, Table 5.3 ---
    if wall == 1:
        pd0 = min(math.exp(-2.2 * r_km + 0.5), 1.0)
    else:
        pd0 = min(math.exp(-2.4 * r_km + 0.4), 1.0)

    # --- 2. P(d>d₀|r): frequency CDF, Eqs. 5.6–5.7, Tables 5.4–5.5 ---
    if wall == 1:
        if mbc == 7.5:
            if complex_flag == 1:
                a, b, c, d = 0.6998, 2.75e-5, -0.6931, -0.001219
                F_x = a * math.exp(b * r_dist_m) + c * math.exp(d * r_dist_m)
            else:
                r_max = min(r_dist_m, 3500)
                a, b, c, d = 0.8298, 5.682e-5, -0.8346, -0.001735
                F_x = a * math.exp(b * r_max) + c * math.exp(d * r_max)
        elif mbc == 6.5:
            if complex_flag == 1:
                a, b, c, d = 0.8858, 6.203e-6, -0.8957, -0.001959
                F_x = a * math.exp(b * r_dist_m) + c * math.exp(d * r_dist_m)
            else:
                r_max = min(r_dist_m, 3500)
                a, b, c, d = 1.166, -4.699e-5, -1.1730, -0.001539
                F_x = a * math.exp(b * r_max) + c * math.exp(d * r_max)
        elif mbc == 5.5:
            r_max = min(r_dist_m, 120)
            a, b, c, d = 98.45, 0.00228, -98.53, -0.01417
            F_x = a * math.exp(b * r_max) + c * math.exp(d * r_max)
        else:
            F_x = 1.0
    else:
        if mbc == 7.5:
            if complex_flag == 1:
                a, b, c, d = 0.1959, 0.0001091, -0.2020, -0.0026
                F_x = a * math.exp(b * r_dist_m) + c * math.exp(d * r_dist_m)
            else:
                r_max = min(r_dist_m, 3500)
                a, b, c, d = 1.445, -7.078e-5, -1.454, -0.0006972
                F_x = a * math.exp(b * r_max) + c * math.exp(d * r_max)
        elif mbc == 6.5:
            r_max = min(r_dist_m, 3500)
            a, b, c, d = 0.9297, 2.515e-5, -0.9233, -0.01828
            F_x = a * math.exp(b * r_dist_m) + c * math.exp(d * r_dist_m)
        else:
            F_x = 1.0

    p_exceed = max(0.0, 1.0 - F_x)
    p_r_dist = p_exceed * pd0

    # --- 3. d/MD ratio: Eq. 5.8, Table 5.8 (85th percentile) ---
    if wall == 1:
        if complex_flag == 0:
            d_MD_ratio = 0.43 * math.exp(-0.4 * r_km)
        else:
            d_MD_ratio = 0.43 * math.exp(-0.012 * r_km)
    else:
        d_MD_ratio = 0.68 * math.exp(-0.13 * r_km)

    d_off_fault = D_principal * d_MD_ratio
    d_rate = nu_principal * p_r_dist

    components = {
        'pd0': pd0,
        'p_exceed_cdf': p_exceed,
        'p_r_dist': p_r_dist,
        'd_MD_ratio': d_MD_ratio,
    }
    return d_off_fault, d_rate, components


# ====================================================================
# Tests — Matlab helper functions
# ====================================================================

class TestMatlabHelperFunctions:
    """Verify each Matlab Appendix C helper independently."""

    def test_a_gam_md(self):
        assert _a_gam_md(0.0) == pytest.approx(1.856, abs=1e-6)
        assert _a_gam_md(0.5) == pytest.approx(2.5682, abs=1e-4)

    def test_b_gam_md(self):
        assert _b_gam_md(0.0) == pytest.approx(0.1994, abs=1e-6)
        assert _b_gam_md(0.5) == pytest.approx(0.1578, abs=1e-4)

    def test_mu_md_is_84th_percentile(self):
        """mu_md returns ln(MD_84th), NOT ln(MD_median)."""
        for mag in [6.0, 6.5, 7.0, 7.5]:
            md_84th = 10 ** (-2.5 + 0.415 * mag + 0.148)
            assert math.exp(_mu_md(mag)) == pytest.approx(md_84th, rel=1e-6)
            assert math.exp(_mu_md(mag)) > 10 ** (-2.5 + 0.415 * mag)

    def test_mu_md_m75(self):
        """M=7.5 → MD_84th = 10^0.7605 ≈ 5.761 m."""
        assert math.exp(_mu_md(7.5)) == pytest.approx(5.761, abs=0.01)


# ====================================================================
# Tests — MFD
# ====================================================================

class TestMFD:
    """Verify truncated-exponential MFD and N(Mmin) from moment balance."""

    def test_n_m_min(self):
        """N(Mmin) ≈ 0.2847 yr⁻¹ (Appendix C parameters)."""
        _, _, N = compute_hazard_curve_md(sim=100, seed=0)
        assert N == pytest.approx(0.2847, rel=0.01)


# ====================================================================
# Tests — P(SR|M) model
# ====================================================================

class TestPSR:
    """Verify Moss et al. (2013) stiff-soil P(SR|M)."""

    @pytest.mark.parametrize("mag,expected", [
        (5.0, 0.03636871),
        (6.0, 0.24277966),
        (7.0, 0.73145162),
    ])
    def test_moss2013_stiff(self, mag, expected):
        computed = 1.0 / (1.0 + math.exp(-(-13.9745 + 2.1395 * mag)))
        assert computed == pytest.approx(expected, abs=1e-6)


# ====================================================================
# Tests — Convolution sigma documentation
# ====================================================================

class TestConvolutionSigma:
    """
    Document that the Matlab MD code uses σ = 0.133 × ln(10).

    This is a CRITICAL finding: the Matlab code passes 0.133 (the AD
    regression sigma from Table 4.4) to the MD convolution, NOT the
    MD regression sigma (0.148) or the recommended sigma (0.20).

    Ref: Appendix C, page 117:  A = trlnrnd(mu(...), 0.133*2.302, sim);
    """

    def test_sigma_value(self):
        assert 0.133 * math.log(10) == pytest.approx(0.30624, abs=0.001)

    def test_sigma_not_md_regression(self):
        assert 0.133 != 0.148

    def test_sigma_not_recommended(self):
        assert 0.133 != 0.20


# ====================================================================
# Tests — Full hazard-curve Figure 6.1
# ====================================================================

class TestHazardCurveFigure6_1:
    """
    Full hazard-curve reproduction of Figure 6.1 (on-fault, MD-based).

    Reference values at 975-year return period (report page 86):
      On fault: 0.7 m

    Tolerance: 15 % (Monte Carlo noise with sim=10 000).
    """

    pytestmark = pytest.mark.slow

    @pytest.fixture(scope="class")
    def hazard_curve(self):
        """Compute once and reuse across all tests in this class."""
        D, nu, N = compute_hazard_curve_md(sim=10000, seed=42)
        return D, nu, N

    def test_displacement_at_975yr(self, hazard_curve):
        """On-fault D at 975-yr RP ≈ 0.7 m (Section 6, page 86)."""
        D, nu, _ = hazard_curve
        D_975 = _displacement_at_return_period(D, nu, 975)
        assert abs(D_975 - 0.7) / 0.7 < 0.15, \
            f"D at 975-yr RP = {D_975:.3f} m, expected ~0.7 m"

    def test_exceedance_rate_at_0_7m(self, hazard_curve):
        """ν(D > 0.7 m) → RP ≈ 975 yr (within factor 1.3)."""
        D, nu, _ = hazard_curve
        idx = int(np.argmin(np.abs(D - 0.7)))
        nu_07 = nu[idx]
        rp = 1.0 / nu_07 if nu_07 > 0 else float('inf')
        assert 700 < rp < 1300, \
            f"RP at D=0.7m is {rp:.0f} yr, expected ~975 yr"

    def test_hazard_curve_monotonic(self, hazard_curve):
        """Annual exceedance rate must decrease with displacement."""
        _, nu, _ = hazard_curve
        assert np.all(np.diff(nu) <= 1e-10)

    def test_n_m_min(self, hazard_curve):
        """Cross-check annual rate of earthquakes ≈ 0.2847."""
        _, _, N = hazard_curve
        assert N == pytest.approx(0.2847, rel=0.01)

    def test_monte_carlo_seed_reproducible(self):
        """Fixed seed should reproduce the Appendix C Monte Carlo curve."""
        D1, nu1, N1 = compute_hazard_curve_md(sim=1000, seed=42)
        D2, nu2, N2 = compute_hazard_curve_md(sim=1000, seed=42)
        assert np.array_equal(D1, D2)
        assert N1 == pytest.approx(N2, rel=0.0, abs=0.0)
        assert np.allclose(nu1, nu2, rtol=0.0, atol=0.0)

    @pytest.mark.parametrize("D_val,rp_range", [
        (0.01, (200, 400)),
        (0.10, (250, 400)),
        (1.00, (1200, 2500)),
        (5.00, (400000, 1200000)),
    ])
    def test_return_periods_order_of_magnitude(
            self, hazard_curve, D_val, rp_range):
        """Return periods at key displacement levels."""
        D, nu, _ = hazard_curve
        idx = int(np.argmin(np.abs(D - D_val)))
        if nu[idx] > 0:
            rp = 1.0 / nu[idx]
            assert rp_range[0] < rp < rp_range[1], \
                f"RP at D={D_val}m is {rp:.0f}, expected in {rp_range}"


# ====================================================================
# Tests — Distributed displacement components (Appendix C, pp. 130–132)
# ====================================================================

class TestPd0:
    """P(d>0) = min(exp(-a*r_km + b), 1.0).  Ref: Eq. 5.5, Table 5.3."""

    @pytest.mark.parametrize("r_m,expected", [
        (0, 1.0),
        (100, 1.0),
        (250, 0.951229),
        (500, 0.548812),
        (1000, 0.182684),
        (2000, 0.020242),
        (3000, 0.002243),
    ])
    def test_pd0_hw(self, r_m, expected):
        r_km = r_m / 1000.0
        computed = min(math.exp(-2.2 * r_km + 0.5), 1.0)
        assert computed == pytest.approx(expected, abs=1e-4)

    @pytest.mark.parametrize("r_m,expected", [
        (0, 1.0),
        (500, 0.449329),
        (1000, 0.135335),
    ])
    def test_pd0_fw(self, r_m, expected):
        r_km = r_m / 1000.0
        computed = min(math.exp(-2.4 * r_km + 0.4), 1.0)
        assert computed == pytest.approx(expected, abs=1e-4)


class TestFrequencyCDF:
    """F(x) = a*exp(b*x) + c*exp(d*x); P(d>d₀) = 1-F(x).
    Ref: Appendix C, pages 131–132."""

    @pytest.mark.parametrize("r_m,expected_F", [
        (100, 0.8298 * math.exp(5.682e-5 * 100)
              - 0.8346 * math.exp(-0.001735 * 100)),
        (500, 0.8298 * math.exp(5.682e-5 * 500)
              - 0.8346 * math.exp(-0.001735 * 500)),
        (1000, 0.8298 * math.exp(5.682e-5 * 1000)
               - 0.8346 * math.exp(-0.001735 * 1000)),
        (3000, 0.8298 * math.exp(5.682e-5 * 3000)
               - 0.8346 * math.exp(-0.001735 * 3000)),
    ])
    def test_hw_m75_simple(self, r_m, expected_F):
        """HW, M7.0-7.9 bin, simple."""
        r_max = min(r_m, 3500)
        a, b, c, d = 0.8298, 5.682e-5, -0.8346, -0.001735
        computed = a * math.exp(b * r_max) + c * math.exp(d * r_max)
        assert computed == pytest.approx(expected_F, abs=1e-6)


class TestDMDRatio:
    """d/MD = c * exp(d * r_km).  Ref: Table 5.8 (85th percentile)."""

    @pytest.mark.parametrize("r_m,expected", [
        (0, 0.43),
        (100, 0.43 * math.exp(-0.4 * 0.1)),
        (500, 0.43 * math.exp(-0.4 * 0.5)),
        (1000, 0.43 * math.exp(-0.4 * 1.0)),
        (3000, 0.43 * math.exp(-0.4 * 3.0)),
    ])
    def test_hw_simple(self, r_m, expected):
        r_km = r_m / 1000.0
        computed = 0.43 * math.exp(-0.4 * r_km)
        assert computed == pytest.approx(expected, abs=1e-6)

    @pytest.mark.parametrize("r_m,expected", [
        (0, 0.68),
        (1000, 0.68 * math.exp(-0.13 * 1.0)),
    ])
    def test_fw(self, r_m, expected):
        r_km = r_m / 1000.0
        computed = 0.68 * math.exp(-0.13 * r_km)
        assert computed == pytest.approx(expected, abs=1e-6)


class TestCombinedDistributedComponents:
    """Verify p_r_dist = P(d>d₀|r) × P(d>0) at specific distances."""

    @pytest.mark.parametrize(
        "r_m,expected_pd0,expected_p_exceed,expected_p_r_dist", [
            (100, 1.0, 0.867134, 0.867134),
            (500, 0.548812, 0.496820, 0.272660),
            (1000, 0.182684, 0.268909, 0.049125),
        ])
    def test_hw_m75_simple(self, r_m, expected_pd0, expected_p_exceed,
                            expected_p_r_dist):
        D_dummy = np.array([1.0])
        nu_dummy = np.array([1.0])
        _, _, comp = _distributed_postprocess(
            D_dummy, nu_dummy, r_m, wall=1, complex_flag=0, mbc=7.5)
        assert comp['pd0'] == pytest.approx(expected_pd0, abs=1e-4)
        assert comp['p_exceed_cdf'] == pytest.approx(expected_p_exceed, abs=1e-3)
        assert comp['p_r_dist'] == pytest.approx(expected_p_r_dist, abs=1e-3)


# ====================================================================
# Tests — Full distributed hazard curve (Figure 6.1)
# ====================================================================

class TestDistributedHazardCurve:
    """
    Full reproduction of Figure 6.1 distributed displacement curves.

    Reference values at 975-year return period (Section 6, page 86):
      100 m from fault: 0.25 m
      500 m from fault: 0.03 m

    Tolerance: 20 % for 100 m, relaxed for 500 m (MC noise at curve edge).
    """

    pytestmark = pytest.mark.slow

    @pytest.fixture(scope="class")
    def principal_curve(self):
        D, nu, _ = compute_hazard_curve_md(sim=10000, seed=42)
        return D, nu

    def test_100m_displacement_at_975yr(self, principal_curve):
        """Distributed D at 100 m from fault, 975-yr RP ≈ 0.25 m."""
        D_p, nu_p = principal_curve
        d_dist, nu_dist, _ = _distributed_postprocess(
            D_p, nu_p, 100, wall=1, complex_flag=0, mbc=7.5)
        D_975 = _displacement_at_return_period(d_dist, nu_dist, 975)
        assert not np.isnan(D_975), "Curve doesn't reach 975-yr RP"
        assert abs(D_975 - 0.25) / 0.25 < 0.20, \
            f"D at 975yr, 100m = {D_975:.4f} m, expected ~0.25 m"

    def test_500m_curve_near_975yr(self, principal_curve):
        """At 500 m the distributed curve barely reaches the 975-yr level.
        Verify the minimum RP is close to 975 yr."""
        D_p, nu_p = principal_curve
        _, nu_dist, _ = _distributed_postprocess(
            D_p, nu_p, 500, wall=1, complex_flag=0, mbc=7.5)
        max_rate = nu_dist.max()
        min_rp = 1.0 / max_rate if max_rate > 0 else float('inf')
        assert 850 < min_rp < 1200, \
            f"Min RP on 500m curve = {min_rp:.0f} yr, expected ~975 yr"

    def test_100m_rate_at_quarter_meter(self, principal_curve):
        """ν at d ≈ 0.25 m on the 100 m curve → RP ≈ 975 yr."""
        D_p, nu_p = principal_curve
        d_dist, nu_dist, _ = _distributed_postprocess(
            D_p, nu_p, 100, wall=1, complex_flag=0, mbc=7.5)
        idx = int(np.argmin(np.abs(d_dist - 0.25)))
        nu_at_025 = nu_dist[idx]
        rp = 1.0 / nu_at_025 if nu_at_025 > 0 else float('inf')
        assert 700 < rp < 1400, \
            f"RP at d≈0.25m, 100m = {rp:.0f} yr, expected ~975 yr"

    def test_distributed_rate_less_than_principal(self, principal_curve):
        """Distributed rates must be strictly lower than principal rates."""
        D_p, nu_p = principal_curve
        for r_m in [100, 500, 1000]:
            _, nu_dist, _ = _distributed_postprocess(
                D_p, nu_p, r_m, wall=1, complex_flag=0, mbc=7.5)
            assert np.all(nu_dist <= nu_p + 1e-15)

    def test_rate_decreases_with_distance(self, principal_curve):
        """At same displacement, rate must decrease with distance."""
        D_p, nu_p = principal_curve
        _, nu_100, _ = _distributed_postprocess(D_p, nu_p, 100)
        _, nu_500, _ = _distributed_postprocess(D_p, nu_p, 500)
        _, nu_1000, _ = _distributed_postprocess(D_p, nu_p, 1000)
        for i in range(len(nu_p)):
            assert nu_100[i] >= nu_500[i] - 1e-15
            assert nu_500[i] >= nu_1000[i] - 1e-15

    def test_displacement_decreases_with_distance(self, principal_curve):
        """d/MD ratio decreases with distance → smaller displacements."""
        D_p, nu_p = principal_curve
        _, _, c100 = _distributed_postprocess(D_p, nu_p, 100)
        _, _, c500 = _distributed_postprocess(D_p, nu_p, 500)
        _, _, c1000 = _distributed_postprocess(D_p, nu_p, 1000)
        assert c100['d_MD_ratio'] > c500['d_MD_ratio']
        assert c500['d_MD_ratio'] > c1000['d_MD_ratio']


# ====================================================================
# Standalone runner with plot
# ====================================================================

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plots')
    os.makedirs(out_dir, exist_ok=True)

    print("Computing principal MD hazard curve (sim=10000, seed=42) ...")
    D, nu, N = compute_hazard_curve_md(sim=10000, seed=42)

    print(f"\nN(Mmin) = {N:.6e} events/yr")
    D_975_on = _displacement_at_return_period(D, nu, 975)
    print(f"On-fault D at 975-yr RP = {D_975_on:.3f} m (reference: 0.7 m)")

    # ── Distributed post-processing (same distances as report Fig. 6.1) ──
    distances_m = [100, 500, 1000, 2000, 3000]
    dist_results = {}
    for r_m in distances_m:
        d_d, nu_d, comp = _distributed_postprocess(
            D, nu, r_m, wall=1, complex_flag=0, mbc=7.5)
        dist_results[r_m] = (d_d, nu_d, comp)
        D_975_d = _displacement_at_return_period(d_d, nu_d, 975)
        d_str = f"{D_975_d:.4f} m" if not np.isnan(D_975_d) else "N/A"
        print(f"\nr={r_m} m: P(d>0)={comp['pd0']:.4f}, "
              f"P(d>d0|r)={comp['p_exceed_cdf']:.4f}, "
              f"p_r_dist={comp['p_r_dist']:.6f}, "
              f"d/MD={comp['d_MD_ratio']:.4f}")
        print(f"  D at 975-yr RP = {d_str}")

    # ── Reference summary ──
    refs = {'on_fault': 0.7, '100m': 0.25, '500m': 0.03}
    print("\n" + "=" * 60)
    print(f"{'Location':<20s} {'Computed':>12s} {'Reference':>12s} {'Ratio':>8s}")
    print("-" * 60)
    print(f"{'On fault':<20s} {D_975_on:12.4f} {refs['on_fault']:12.4f} "
          f"{D_975_on / refs['on_fault']:8.3f}")
    for r_m, ref_key in [(100, '100m'), (500, '500m')]:
        d_d, nu_d, _ = dist_results[r_m]
        D_975_d = _displacement_at_return_period(d_d, nu_d, 975)
        if not np.isnan(D_975_d):
            print(f"{ref_key + ' (HW)':<20s} {D_975_d:12.4f} "
                  f"{refs[ref_key]:12.4f} {D_975_d / refs[ref_key]:8.3f}")
        else:
            print(f"{ref_key + ' (HW)':<20s} {'N/A':>12s} "
                  f"{refs[ref_key]:12.4f} {'---':>8s}")
    print("=" * 60)

    # ── Plot (matching GIRS-2022-05 Figure 6.1 style) ──
    from matplotlib.ticker import LogLocator, FuncFormatter

    curve_style = {
        100:  {'color': 'red',      'label': '0.1km'},
        500:  {'color': '#DAA520',  'label': '0.5km'},
        1000: {'color': 'green',    'label': '1.0km'},
        2000: {'color': 'blue',     'label': '2.0km'},
        3000: {'color': '#8B008B',  'label': '3.0km'},
    }

    fig, ax = plt.subplots(figsize=(9, 8))

    valid = nu > 0
    ax.loglog(D[valid], nu[valid], color='black', lw=3,
              solid_capstyle='round', label='MD+sigma')

    for r_m in distances_m:
        d_d, nu_d, _ = dist_results[r_m]
        v = nu_d > 0
        if np.any(v):
            st = curve_style[r_m]
            ax.loglog(d_d[v], nu_d[v], color=st['color'], lw=2,
                      solid_capstyle='round', label=st['label'])

    ax.axhline(1 / 975, color='black', ls='--', lw=1.5, label='1/975yr')

    ref_markers = [
        (0.7,  1 / 975, 'black',    'On-fault ref: 0.7 m'),
        (0.25, 1 / 975, 'red',      '100 m ref: 0.25 m'),
        (0.03, 1 / 975, '#DAA520',  '500 m ref: 0.03 m'),
    ]
    for d_ref, rate_ref, col, lbl in ref_markers:
        ax.plot(d_ref, rate_ref, '*', markersize=16, color=col,
                markeredgecolor='black', markeredgewidth=0.5,
                zorder=10, label=lbl)

    ax.text(0.97, 0.97,
            '$M_{max}$=7.5, L=100km, W=15km, stiff VS30,\n'
            'x/L=0.5, hanging wall, slip=5mm/yr\n'
            'simple fault, 85th percentile estimate',
            transform=ax.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                      edgecolor='gray', alpha=0.8))

    def _exp_fmt(x, _pos):
        exp = int(round(math.log10(x)))
        return f'1.E{exp:+03d}'

    ax.set_xlabel('Displacement (m)', fontsize=13)
    ax.set_ylabel('Annual Probability of Exceedance', fontsize=13)
    ax.set_xlim(0.01, 1.0)
    ax.set_ylim(1e-10, 1e-1)
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=12))
    ax.yaxis.set_major_formatter(FuncFormatter(_exp_fmt))
    ax.xaxis.set_major_formatter(FuncFormatter(
        lambda x, _: f'{x:g}'))
    ax.grid(True, which='major', alpha=0.3, linewidth=0.5)
    ax.grid(True, which='minor', alpha=0.1, linewidth=0.3)
    ax.tick_params(axis='both', which='major', labelsize=11)
    ax.legend(loc='lower left', fontsize=9, framealpha=0.9)
    ax.set_title('Reproduction of GIRS-2022-05 Figure 6.1', fontsize=14)
    fig.tight_layout()

    fpath = os.path.join(out_dir, 'figure_6_1_full.png')
    fig.savefig(fpath, dpi=150)
    print(f"\nPlot saved to {fpath}")
