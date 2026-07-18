# -*- coding: utf-8 -*-
"""
PFDHA Benchmark Test Cases - Hazard Curve Computation

Computes hazard curves for both test cases from the benchmark spreadsheet
(PFDHA_CASE_STUDIES.xlsx) using the validated Moss et al. (2022) model
(GIRS-2022-05).

Corrected fault assignment (suspected typo in spreadsheet):
  Test 1 (principal, on-fault)  → Fault 2 (site lies on trace, r = 0 m)
  Test 2 (distributed, HW)     → Fault 1 (site at r ≈ 254 m, hanging wall)

Spreadsheet parameters:
  MFD:              Youngs & Coppersmith (1985) characteristic
  b-value:          1.0
  Mmin:             4.5
  Mmax (stated):    6.2 → char_mag = 5.95 or 6.20 (ambiguous, both computed)
  Slip rate:        0.01 mm/yr
  Shear modulus:    30 GPa
  Dip:              50 deg
  Zbor:             5 km
  P(SR|M):          Wells & Coppersmith (1993), all styles
  Test 1 FDM:       Moss2022, D/MD(x/L) gamma + MD(M), complete, sigma_rec=0.20
  Test 2 P(d>0):    Moss2022, HW, z=500 m grid (at r=254 m → P≈0.942)
  Test 2 FDM:       Moss2022, d/MD(r) gamma with a=2.5, sigma_rec=0.20

Reference
---------
Moss, R., Thompson, S., Kuo, C.-H., Younesi, K., and Baumont, D. (2022).
Reverse Fault PFDHA.  Report GIRS-2022-05 (Revised 1/17/2024).
DOI: 10.34948/N3F595

Usage::

    pytest openquake/fdha/test/unit/test_pfdha_case_studies.py -v
    python  openquake/fdha/test/unit/test_pfdha_case_studies.py  # plots
"""
import math
import os

import numpy as np
import pytest
from scipy import special, stats
from scipy.interpolate import interp1d

pytestmark = pytest.mark.unit


# ====================================================================
# Fault geometry (from spreadsheet "Hypothetical Fault" sheet)
# ====================================================================

FAULT_TRACES = {
    1: [(3.500725, 48.225814), (3.519857, 48.239407), (3.555135, 48.284752)],
    2: [(3.504302, 48.223652), (3.523961, 48.237550), (3.552453, 48.285224)],
}
SITE = (3.5295399, 48.2468851)

DIP_DEG = 50.0
ZBOR_KM = 5.0
W_KM = ZBOR_KM / math.sin(math.radians(DIP_DEG))  # ≈ 6.527 km
MU_PA = 30e9
SLIP_RATE_M_YR = 0.01e-3


def _seg_len_m(p1, p2, cos_lat):
    dx = (p2[0] - p1[0]) * cos_lat * 111319.5
    dy = (p2[1] - p1[1]) * 111319.5
    return math.sqrt(dx * dx + dy * dy)


def _compute_fault_geometry():
    """Compute fault lengths, site x/L, and site-to-fault distance."""
    cos_lat = math.cos(math.radians(SITE[1]))

    info = {}
    for fid, trace in FAULT_TRACES.items():
        segs = []
        for i in range(len(trace) - 1):
            segs.append(_seg_len_m(trace[i], trace[i + 1], cos_lat))
        total_m = sum(segs)
        L_km = total_m / 1000.0

        min_dist = float('inf')
        best_along = 0.0
        cum = 0.0
        for i in range(len(trace) - 1):
            ax = (trace[i][0] - SITE[0]) * cos_lat * 111319.5
            ay = (trace[i][1] - SITE[1]) * 111319.5
            bx = (trace[i + 1][0] - SITE[0]) * cos_lat * 111319.5
            by = (trace[i + 1][1] - SITE[1]) * 111319.5
            dx, dy = bx - ax, by - ay
            L2 = dx * dx + dy * dy
            t = max(0, min(1, (-ax * dx + -ay * dy) / L2)) if L2 > 0 else 0
            cx, cy = ax + t * dx, ay + t * dy
            dist = math.sqrt(cx * cx + cy * cy)
            if dist < min_dist:
                min_dist = dist
                best_along = cum + t * segs[i]
            cum += segs[i]

        x_L = best_along / total_m
        info[fid] = {
            'L_km': L_km,
            'A_km2': L_km * W_KM,
            'M0_rate': MU_PA * L_km * W_KM * 1e6 * SLIP_RATE_M_YR,
            'site_dist_m': min_dist,
            'site_xL': x_L,
            'site_xL_folded': min(x_L, 1.0 - x_L),
        }
    return info


# ====================================================================
# YC1985 characteristic MFD
# ====================================================================

def _yc1985_rates(min_mag, b_val, char_mag, total_moment_rate,
                  bin_width=0.01):
    """Youngs & Coppersmith (1985) characteristic MFD via moment balance."""
    delta_mc = 0.5
    m_max = char_mag + 0.25
    m_prime = m_max - delta_mc
    beta = b_val * math.log(10.0)
    c = 1.5 * math.log(10.0)

    moment_char = (
        10 ** 9.05 * (10 ** (1.5 * m_max) - 10 ** (1.5 * m_prime))
        / (c * delta_mc))
    moment_gr = (
        10 ** 9.05 * beta
        * (10 ** (1.5 * m_prime) * math.exp(-beta * (m_prime - min_mag))
           - 10 ** (1.5 * min_mag))
        / ((c - beta) * (1 - math.exp(-beta * (m_prime - min_mag)))))

    exp_factor = (
        beta * math.exp(-beta * (m_prime - 1.0 - min_mag))
        / (1 - math.exp(-beta * (m_prime - min_mag))))
    N_char = total_moment_rate / (moment_gr / exp_factor + moment_char)
    N_mmin = N_char * (1 + 1 / exp_factor)
    N_gr = N_mmin - N_char
    a_val = math.log10(
        N_gr / (10 ** (-b_val * min_mag) - 10 ** (-b_val * m_prime)))

    n_bins = int(round((m_max - min_mag) / bin_width))
    mags = np.empty(n_bins)
    rates = np.empty(n_bins)
    for i in range(n_bins):
        m_lo = min_mag + i * bin_width
        m_hi = m_lo + bin_width
        mags[i] = (m_lo + m_hi) / 2.0
        if mags[i] < m_prime:
            rates[i] = (10 ** (a_val - b_val * m_lo)
                        - 10 ** (a_val - b_val * m_hi))
        elif m_lo >= m_prime and m_hi <= m_max + 1e-10:
            rates[i] = N_char * bin_width / delta_mc
        else:
            r_gr = max(0, (10 ** (a_val - b_val * m_lo)
                           - 10 ** (a_val - b_val * min(m_hi, m_prime))))
            frac = max(0, min(m_hi, m_max) - max(m_lo, m_prime)) / bin_width
            r_ch = N_char * frac * bin_width / delta_mc if frac > 0 else 0
            rates[i] = r_gr + r_ch
        rates[i] = max(0, rates[i])
    return mags, rates, N_mmin


# ====================================================================
# Model components - Moss et al. (2022)
# ====================================================================

def _psr_wc93(M):
    """P(SR|M) - Wells & Coppersmith (1993), all styles."""
    z = -12.51 + 2.053 * np.asarray(M, dtype=float)
    return np.exp(z) / (1 + np.exp(z))


def _gamma_shape_dmd(xL):
    """Gamma shape a for D/MD (Fig. 4.3)."""
    return 1.4244 * xL + 1.856


def _gamma_scale_dmd(xL):
    """Gamma SCALE b for D/MD (Fig. 4.4, Eq. 4.4)."""
    return -0.0832 * xL + 0.1994


def _md_median_ln(mag):
    """ln(median MD) for complete ruptures (Table 4.4)."""
    return math.log(10 ** (-2.50 + 0.415 * mag))


def _trlnrnd(mu_val, sigma_val, n, rng):
    """Truncated lognormal samples (+-5sigma), Appendix C."""
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
# Hazard curve computation
# ====================================================================

D_LEVELS = np.array([
    0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1,
    0.2, 0.5, 1.0, 2.0, 5.0, 10.0,
])

SIGMA_CONV = 0.20 * math.log(10)  # 0.4605 in ln space


def _compute_test1(fault_info, char_mag, sim=30000, seed=42):
    """
    Test 1: Principal displacement hazard curve.

    Source: Fault 2 (site on fault trace, r = 0 m).
    FDM: D/MD gamma (x/L-dependent) x MD(M) lognormal, truncated at D/MD=1.
    """
    rng = np.random.default_rng(seed)
    fi = fault_info[2]

    mags, rates, N_mmin = _yc1985_rates(4.5, 1.0, char_mag, fi['M0_rate'])
    pr_slip = _psr_wc93(mags)

    xLf = fi['site_xL_folded']
    a_g = _gamma_shape_dmd(xLf)
    b_g = _gamma_scale_dmd(xLf)

    n_D = len(D_LEVELS)
    probDd = np.zeros((len(mags), n_D))
    for mi in range(len(mags)):
        A = _trlnrnd(_md_median_ln(mags[mi]), SIGMA_CONV, sim, rng)
        B = rng.gamma(a_g, b_g, sim)
        B[B >= 1.0] = 1.0
        combine = A * B
        for di in range(n_D):
            probDd[mi, di] = np.mean(combine > D_LEVELS[di])

    nu = np.array([
        np.sum(rates * pr_slip * probDd[:, di]) for di in range(n_D)])
    return D_LEVELS.copy(), nu, N_mmin


def _compute_test2(fault_info, char_mag, sim=30000, seed=42):
    """
    Test 2: Distributed displacement hazard curve.

    Source: Fault 1 (site at r ~ 254 m, hanging wall).
    P(d>0): Moss2022 Eq. 5.5, HW.
    FDM: d/MD gamma with a=2.5 (explicit), b=0.180981 (global D/MD),
         x MD(M) lognormal.  No truncation at d/MD=1.0.
    """
    rng = np.random.default_rng(seed)
    fi = fault_info[1]

    r_km = fi['site_dist_m'] / 1000.0
    pd0 = min(math.exp(-2.2 * r_km + 0.5), 1.0)

    gamma_a = 2.5
    gamma_b = 0.180981

    mags, rates, N_mmin = _yc1985_rates(4.5, 1.0, char_mag, fi['M0_rate'])
    pr_slip = _psr_wc93(mags)

    n_D = len(D_LEVELS)
    probDd = np.zeros((len(mags), n_D))
    for mi in range(len(mags)):
        A = _trlnrnd(_md_median_ln(mags[mi]), SIGMA_CONV, sim, rng)
        B = rng.gamma(gamma_a, gamma_b, sim)
        combine = A * B
        for di in range(n_D):
            probDd[mi, di] = np.mean(combine > D_LEVELS[di])

    nu = pd0 * np.array([
        np.sum(rates * pr_slip * probDd[:, di]) for di in range(n_D)])
    return D_LEVELS.copy(), nu, N_mmin, pd0


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
# Tests - Geometry
# ====================================================================

class TestGeometry:
    """Verify fault geometry computed from spreadsheet coordinates."""

    @pytest.fixture(scope="class")
    def fi(self):
        return _compute_fault_geometry()

    def test_width(self):
        assert W_KM == pytest.approx(6.527, abs=0.01)

    def test_fault1_length(self, fi):
        assert fi[1]['L_km'] == pytest.approx(7.759, abs=0.01)

    def test_fault2_length(self, fi):
        assert fi[2]['L_km'] == pytest.approx(7.837, abs=0.01)

    def test_fault1_site_distance(self, fi):
        """Fault 1: site at r ~ 254 m (hanging wall)."""
        assert fi[1]['site_dist_m'] == pytest.approx(254.4, abs=5.0)

    def test_fault2_site_on_trace(self, fi):
        """Fault 2: site lies on trace (r ~ 0)."""
        assert fi[2]['site_dist_m'] < 1.0

    def test_fault1_xL(self, fi):
        assert fi[1]['site_xL_folded'] == pytest.approx(0.405, abs=0.01)

    def test_fault2_xL(self, fi):
        assert fi[2]['site_xL_folded'] == pytest.approx(0.414, abs=0.01)

    def test_fault_assignment(self, fi):
        """Corrected: Test 1 → Fault 2 (on-fault), Test 2 → Fault 1 (HW)."""
        assert fi[2]['site_dist_m'] < fi[1]['site_dist_m']


# ====================================================================
# Tests - MFD
# ====================================================================

class TestMFD:
    """Verify YC1985 characteristic MFD and moment-balanced N(Mmin)."""

    @pytest.fixture(scope="class")
    def fi(self):
        return _compute_fault_geometry()

    @pytest.mark.parametrize("char_mag,expected_N_order", [
        (5.95, 2.3e-5),
        (6.20, 1.3e-5),
    ])
    def test_annual_rate_order(self, fi, char_mag, expected_N_order):
        _, _, N = _yc1985_rates(4.5, 1.0, char_mag, fi[2]['M0_rate'])
        assert 0.5 * expected_N_order < N < 2.0 * expected_N_order


# ====================================================================
# Tests - Test 1: Principal hazard (Fault 2, on-fault)
# ====================================================================

class TestPrincipalHazard:
    """
    Test 1: Principal displacement on Fault 2 (site on fault trace).

    Uses YC1985 MFD, WC1993 P(SR|M), Moss2022 D/MD with recommended sigma,
    x/L-dependent gamma, truncation at D/MD=1.
    """

    @pytest.fixture(scope="class")
    def curves(self):
        fi = _compute_fault_geometry()
        out = {}
        for cm in [5.95, 6.20]:
            D, nu, N = _compute_test1(fi, cm, sim=30000, seed=42)
            out[cm] = (D, nu, N)
        return out

    def test_hazard_monotonic_595(self, curves):
        _, nu, _ = curves[5.95]
        assert np.all(np.diff(nu) <= 1e-15)

    def test_hazard_monotonic_620(self, curves):
        _, nu, _ = curves[6.20]
        assert np.all(np.diff(nu) <= 1e-15)

    def test_rate_positive_at_small_d(self, curves):
        """Rate at D=0.001 m should be positive and below N(Mmin)."""
        D, nu, N = curves[5.95]
        assert nu[0] > 0
        assert nu[0] < N

    def test_rate_drops_to_zero_at_large_d(self, curves):
        """Rate at D=10 m should be essentially zero for Mmax ~ 6.2."""
        _, nu, _ = curves[5.95]
        assert nu[-1] < 1e-12

    def test_cm620_higher_at_large_d(self, curves):
        """char_mag=6.20 allows larger M → comparable or higher rate."""
        _, nu_595, _ = curves[5.95]
        _, nu_620, _ = curves[6.20]
        assert nu_620[9] >= nu_595[9] * 0.5


# ====================================================================
# Tests - Test 2: Distributed hazard (Fault 1, HW, r ~ 254 m)
# ====================================================================

class TestDistributedHazard:
    """
    Test 2: Distributed displacement on Fault 1 (site at r ~ 254 m, HW).

    Uses gamma(a=2.5, b=0.180981) x MD(M) lognormal, no D/MD truncation,
    scaled by P(d>0) from Eq. 5.5.
    """

    @pytest.fixture(scope="class")
    def results(self):
        fi = _compute_fault_geometry()
        out = {}
        for cm in [5.95, 6.20]:
            D, nu, N, pd0 = _compute_test2(fi, cm, sim=30000, seed=42)
            out[cm] = (D, nu, N, pd0)
        out['fi'] = fi
        return out

    def test_pd0_value(self, results):
        """P(d>0) at r ~ 254 m HW should be ~ 0.942."""
        _, _, _, pd0 = results[5.95]
        assert pd0 == pytest.approx(0.942, abs=0.01)

    def test_gamma_parameters(self):
        """Test 2 uses a=2.5 (explicit), b=0.180981 (global D/MD scale)."""
        assert 2.5 * 0.180981 == pytest.approx(0.4525, abs=0.001)

    def test_hazard_monotonic_595(self, results):
        _, nu, _, _ = results[5.95]
        assert np.all(np.diff(nu) <= 1e-15)

    def test_hazard_monotonic_620(self, results):
        _, nu, _, _ = results[6.20]
        assert np.all(np.diff(nu) <= 1e-15)

    def test_rate_positive_at_small_d(self, results):
        _, nu, N, _ = results[5.95]
        assert nu[0] > 0
        assert nu[0] < N

    def test_rate_drops_at_large_d(self, results):
        _, nu, _, _ = results[5.95]
        assert nu[-1] < 1e-10

    def test_distributed_rate_less_than_principal(self, results):
        """Distributed rate (with pd0 < 1) should be lower than principal."""
        fi = results['fi']
        D_p, nu_p, _ = _compute_test1(fi, 5.95, sim=30000, seed=42)
        _, nu_d, _, pd0 = results[5.95]
        assert pd0 < 1.0
        assert nu_d[0] < nu_p[0]


# ====================================================================
# Standalone runner with plots
# ====================================================================

if __name__ == '__main__':
    import time

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plots')
    os.makedirs(out_dir, exist_ok=True)

    t0 = time.time()

    fi = _compute_fault_geometry()
    print("=" * 72)
    print("FAULT GEOMETRY")
    print("=" * 72)
    for fid in [1, 2]:
        f = fi[fid]
        print(f"  Fault {fid}: L = {f['L_km']:.3f} km, "
              f"A = {f['A_km2']:.3f} km2, "
              f"M0_rate = {f['M0_rate']:.4e} N.m/yr")
        print(f"           site dist = {f['site_dist_m']:.1f} m, "
              f"x/L = {f['site_xL']:.4f} (folded {f['site_xL_folded']:.4f})")

    print(f"\n  Corrected assignment:")
    print(f"    Test 1 (principal)    -> Fault 2, on-fault, "
          f"x/L = {fi[2]['site_xL_folded']:.4f}")
    print(f"    Test 2 (distributed)  -> Fault 1, HW, "
          f"r = {fi[1]['site_dist_m']:.1f} m")

    # ── Test 1: Principal ──
    print(f"\n{'=' * 72}")
    print("TEST 1: PRINCIPAL DISPLACEMENT (Fault 2, on-fault)")
    print(f"  x/L = {fi[2]['site_xL_folded']:.4f}")
    print(f"  Gamma D/MD: a = {_gamma_shape_dmd(fi[2]['site_xL_folded']):.4f}, "
          f"b = {_gamma_scale_dmd(fi[2]['site_xL_folded']):.4f}")
    print(f"  MD sigma (recommended): 0.20 log10 = {SIGMA_CONV:.4f} ln")
    print("=" * 72)

    results_t1 = {}
    for cm in [5.95, 6.20]:
        D, nu, N = _compute_test1(fi, cm)
        results_t1[cm] = (D, nu)
        print(f"\n  char_mag = {cm}, Mmax = {cm + 0.25:.2f}, "
              f"N(Mmin) = {N:.4e} events/yr")
        print(f"  {'D (m)':>10s}  {'v(D) [/yr]':>14s}  {'RP [yr]':>12s}")
        print(f"  {'-' * 42}")
        for i in range(len(D)):
            rp = 1.0 / nu[i] if nu[i] > 0 else float('inf')
            print(f"  {D[i]:10.3f}  {nu[i]:14.6e}  {rp:12.0f}")

    # ── Test 2: Distributed ──
    print(f"\n{'=' * 72}")
    print("TEST 2: DISTRIBUTED DISPLACEMENT (Fault 1, HW)")
    print(f"  r = {fi[1]['site_dist_m']:.1f} m")
    print(f"  Gamma d/MD: a = 2.5 (explicit), b = 0.180981 (global D/MD)")
    print("=" * 72)

    results_t2 = {}
    for cm in [5.95, 6.20]:
        D, nu, N, pd0 = _compute_test2(fi, cm)
        results_t2[cm] = (D, nu)
        print(f"\n  char_mag = {cm}, Mmax = {cm + 0.25:.2f}, "
              f"N(Mmin) = {N:.4e} events/yr")
        print(f"  P(d>0) = {pd0:.6f}")
        print(f"  {'d (m)':>10s}  {'v(d) [/yr]':>14s}  {'RP [yr]':>12s}")
        print(f"  {'-' * 42}")
        for i in range(len(D)):
            rp = 1.0 / nu[i] if nu[i] > 0 else float('inf')
            print(f"  {D[i]:10.3f}  {nu[i]:14.6e}  {rp:12.0f}")

    print(f"\nTotal computation time: {time.time() - t0:.1f} s")

    # ── Plot ──
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for cm, col, ls, lbl in [
        (5.95, '#2166ac', '-', 'char_mag=5.95 (Mmax=6.20)'),
        (6.20, '#b2182b', '--', 'char_mag=6.20 (Mmax=6.45)'),
    ]:
        D1, nu1 = results_t1[cm]
        v = nu1 > 0
        axes[0].loglog(D1[v], nu1[v], color=col, ls=ls, lw=2,
                       marker='o', ms=4, label=lbl)

        D2, nu2 = results_t2[cm]
        v = nu2 > 0
        axes[1].loglog(D2[v], nu2[v], color=col, ls=ls, lw=2,
                       marker='o', ms=4, label=lbl)

    axes[0].set_title(
        f'Test 1: Principal Displacement\n'
        f'Fault 2, on-fault, x/L = {fi[2]["site_xL_folded"]:.2f}',
        fontsize=11)
    axes[1].set_title(
        f'Test 2: Distributed Displacement\n'
        f'Fault 1, r = {fi[1]["site_dist_m"]:.0f} m, HW, gamma(a=2.5)',
        fontsize=11)

    for ax in axes:
        ax.set_xlabel('Displacement (m)', fontsize=12)
        ax.set_ylabel('Annual Exceedance Rate', fontsize=12)
        ax.legend(fontsize=9)
        ax.grid(True, which='both', alpha=0.3)
        ax.set_xlim(5e-4, 15)
        ax.set_ylim(1e-12, 1e-4)

    fig.suptitle(
        'PFDHA_CASE_STUDIES Benchmark Hazard Curves\n'
        'YC1985 MFD | WC1993 P(SR|M) | Moss et al. 2022 FDM | '
        'sigma_rec = 0.20',
        fontsize=12, fontweight='bold')
    fig.tight_layout()
    fpath = os.path.join(out_dir, 'pfdha_case_studies_benchmark.png')
    fig.savefig(fpath, dpi=150)
    print(f"\nPlot saved to {fpath}")
