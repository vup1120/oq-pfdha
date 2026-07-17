# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#
# Copyright (C) 2024-2026 Yen-Shin Chen, OGS
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""Benchmark: Petersen et al. (2011) Fig. 9c (paper p. 820) vs the tool's
rupture-location kernel (roadmap stage C5).

Setup (paper's worked example, p. 819): characteristic M 7.0 every 140 yr
(alpha = 1/140/yr), P(sr != 0 | m) from eq. 5 (``WC1993PrimarySR``), bilinear
principal displacement at the mid-rupture branch l/L = 0.5 >= 0.3
(``Petersen2011PrimaryFD_bilinear``, eq. 8), distributed occurrence + eq. 18
displacement at 25-m cells (``Petersen2011SecondarySR``,
``Petersen2011SecondaryFD``, near-field floored per D7 at
``NEAR_FIELD_FLOOR_KM``), rupture-location weight ``W_p(r)`` from
``openquake.fdha.calc.location_weight.location_weight`` on the sigma > 0
(pure Gaussian, pinned, +/-2 sigma truncated) path, additively combined with
the distributed term per the kernel spec (D1, sigma > 0 -> G(r) = 1):

    lambda_total(D0, r) = alpha * P_sr * [W_p(r) * P(D_p > D0)
                                          + P(d != 0 | r) * P(D_d > D0 | r)]

D0(r) is solved by bisection for lambda_total(D0, r) == lambda* = -ln(0.9)/50
(10% in 50 yr). This is exactly ``v4_fig9c.py`` / ``peak_audit.py``'s
"pin+comp" analytic construction, but built from the tool's own model classes
and kernel (``location_weight``, ``WC1993PrimarySR``,
``Petersen2011PrimaryFD_bilinear``, ``Petersen2011SecondarySR``,
``Petersen2011SecondaryFD``) instead of re-implementing the formulas locally.

See README.md in this directory for the digitization provenance and the
documented text-vs-figure inconsistency (0.90 on-trace pin / 1.65x sigma
width, neither derivable from the paper -- ``docs/design/
rupture_location_uncertainty.md`` section 2). Policy: the tool implements
the paper's *stated* method; this benchmark compares that stated method
against the digitized figure on an explicitly documented-offset basis, and
does NOT tune anything to fit the printed curves.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from openquake.fdha.calc.location_weight import location_weight
from openquake.fdha.calc.model_adapter import NEAR_FIELD_FLOOR_KM
from openquake.fdha.primary_surf_displ import Petersen2011PrimaryFD_bilinear
from openquake.fdha.primary_surf_rup import WC1993PrimarySR
from openquake.fdha.secondary_surf_displ import Petersen2011SecondaryFD
from openquake.fdha.secondary_surf_rup import Petersen2011SecondarySR

pytestmark = [pytest.mark.benchmark, pytest.mark.petersen2011]

HERE = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# Worked-example setup (paper p. 819)
# ---------------------------------------------------------------------------
ALPHA = 1.0 / 140.0                  # characteristic M7 rate, /yr
MAG = 7.0
X_L = 0.5                            # mid-rupture, l/L >= 0.3 (bilinear branch)
LAMBDA_TARGET = -np.log(0.9) / 50.0  # 10% in 50 yr, Poisson
R_SIGMA_TRUNCATION = 2.0             # "within 2 standard deviations", p. 819
PIXEL_SIZE = 25                      # 25-m cells (Fig. 9c)

# Petersen (2011) Tables 2-3 two-sided mapping-accuracy sigma, km. Matches
# openquake/fdha/test/unit/test_location_weight.py's PETERSEN_SIGMAS and the
# four classes plotted in Fig. 9c (Table 3's fifth value, 0.116, covers the
# complex-inferred/complex-concealed classes and is not one of the four
# plotted curves -- see README.md and 09-Uncertainty.md for the full table).
SIGMA_KM = {
    "Accurate": 0.02689,
    "Approximate": 0.04382,
    "Concealed": 0.06552,
    "Inferred": 0.07269,
}

# Reverse-engineered figure offset (docs/design/rupture_location_uncertainty.md
# section 2, "Validation vs Fig. 9c -- settled by digitization"): the printed
# curves collapse onto ONE universal Gaussian in r/sigma with an effective
# width ~1.65x the table sigma and an on-trace pin ~0.90, not derivable from
# anything stated in the paper.
FIGURE_OFFSET_PIN = 0.90
FIGURE_OFFSET_KSIG = 1.65

# ---------------------------------------------------------------------------
# Tool model instances (real kernel + real Petersen model pieces)
# ---------------------------------------------------------------------------
_PSR_MODEL = WC1993PrimarySR()
_PFD_MODEL = Petersen2011PrimaryFD_bilinear()
_SSR_MODEL = Petersen2011SecondarySR()
_SFD_MODEL = Petersen2011SecondaryFD()

P_SR = float(_PSR_MODEL.get_prob(mag=MAG))  # eq. 5


def _p_primary_exceed(d0_cm: np.ndarray) -> np.ndarray:
    """P(D_principal > d0) at X/L = 0.5, mag 7 (eq. 8, mid-rupture branch)."""
    d0_m = np.atleast_1d(d0_cm).astype(float) / 100.0
    prob = _PFD_MODEL.get_prob(d=d0_m, X_L_ratio=np.array([X_L]), mag=MAG)
    return np.asarray(prob)[:, 0]


def _p_dist_occurrence(r_km: float) -> float:
    """P(d != 0 | r) -- distributed occurrence, near-field Table 5 variant."""
    return float(_SSR_MODEL.get_prob(r=np.array([abs(r_km)]),
                                      pixel_size=PIXEL_SIZE,
                                      version="near_field")[0])


def _p_dist_exceed(d0_cm: np.ndarray, r_km: float) -> np.ndarray:
    """P(D_distributed > d0 | r) -- eq. 18, near-field floored per D7.

    ``Petersen2011SecondaryFD.get_prob`` collapses to a bare ``float`` when
    both ``d`` and ``r`` are length-1 (matches the Youngs2003-style calling
    convention used throughout the model package) -- ``atleast_1d`` restores
    a consistent array shape for callers here.
    """
    r_floored_km = max(abs(r_km), NEAR_FIELD_FLOOR_KM)
    d0_m = np.atleast_1d(d0_cm).astype(float) / 100.0
    prob = _SFD_MODEL.get_prob(d=d0_m, mag=MAG, r=np.array([r_floored_km]))
    return np.atleast_1d(prob)


def lambda_total(d0_cm: float, r_km: float, sigma_km: float) -> float:
    """Annual rate of exceeding ``d0_cm`` at across-strike distance ``r_km``.

    Additive kernel (D1, sigma > 0 path): lambda_principal + lambda_distributed,
    W_p(r) from the tool's real ``location_weight`` kernel, G(r) = 1 (no
    complementary masking on the Gaussian path).
    """
    wp = float(location_weight(np.array([r_km]), r_threshold_km=0.1,
                                r_sigma_km=sigma_km,
                                r_sigma_truncation=R_SIGMA_TRUNCATION)[0])
    p_p = float(_p_primary_exceed(np.array([d0_cm]))[0])
    p_d_occ = _p_dist_occurrence(r_km)
    p_d_exc = float(_p_dist_exceed(np.array([d0_cm]), r_km)[0])
    return ALPHA * P_SR * (wp * p_p + p_d_occ * p_d_exc)


def solve_d0(r_km: float, sigma_km: float,
             lo_cm: float = 1e-3, hi_cm: float = 5000.0,
             iters: int = 70) -> float:
    """Displacement (cm) at exceedance rate ``LAMBDA_TARGET``; 0 if unreachable."""
    if lambda_total(lo_cm, r_km, sigma_km) < LAMBDA_TARGET:
        return 0.0
    for _ in range(iters):
        mid = np.sqrt(lo_cm * hi_cm)
        if lambda_total(mid, r_km, sigma_km) > LAMBDA_TARGET:
            lo_cm = mid
        else:
            hi_cm = mid
    return float(np.sqrt(lo_cm * hi_cm))


def stated_method_curve(r_grid_m: np.ndarray, sigma_km: float) -> np.ndarray:
    """D0(r) (cm) on the tool's stated method (real kernel + real models)."""
    r_grid_km = r_grid_m / 1000.0
    return np.array([solve_d0(r, sigma_km) for r in r_grid_km])


def figure_offset_curve(r_grid_m: np.ndarray, sigma_km: float) -> np.ndarray:
    """D0(r) (cm) using the documented figure-equivalent weight
    ``0.90 * exp(-r^2 / 2*(1.65*sigma)^2)`` in place of the stated W_p(r), to
    test the documented offset -- NOT the tool's default behaviour.
    """
    sigma_eff_km = FIGURE_OFFSET_KSIG * sigma_km
    r_grid_km = r_grid_m / 1000.0
    out = []
    for r_km in r_grid_km:
        wp = FIGURE_OFFSET_PIN * float(
            location_weight(np.array([r_km]), r_threshold_km=0.1,
                             r_sigma_km=sigma_eff_km,
                             r_sigma_truncation=R_SIGMA_TRUNCATION)[0])

        def lam(d0_cm, wp=wp, r_km=r_km):
            p_p = float(_p_primary_exceed(np.array([d0_cm]))[0])
            p_d_occ = _p_dist_occurrence(r_km)
            p_d_exc = float(_p_dist_exceed(np.array([d0_cm]), r_km)[0])
            return ALPHA * P_SR * (wp * p_p + p_d_occ * p_d_exc)

        if lam(1e-3) < LAMBDA_TARGET:
            out.append(0.0)
            continue
        lo, hi = 1e-3, 5000.0
        for _ in range(70):
            mid = np.sqrt(lo * hi)
            lo, hi = (mid, hi) if lam(mid) > LAMBDA_TARGET else (lo, mid)
        out.append(float(np.sqrt(lo * hi)))
    return np.array(out)


ANALYTIC_PEAK_CM = 128.9  # peak_audit.py, principal-only + eq.18 union bump

# ---------------------------------------------------------------------------
# Shared fixtures: build all four class curves once per test session
# ---------------------------------------------------------------------------
_R_GRID_M = np.arange(0.0, 301.0, 4.0)  # matches digitized r-range (+/-300 m)


@pytest.fixture(scope="module")
def digitized():
    data = np.load(HERE / "fig9c_digitized.npz")
    return {name: data[name] for name in SIGMA_KM}


@pytest.fixture(scope="module")
def stated_curves():
    return {name: stated_method_curve(_R_GRID_M, sigma)
            for name, sigma in SIGMA_KM.items()}


@pytest.fixture(scope="module")
def transformed_curves():
    return {name: figure_offset_curve(_R_GRID_M, sigma)
            for name, sigma in SIGMA_KM.items()}


def _rms_against_digitized(r_grid_m, curve_cm, r_dig, d_dig, min_d=2.0):
    """RMS (cm) of ``curve_cm`` (defined on ``r_grid_m``, |r| >= 0) against
    the digitized (r, d) samples, interpolated onto |r_dig|, restricted to
    d_dig > min_d cm (foot noise / zero-baseline exclusion, matches
    roundtrip.py)."""
    sel = d_dig > min_d
    r_abs = np.abs(r_dig[sel])
    d_sel = d_dig[sel]
    recon = np.interp(r_abs, r_grid_m, curve_cm, left=0.0, right=0.0)
    return float(np.sqrt(np.mean((recon - d_sel) ** 2)))


# ---------------------------------------------------------------------------
# Assertions
# ---------------------------------------------------------------------------
def test_pinned_equal_peaks(stated_curves):
    """The four class curves peak at the same value at r=0 (W_p pinned to 1
    regardless of sigma), and that peak matches the closed-form analytic
    peak (peak_audit.py) within +/-10%."""
    peaks = {name: curve[0] for name, curve in stated_curves.items()}
    values = np.array(list(peaks.values()))
    spread = (values.max() - values.min()) / values.mean()
    print(f"peaks (cm): {peaks}  relative spread={spread:.2e}  "
          f"analytic={ANALYTIC_PEAK_CM:.1f} cm")
    assert spread < 1e-9, (
        f"class peaks should be identical (W_p pinned to 1 at r=0, "
        f"independent of sigma); relative spread={spread:.2e}: {peaks}")

    mean_peak = float(values.mean())
    ratio = mean_peak / ANALYTIC_PEAK_CM
    assert 0.90 <= ratio <= 1.10, (
        f"stated-method peak {mean_peak:.1f} cm vs analytic "
        f"{ANALYTIC_PEAK_CM:.1f} cm, ratio {ratio:.3f} outside +/-10%")


def test_tails_widen_and_reach_zero_beyond_2sigma(stated_curves):
    """Tails widen monotonically with sigma; each curve's W_p toe (and hence
    the principal contribution) is exactly zero beyond 2*sigma, so the
    *reach* of appreciable displacement grows with the class's sigma."""
    ordered = sorted(SIGMA_KM.items(), key=lambda kv: kv[1])
    reach_m = {}
    for name, sigma_km in ordered:
        curve = stated_curves[name]
        above = _R_GRID_M[curve > 1.0]  # 1 cm threshold, ignores float noise
        reach_m[name] = float(above.max()) if above.size else 0.0
    print(f"reach (m, D0>1cm) by class: {reach_m}")

    reaches = [reach_m[name] for name, _ in ordered]
    assert all(a <= b + 1e-9 for a, b in zip(reaches, reaches[1:])), (
        f"tail reach should be non-decreasing with sigma: {reach_m}")
    # strictly widening across the full accurate->inferred span
    assert reach_m["Inferred"] > reach_m["Accurate"], (
        f"Inferred tail ({reach_m['Inferred']} m) should reach farther than "
        f"Accurate ({reach_m['Accurate']} m)")

    # W_p toe: the kernel's documented contract (location_weight.py) is
    # "exactly 0 for |r| > n*sigma" -- a CLOSED interval, so r == 2*sigma is
    # still (barely) inside and only strictly-beyond vanishes.
    for name, sigma_km in SIGMA_KM.items():
        wp_beyond = float(location_weight(
            np.array([2.0 * sigma_km + 1e-6]), r_threshold_km=0.1,
            r_sigma_km=sigma_km, r_sigma_truncation=R_SIGMA_TRUNCATION)[0])
        assert wp_beyond == 0.0, (
            f"{name}: W_p should vanish strictly beyond 2*sigma "
            f"(got {wp_beyond})")


def test_documented_offset_fits_digitized_better(digitized, stated_curves,
                                                   transformed_curves):
    """The documented figure offset (sigma -> 1.65*sigma, pin 0.90) must fit
    the digitized curves substantially better than the raw stated method.
    Per-class: rms(transformed vs digitized) < rms(stated vs digitized), and
    the transformed rms stays below ~0.15 m (measured 7.5-12 cm vs 37-49 cm
    in the design-doc digitization; headroom, not overfit)."""
    summary = {}
    for name in SIGMA_KM:
        r_dig, d_dig = digitized[name]
        rms_stated = _rms_against_digitized(_R_GRID_M, stated_curves[name],
                                             r_dig, d_dig)
        rms_transformed = _rms_against_digitized(
            _R_GRID_M, transformed_curves[name], r_dig, d_dig)
        summary[name] = (rms_stated, rms_transformed)
        print(f"{name:12s} rms stated={rms_stated:6.1f} cm  "
              f"rms transformed={rms_transformed:6.1f} cm")

        assert rms_transformed < rms_stated, (
            f"{name}: documented-offset transform should fit the digitized "
            f"curve better than the raw stated method "
            f"(transformed={rms_transformed:.1f} cm, stated={rms_stated:.1f} "
            f"cm)")
        assert rms_transformed < 15.0, (
            f"{name}: transformed-vs-digitized rms {rms_transformed:.1f} cm "
            f"exceeds the 15 cm headroom bound")

    print("summary (class: rms_stated cm, rms_transformed cm):", summary)


def test_peak_ratio_digitized_over_stated(digitized, stated_curves):
    """Peak-of-digitized / peak-of-stated-method ~ 0.90 (the reverse-engineered
    on-trace pin), within a reasonable band."""
    stated_peak = float(np.mean([c[0] for c in stated_curves.values()]))
    digitized_peaks = []
    for name in SIGMA_KM:
        _r_dig, d_dig = digitized[name]
        digitized_peaks.append(float(d_dig.max()))
    digitized_peak = float(np.mean(digitized_peaks))
    ratio = digitized_peak / stated_peak
    print(f"digitized peak (mean of 4 classes) = {digitized_peak:.1f} cm, "
          f"stated peak = {stated_peak:.1f} cm, ratio = {ratio:.3f}")
    assert 0.85 <= ratio <= 0.95, (
        f"digitized/stated peak ratio {ratio:.3f} outside [0.85, 0.95] "
        f"(expected ~{FIGURE_OFFSET_PIN})")


if __name__ == "__main__":
    # Manual run: print the full summary table without pytest.
    stated = {name: stated_method_curve(_R_GRID_M, sigma)
              for name, sigma in SIGMA_KM.items()}
    transformed = {name: figure_offset_curve(_R_GRID_M, sigma)
                   for name, sigma in SIGMA_KM.items()}
    data = {name: np.load(HERE / "fig9c_digitized.npz")[name]
            for name in SIGMA_KM}
    print(f"P_sr(M{MAG}) = {P_SR:.4f}  lambda* = {LAMBDA_TARGET:.6f}/yr  "
          f"analytic peak = {ANALYTIC_PEAK_CM:.1f} cm")
    for name, sigma in SIGMA_KM.items():
        peak = stated[name][0]
        r_dig, d_dig = data[name]
        rms_s = _rms_against_digitized(_R_GRID_M, stated[name], r_dig, d_dig)
        rms_t = _rms_against_digitized(_R_GRID_M, transformed[name], r_dig,
                                        d_dig)
        print(f"{name:12s} sigma={sigma*1000:5.1f} m  peak={peak:6.1f} cm  "
              f"rms_stated={rms_s:6.1f} cm  rms_transformed={rms_t:6.1f} cm")
