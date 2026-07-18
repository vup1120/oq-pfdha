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

"""Unit tests for the rupture-location weight ``W_p(r)`` (verification V2 of
``docs/design/rupture_location_uncertainty.md``).

``W_p`` has two separate, mutually exclusive paths (design doc section 2):

* ``sigma == 0`` -- the legacy boxcar ``1{|r| <= h}`` (``h`` = r_threshold_km),
  combined complementarily with distributed. ``h`` is used only here.
* ``sigma > 0`` -- Petersen's pinned, +/-n-sigma-truncated Gaussian
  mapping-error weight ``exp(-r^2/2 sigma^2)``. ``h`` plays NO role; the weight
  is a pure bell with toe at ``|r| = n*sigma``.

The two are deliberately NOT a smooth limit of one another -- ``h`` (a fixed
rupture-zone half-width) and ``sigma`` (a mapping uncertainty) describe
different things.
"""

import numpy as np
import pytest

from openquake.fdha.calc.location_weight import location_weight


# Petersen (2011) two-sided mapping accuracy classes, km (Tables 2-3).
PETERSEN_SIGMAS = [0.0269, 0.0438, 0.0655, 0.0727, 0.116]


# --------------------------------------------------------------------------
# Path 1: boxcar (sigma == 0). h is used here and only here.
# --------------------------------------------------------------------------
def test_sigma0_recovers_boxcar_exactly():
    """sigma == 0 must reproduce the legacy boxcar |r| <= h bit-for-bit."""
    h = 0.1
    r = np.linspace(-0.5, 0.5, 4001)
    wp = location_weight(r, r_threshold_km=h, r_sigma_km=0.0)
    boxcar = (np.abs(r) <= h).astype(np.float64)
    # Byte-identical: the kernel's V1 byte-identity gate depends on this.
    assert wp.dtype == np.float64
    assert wp.tobytes() == boxcar.tobytes()


def test_sigma0_boxcar_edge_inclusive():
    """The boxcar is closed at |r| == h (matches np.abs(r) <= h)."""
    h = 0.1
    wp = location_weight(np.array([h, -h, h + 1e-9]), r_threshold_km=h,
                         r_sigma_km=0.0)
    assert wp[0] == 1.0 and wp[1] == 1.0 and wp[2] == 0.0


# --------------------------------------------------------------------------
# Path 2: Gaussian (sigma > 0). h must NOT influence the result.
# --------------------------------------------------------------------------
@pytest.mark.parametrize("sigma", PETERSEN_SIGMAS)
def test_pinned_at_zero(sigma):
    """W_p(0) == 1 exactly (Gaussian density pinned to its peak)."""
    wp0 = location_weight(0.0, r_threshold_km=0.1, r_sigma_km=sigma,
                          r_sigma_truncation=2.0)
    assert float(wp0) == 1.0


@pytest.mark.parametrize("sigma", PETERSEN_SIGMAS)
def test_gaussian_shape(sigma):
    """On the sigma>0 path W_p is exactly exp(-r^2/2 sigma^2) inside the
    truncation, independent of h."""
    r = np.linspace(0.0, 1.9 * sigma, 40)
    wp = location_weight(r, r_threshold_km=0.1, r_sigma_km=sigma,
                         r_sigma_truncation=2.0)
    expected = np.exp(-(r ** 2) / (2.0 * sigma ** 2))
    np.testing.assert_allclose(wp, expected, rtol=0, atol=1e-12)


@pytest.mark.parametrize("h", [0.0, 0.02, 0.1, 5.0])
def test_gaussian_independent_of_h(h):
    """Changing h must not change the sigma>0 weight at all."""
    sigma = 0.0655
    r = np.linspace(-0.3, 0.3, 601)
    ref = location_weight(r, r_threshold_km=0.1, r_sigma_km=sigma)
    got = location_weight(r, r_threshold_km=h, r_sigma_km=sigma)
    np.testing.assert_array_equal(got, ref)


@pytest.mark.parametrize("sigma", PETERSEN_SIGMAS)
def test_monotone_non_increasing_in_abs_r(sigma):
    """W_p is monotone non-increasing as |r| grows."""
    r = np.linspace(0.0, 1.0, 2001)
    wp = location_weight(r, r_threshold_km=0.1, r_sigma_km=sigma,
                         r_sigma_truncation=2.0)
    diffs = np.diff(wp)
    assert np.all(diffs <= 1e-15)


@pytest.mark.parametrize("sigma", PETERSEN_SIGMAS)
def test_exactly_zero_beyond_support(sigma):
    """W_p == 0 for |r| > n*sigma, and > 0 just inside it.

    The toe is n*sigma (Petersen p. 819, "within 2 standard deviations"),
    independent of h -- the Gaussian path has no boxcar term.
    """
    n = 2.0
    edge = n * sigma
    inside = location_weight(edge - 1e-6, r_threshold_km=0.1, r_sigma_km=sigma,
                             r_sigma_truncation=n)
    beyond = location_weight(
        np.array([edge + 1e-9, edge + 0.01, 2.0]),
        r_threshold_km=0.1, r_sigma_km=sigma, r_sigma_truncation=n)
    assert float(inside) > 0.0
    assert np.all(beyond == 0.0)


def test_truncation_at_n_sigma_value():
    """Just inside the toe W_p equals exp(-n^2/2); just outside it is 0."""
    sigma, n = 0.05, 2.0
    inside = float(location_weight(n * sigma - 1e-9, r_threshold_km=0.1,
                                   r_sigma_km=sigma, r_sigma_truncation=n))
    assert inside == pytest.approx(np.exp(-n * n / 2.0), rel=1e-6)


def test_symmetric_about_trace():
    """W_p(r) == W_p(-r): only |r| matters."""
    r = np.linspace(0.0, 0.4, 501)
    wp_pos = location_weight(r, r_threshold_km=0.1, r_sigma_km=0.0655)
    wp_neg = location_weight(-r, r_threshold_km=0.1, r_sigma_km=0.0655)
    np.testing.assert_array_equal(wp_pos, wp_neg)


def test_bounded_between_zero_and_one():
    """0 <= W_p <= 1 everywhere (a peak-normalised density)."""
    r = np.linspace(-0.6, 0.6, 4001)
    for sigma in PETERSEN_SIGMAS:
        wp = location_weight(r, r_threshold_km=0.1, r_sigma_km=sigma)
        assert wp.min() >= 0.0
        assert wp.max() <= 1.0


def test_sigma_to_zero_is_a_spike_not_the_h_boxcar():
    """sigma -> 0+ collapses onto a spike at r = 0, NOT the h boxcar.

    The two paths describe different things and are deliberately not a smooth
    limit of one another: as sigma vanishes the Gaussian narrows to a point at
    the trace. In particular, points inside |r| < h that the boxcar would keep
    at 1 go to 0 here (unless they are also within ~sigma of the trace).
    """
    h = 0.1
    tiny = 1e-7
    on_trace = location_weight(0.0, r_threshold_km=h, r_sigma_km=tiny)
    inside_boxcar = location_weight(np.array([0.5 * h, h - 1e-4]),
                                    r_threshold_km=h, r_sigma_km=tiny)
    assert float(on_trace) == 1.0
    np.testing.assert_allclose(inside_boxcar, 0.0, atol=1e-9)


def test_shape_and_dtype_preserved():
    """Array shape is preserved; scalar in -> 0-d float64 out."""
    r2d = np.zeros((3, 4))
    wp2d = location_weight(r2d, r_threshold_km=0.1, r_sigma_km=0.05)
    assert wp2d.shape == (3, 4)
    assert wp2d.dtype == np.float64
    scalar = location_weight(0.0, r_threshold_km=0.1, r_sigma_km=0.05)
    assert np.ndim(scalar) == 0
