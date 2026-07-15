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

The properties exercised are, per the design doc section 6:

* exact boxcar recovery at ``sigma == 0``;
* ``W_p(0) == 1`` pinned (``sigma > 0``);
* monotone non-increasing in ``|r|``;
* exactly ``0`` beyond ``z/2 + n*sigma``;
* ``sigma -> 0+`` collapses to the ``z/2`` boxcar -- which differs from the
  ``h`` boxcar unless ``h == z/2`` (a documented discontinuity between the two
  limits).
"""

import numpy as np
import pytest

from openquake.fdha.calc.location_weight import location_weight


# Petersen (2011) two-sided mapping accuracy classes, km (Tables 2-3).
PETERSEN_SIGMAS = [0.0269, 0.0438, 0.0655, 0.0727, 0.116]


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


@pytest.mark.parametrize("sigma", PETERSEN_SIGMAS)
def test_pinned_at_zero(sigma):
    """W_p(0) == 1 exactly (pinning subsumes truncated-normal renorm)."""
    wp0 = location_weight(0.0, r_threshold_km=0.1, r_sigma_km=sigma,
                          site_footprint_m=25.0, r_sigma_truncation=2.0)
    assert float(wp0) == 1.0


@pytest.mark.parametrize("sigma", PETERSEN_SIGMAS)
def test_monotone_non_increasing_in_abs_r(sigma):
    """W_p is monotone non-increasing as |r| grows."""
    r = np.linspace(0.0, 1.0, 2001)
    wp = location_weight(r, r_threshold_km=0.1, r_sigma_km=sigma,
                         site_footprint_m=25.0, r_sigma_truncation=2.0)
    diffs = np.diff(wp)
    assert np.all(diffs <= 1e-15)


@pytest.mark.parametrize("sigma", PETERSEN_SIGMAS)
def test_exactly_zero_beyond_support(sigma):
    """W_p == 0 for |r| > z/2 + n*sigma, and > 0 just inside it."""
    z_m, n = 25.0, 2.0
    half = (z_m / 1000.0) / 2.0
    edge = half + n * sigma
    inside = location_weight(edge - 1e-6, r_threshold_km=0.1, r_sigma_km=sigma,
                             site_footprint_m=z_m, r_sigma_truncation=n)
    beyond = location_weight(
        np.array([edge + 1e-9, edge + 0.01, 1.0]),
        r_threshold_km=0.1, r_sigma_km=sigma, site_footprint_m=z_m,
        r_sigma_truncation=n)
    assert float(inside) > 0.0
    assert np.all(beyond == 0.0)


def test_symmetric_about_trace():
    """W_p(r) == W_p(-r): only |r| matters."""
    r = np.linspace(0.0, 0.4, 501)
    wp_pos = location_weight(r, r_threshold_km=0.1, r_sigma_km=0.0655,
                             site_footprint_m=25.0)
    wp_neg = location_weight(-r, r_threshold_km=0.1, r_sigma_km=0.0655,
                             site_footprint_m=25.0)
    np.testing.assert_array_equal(wp_pos, wp_neg)


def test_bounded_between_zero_and_one():
    """0 <= W_p <= 1 everywhere (pinned mass can never exceed M(0))."""
    r = np.linspace(-0.6, 0.6, 4001)
    for sigma in PETERSEN_SIGMAS:
        wp = location_weight(r, r_threshold_km=0.1, r_sigma_km=sigma,
                             site_footprint_m=25.0)
        assert wp.min() >= 0.0
        assert wp.max() <= 1.0


def test_sigma_to_zero_collapses_to_half_z_boxcar():
    """sigma -> 0+ collapses to the z/2 boxcar, NOT the h boxcar.

    This is the documented discontinuity between the two limits: the sigma == 0
    branch uses the h = r_threshold_km boxcar, but the sigma > 0 branch shrinks
    onto the footprint half-width z/2 as sigma vanishes. They coincide only when
    h == z/2. Here z/2 = 12.5 m = 0.0125 km, deliberately != h = 0.1 km.
    """
    z_m = 25.0
    half = (z_m / 1000.0) / 2.0  # 0.0125 km
    h = 0.1                       # != half
    tiny = 1e-7
    # Well inside z/2 -> ~1; well outside z/2 but inside h -> ~0.
    inside_half = location_weight(0.5 * half, r_threshold_km=h, r_sigma_km=tiny,
                                  site_footprint_m=z_m)
    between = location_weight(0.5 * (half + h), r_threshold_km=h,
                              r_sigma_km=tiny, site_footprint_m=z_m)
    assert float(inside_half) == pytest.approx(1.0, abs=1e-6)
    assert float(between) == pytest.approx(0.0, abs=1e-6)
    # The h boxcar would have returned 1.0 at `between` (|r| < h); it does not.
    assert float(between) < 0.5


def test_shape_and_dtype_preserved():
    """Array shape is preserved; scalar in -> 0-d float64 out."""
    r2d = np.zeros((3, 4))
    wp2d = location_weight(r2d, r_threshold_km=0.1, r_sigma_km=0.05,
                           site_footprint_m=25.0)
    assert wp2d.shape == (3, 4)
    assert wp2d.dtype == np.float64
    scalar = location_weight(0.0, r_threshold_km=0.1, r_sigma_km=0.05)
    assert np.ndim(scalar) == 0
