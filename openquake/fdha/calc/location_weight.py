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

"""
Rupture-location weight ``W_p(r)`` for the principal fault-displacement
contribution.

See ``docs/design/rupture_location_uncertainty.md`` (section 2) for the
normative specification. ``W_p`` weights the principal contribution by the
probability that the mapped fault trace is the true rupture location at
across-strike distance ``r``. It has TWO SEPARATE, MUTUALLY EXCLUSIVE paths --
``h`` (``r_threshold_km``) enters ONLY the first, ``sigma`` enters ONLY the
second:

* ``sigma == 0`` -- the legacy boxcar ``W_p(r) = 1{|r| <= h}``. The principal
  zone is a hard half-width ``h``; the distributed contribution is combined
  complementarily (``G = 1 - W_p``), the historical principal/distributed
  split. ``h`` is a structure/rupture-zone property, NOT a mapping error.

* ``sigma > 0`` -- Petersen's Gaussian mapping-error weight, with ``h`` playing
  NO role. The true rupture offset from the mapped trace is taken to be normal,
  ``N(0, sigma)`` (Petersen et al. 2011, Tables 2-3 two-sided mapping errors,
  p. 810), pinned to 1 on the trace ("probability of fault rupture on the fault
  to be one", p. 819) and truncated beyond +/-n-sigma ("within 2 standard
  deviations of the mapped fault trace", p. 819, n = 2 by default)::

      W_p(r) = exp(-r^2 / (2 sigma^2))   for |r| <= n*sigma,   else 0

  This is the standard-normal density normalised to its peak (pin at r = 0),
  i.e. a pure bell -- no boxcar, no plateau. Its width scales with sigma and it
  vanishes at the toe ``|r| = n*sigma``.

The two paths are deliberately NOT a smooth limit of one another: ``h`` and
``sigma`` describe different things (a fixed rupture-zone half-width vs a
mapping uncertainty), so choosing ``sigma > 0`` replaces the boxcar with the
Gaussian rather than blurring the boxcar's edges. The site footprint ``z``
enters neither path -- it belongs to the distributed rupture probability
(Petersen electronic supplement / Tables 4-5) and the near-field floor.
"""

import numpy as np

__all__ = ["location_weight"]


def location_weight(
    r,
    r_threshold_km: float,
    r_sigma_km: float,
    r_sigma_truncation: float = 2.0,
) -> np.ndarray:
    """Principal-contribution location weight ``W_p(r)``.

    All distances are in kilometres. The result is a ``float64`` array with the
    shape of ``r``.

    :param r:
        Across-strike distance(s) from the mapped fault trace (``ctx.r``), km.
        Only ``|r|`` matters; the weight is symmetric about the trace.
    :param r_threshold_km:
        Boxcar half-width ``h``. Used ONLY on the ``sigma == 0`` path; ignored
        entirely when ``sigma > 0`` (the Gaussian path has no boxcar).
    :param r_sigma_km:
        Two-sided mapping-accuracy sigma (Petersen Tables 2-3). ``0`` selects
        the legacy boxcar of half-width ``h``; ``> 0`` selects Petersen's
        pinned, +/-n-sigma-truncated Gaussian mapping-error weight (``h`` plays
        no role).
    :param r_sigma_truncation:
        Truncation ``n`` in +/-n-sigma (Petersen p. 819: the principal fault
        occurs "within 2 standard deviations" -- n = 2). Gaussian path only.
    :returns:
        ``W_p`` as a ``float64`` array, ``W_p(0) == 1``. On the boxcar path it
        is ``1{|r| <= h}``; on the Gaussian path it is a pinned bell,
        monotone non-increasing in ``|r|`` and exactly ``0`` for
        ``|r| > n*sigma``.
    """
    r = np.abs(np.asarray(r, dtype=np.float64))
    sigma = float(r_sigma_km)

    if sigma == 0.0:
        # Path 1 -- boxcar. The principal zone is |r| <= h (h = r_threshold_km).
        # Combined complementarily with distributed (G = 1 - W_p): the
        # historical principal/distributed split, reproduced bit-for-bit
        # (docs section 2, verification V1). h is used here and ONLY here.
        h = float(r_threshold_km)
        return (r <= h).astype(np.float64)

    # Path 2 -- Petersen's Gaussian mapping-error weight. The rupture offset
    # from the mapped trace is normal N(0, sigma); W_p is that density pinned to
    # its peak (== 1 on the trace, p. 819) and truncated beyond +/-n-sigma
    # ("within 2 standard deviations", p. 819). h does NOT enter -- this is a
    # pure bell, not a smoothed boxcar (docs section 2).
    n = float(r_sigma_truncation)
    wp = np.exp(-(r * r) / (2.0 * sigma * sigma))
    # np.where (not boolean index-assignment) so a 0-d/scalar r survives with
    # its shape, matching the sigma==0 path and the old kernel's contract.
    return np.where(r > n * sigma, 0.0, wp)
