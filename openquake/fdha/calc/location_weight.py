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
normative specification. In brief, ``W_p`` weights the principal contribution
by the probability that the mapped fault trace is the true rupture location at
across-strike distance ``r``:

* ``sigma == 0`` reproduces the legacy boxcar ``|r| <= h`` exactly (``h`` is
  ``r_threshold_km``) -- the historical principal/distributed split.
* ``sigma > 0`` is the pinned, +/-n-sigma-truncated normal *footprint mass*.
  This carries mapping accuracy: Petersen et al. (2011) give two-sided mapping
  errors per accuracy class (Tables 2-3, p. 810) and use +/-2-sigma normal
  weights in their worked example (p. 819). The mass is pinned so that
  ``W_p(0) == 1`` -- Petersen p. 819 takes "the probability of fault rupture on
  the fault to be one"; this pinning (M(0) normalisation) subsumes the explicit
  truncated-normal renormalisation.
"""

import numpy as np
from scipy.special import ndtr

__all__ = ["location_weight"]


def location_weight(
    r,
    r_threshold_km: float,
    r_sigma_km: float,
    site_footprint_m: float = 25.0,
    r_sigma_truncation: float = 2.0,
) -> np.ndarray:
    """Principal-contribution location weight ``W_p(r)``.

    All distances are in kilometres. The result is a ``float64`` array with the
    shape of ``r``.

    :param r:
        Across-strike distance(s) from the mapped fault trace (``ctx.r``), km.
        Only ``|r|`` matters; the weight is symmetric about the trace.
    :param r_threshold_km:
        Boxcar half-width ``h``. Used only on the ``sigma == 0`` path.
    :param r_sigma_km:
        Two-sided mapping-accuracy sigma (Petersen Tables 2-3). ``0`` selects
        the exact legacy boxcar; ``> 0`` selects the pinned truncated-normal
        footprint mass.
    :param site_footprint_m:
        Site footprint ``z`` (Petersen cell size), metres. Sets the +/-z/2
        window integrated over the mapping-error normal on the ``sigma > 0``
        path. Ignored when ``sigma == 0``.
    :param r_sigma_truncation:
        Truncation ``n`` in +/-n-sigma (Petersen p. 819 uses 2).
    :returns:
        ``W_p`` as a ``float64`` array, ``W_p(0) == 1`` pinned, monotone
        non-increasing in ``|r|``, and exactly ``0`` for ``|r| > z/2 + n*sigma``.
    """
    r = np.abs(np.asarray(r, dtype=np.float64))
    sigma = float(r_sigma_km)

    if sigma == 0.0:
        # Legacy boxcar, bit-for-bit: the principal zone is |r| <= h. The
        # complementary distributed weight 1 - W_p reproduces the historical
        # ~mask_principal split exactly (docs section 2, verification V1).
        return (r <= float(r_threshold_km)).astype(np.float64)

    z = float(site_footprint_m) / 1000.0
    half = z / 2.0
    n = float(r_sigma_truncation)

    def _footprint_mass(rr):
        # Standard-normal mass captured by the +/-z/2 footprint window centred
        # at rr/sigma, with the tails truncated at +/-n (Petersen p. 819).
        upper = np.clip((rr + half) / sigma, -n, n)
        lower = np.clip((rr - half) / sigma, -n, n)
        return ndtr(upper) - ndtr(lower)

    # Pinning: W_p(0) == 1. M(0) normalisation replaces the explicit
    # truncated-normal renormalisation (docs section 2, decision D2).
    m0 = _footprint_mass(0.0)
    return _footprint_mass(r) / m0
