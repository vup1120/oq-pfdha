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
* ``sigma > 0`` is that same boxcar(h) smoothed by the mapping-error normal
  ``N(0, sigma)`` truncated at +/-n-sigma, pinned. This carries mapping
  accuracy: Petersen et al. (2011) give two-sided mapping errors per accuracy
  class (Tables 2-3, p. 810) and place the principal rupture "within 2
  standard deviations of the mapped fault trace" with normal weights in their
  worked example (p. 819) -- hence the +/-n-sigma truncation with n = 2 by
  default. The mass is pinned so that ``W_p(0) == 1`` -- Petersen p. 819 takes
  "the probability of fault rupture on the fault to be one"; the M(0)
  normalisation subsumes the explicit truncated-normal renormalisation.

Because the smoothing window is ``h`` itself (NOT the site footprint ``z`` --
``z`` belongs exclusively to the distributed rupture probability, Petersen
electronic supplement / Tables 4-5), the ``sigma -> 0+`` limit collapses
smoothly onto the ``sigma == 0`` boxcar: there is no discontinuity between the
two paths. Validated against the Petersen Fig. 9c/9d anchors (design doc,
section 2): widths scale with sigma and the profile toe sits at
``h + n*sigma``.
"""

import numpy as np
from scipy.special import ndtr

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
        Principal-zone half-width ``h``: the boxcar half-width at
        ``sigma == 0`` and the smoothing-window half-width at ``sigma > 0``.
        Using the same ``h`` on both paths makes ``sigma -> 0+`` collapse
        smoothly onto the ``sigma == 0`` boxcar.
    :param r_sigma_km:
        Two-sided mapping-accuracy sigma (Petersen Tables 2-3). ``0`` selects
        the exact legacy boxcar; ``> 0`` smooths that boxcar with the
        truncated mapping-error normal, pinned.
    :param r_sigma_truncation:
        Truncation ``n`` in +/-n-sigma (Petersen p. 819: the principal fault
        occurs "within 2 standard deviations" -- n = 2).
    :returns:
        ``W_p`` as a ``float64`` array, ``W_p(0) == 1`` pinned, monotone
        non-increasing in ``|r|``, and exactly ``0`` for ``|r| > h + n*sigma``.
    """
    r = np.abs(np.asarray(r, dtype=np.float64))
    sigma = float(r_sigma_km)
    h = float(r_threshold_km)

    if sigma == 0.0:
        # Legacy boxcar, bit-for-bit: the principal zone is |r| <= h. The
        # complementary distributed weight 1 - W_p reproduces the historical
        # ~mask_principal split exactly (docs section 2, verification V1).
        return (r <= h).astype(np.float64)

    n = float(r_sigma_truncation)

    def _zone_mass(rr):
        # Mass of the mapping-error normal N(0, sigma), truncated at
        # +/-n-sigma (Petersen p. 819), captured by the principal zone
        # [rr - h, rr + h] -- i.e. boxcar(h) convolved with the truncated
        # normal. The window is h, NOT the site footprint z: z belongs
        # exclusively to the distributed rupture probability (Petersen
        # electronic supplement / Tables 4-5).
        upper = np.clip((rr + h) / sigma, -n, n)
        lower = np.clip((rr - h) / sigma, -n, n)
        return ndtr(upper) - ndtr(lower)

    # Pinning: W_p(0) == 1. M(0) normalisation replaces the explicit
    # truncated-normal renormalisation (docs section 2, decision D2).
    m0 = _zone_mass(0.0)
    return _zone_mass(r) / m0
