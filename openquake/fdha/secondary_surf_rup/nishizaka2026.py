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
Module :mod:`openquake.fdha.secondary_surf_rup.nishizaka2026` implements
the far-field surface-rupture probability model of Nishizaka et al. (2026)
in :class:`Nishizaka2026SecondarySR`.

Supported fault styles: strike-slip only.

The model gives the conditional probability of surface rupture ``P2s`` in a
square unit cell, for *all* surface ruptures without classifying them as
principal or distributed (their Equation 2 restricted to the conventional
terms):

    P2s = exp(z) / (1 + exp(z))
    z = c1 + c2*ln(r1 + c3) + c4*ln(r2 + c5)

where ``r1`` is the 3D distance from the subsurface earthquake source fault
and ``r2`` is the shortest surface distance from mapped, pre-existing active
faults (both in km). Coefficients were regressed separately for the 2016
Kumamoto and 2019 Ridgecrest earthquakes and, within each event, for the
"plain" and "mountain" sides of the source fault (Table 1 of the paper,
unit cell size 250 m).

Reference
---------
Nishizaka, N., Onishi, K., Ikeda, M., Si, H., Yamamoto, K., & Tsuji, T.
(2026). Characteristics of far-field surface ruptures caused by two recent
strike-slip earthquakes: insights into fault displacement prediction.
Seismological Research Letters. https://doi.org/10.1785/0220250293
"""

import numpy as np
from scipy.special import expit
from openquake.fdha.secondary_surf_rup.base import BaseSecondarySurfRup


class Nishizaka2026SecondarySR(BaseSecondarySurfRup):
    """
    Far-field surface-rupture probability model of Nishizaka et al. (2026)
    for strike-slip earthquakes, as a function of the distance from the
    earthquake source fault (``r``, km) and the distance from mapped,
    pre-existing active faults (``r2``, km).

    Notes
    -----
    - ``P2s`` covers *all* surface ruptures (principal and distributed are
      not separated). Combining this model with a separate principal-rupture
      model may double-count near-fault hazard; it is intended for far-field
      sites (the paper focuses on distances greater than about 3 km).
    - The model is independent of earthquake magnitude.
    - Coefficients are event-specific (epistemic alternatives): ``region``
      selects the 2016 Kumamoto or 2019 Ridgecrest regression, and ``side``
      selects the plain or mountain side of the source fault.
    - Only the 250 m unit cell coefficients of the paper's Table 1 are
      bundled; the supplemental Table S1 (other cell sizes) is not included.
    """

    #: Regression coefficients (c1, c2, c3, c4, c5) from Table 1 of
    #: Nishizaka et al. (2026), unit cell size 250 m, distances in km.
    COEFFS = {
        250: {
            "kumamoto": {
                "plain": {"c1": 5.05, "c2": -3.32, "c3": 7.66,
                          "c4": -0.569, "c5": 2.20e-3},
                "mountain": {"c1": 13.4, "c2": -7.51, "c3": 7.66,
                             "c4": -0.823, "c5": 2.20e-3},
            },
            "ridgecrest": {
                "plain": {"c1": 0.238, "c2": -2.58, "c3": 1.19e-4,
                          "c4": -0.243, "c5": 9.98e-5},
                "mountain": {"c1": 0.105, "c2": -3.14, "c3": 1.19e-4,
                             "c4": -0.433, "c5": 9.98e-5},
            },
        },
    }

    def get_prob(self, r, r2, region="kumamoto", side="plain",
                 cell_size=250):
        """
        Conditional probability of surface rupture ``P2s`` in a unit cell.

        :param r:
            Distance r1 from the earthquake source fault in km (scalar or
            array). The paper measures r1 as the 3D distance from the
            subsurface source fault; when driven by the hazard calculator
            the rupture-trace distance is used as an approximation.
        :param r2:
            Shortest surface distance from mapped, pre-existing active
            faults in km (scalar or array broadcastable against ``r``).
        :param region:
            Coefficient set: ``"kumamoto"`` (2016 Kumamoto) or
            ``"ridgecrest"`` (2019 Ridgecrest).
        :param side:
            ``"plain"`` or ``"mountain"`` side of the source fault.
        :param cell_size:
            Unit cell size in metres. Currently only 250 is available.
        :return:
            Probability (0-1), scalar or array following the inputs.
        """
        if cell_size not in self.COEFFS:
            raise ValueError(
                f"Invalid cell_size '{cell_size}'. "
                f"Available: {sorted(self.COEFFS.keys())} m")
        region_l = str(region).strip().lower()
        if region_l not in self.COEFFS[cell_size]:
            raise ValueError(
                f"Invalid region '{region}'. "
                f"Accepted: {sorted(self.COEFFS[cell_size].keys())}")
        side_l = str(side).strip().lower()
        if side_l not in self.COEFFS[cell_size][region_l]:
            raise ValueError(
                f"Invalid side '{side}'. Accepted: "
                f"{sorted(self.COEFFS[cell_size][region_l].keys())}")

        c = self.COEFFS[cell_size][region_l][side_l]
        r1_arr = np.asarray(r, dtype=float)
        r2_arr = np.asarray(r2, dtype=float)
        if np.any(r1_arr < 0) or np.any(r2_arr < 0):
            raise ValueError("Distances r and r2 must be non-negative")

        z = (c["c1"] + c["c2"] * np.log(r1_arr + c["c3"])
             + c["c4"] * np.log(r2_arr + c["c5"]))
        prob = expit(z)  # numerically stable logistic
        return prob.item() if prob.shape == () else prob
