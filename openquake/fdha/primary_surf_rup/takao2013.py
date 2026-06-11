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
Module :mod:`openquake.fdha.primary_surf_rup.takao2013` implements
the model of Takao et al. (2013) in :class:`Takao2013PrimarySR`

Supported Fault Styles: Reverse & Strike‑Slip
"""

import numpy as np
from openquake.fdha.primary_surf_rup.base import BasePrimarySurfRup

class Takao2013PrimarySR(BasePrimarySurfRup):
    """Principal surface-rupture probability model of Takao et al. (2013).

    Logistic model of the probability of principal surface rupture as a
    function of magnitude, for reverse and strike-slip faulting.

    References
    ----------
    Takao, M., Tsuchiyama, J., Annaka, T., & Kurita, T. (2013). Application of
    probabilistic fault displacement hazard analysis in Japan. Journal of Japan
    Association for Earthquake Engineering, 13(1), 17-36.
    https://doi.org/10.5610/jaee.13.17
    """

    def get_prob(
        self,
        mag: float
    ) -> float:
        """
        Model of Takao et al. (2013) for the probability of surface
        rupture for rupture with reverse and strike-slip rupturing mechanism.

        :param mag: float or array-like, earthquake magnitude(s)
        :return: probability or array of probabilities
        """
        m = np.asarray(mag, dtype=float)

        fx = -32.03 + 4.9 * m
        prob = np.exp(fx) / (1.0 + np.exp(fx))
        
        return prob.item() if prob.shape == () else prob