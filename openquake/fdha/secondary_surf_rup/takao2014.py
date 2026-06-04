# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#
# Copyright (C) 2012-2024 GEM Foundation
#
# OpenQuake is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenQuake is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenQuake. If not, see <http://www.gnu.org/licenses/>.

"""
Module :mod:`openquake.fdha.secondary_surf_rup.takao2014` 
"""

import numpy as np
from openquake.fdha.secondary_surf_rup.base import BaseSecondarySurfRup


class Takao2014SecondarySR(BaseSecondarySurfRup):
    """Implementation of the Takao et al. (2014) model for reverse and strike-slip faults"""
    
    def get_prob(self, r: float, pixel_size: int = 100):
        """
        Calculates probability of surface rupture for reverse and strike-slip faults
        
        :param r:
            The distance to the principal fault in km
        :param pixel_size:
            Size of the square analysis pixel in meters.
            Valid options are: 500, 250, 100, 50
        :return:
            Probability of surface rupture
        """
        # Coefficients for different pixel sizes
        coefficients = {
            500: (-3.859, -1.499, 0.2),
            250: (-4.903, -1.459, 0.2),
            100: (-6.135, -1.427, 0.2),
            50: (-6.988, -1.410, 0.2)
        }
        
        if pixel_size not in coefficients:
            raise ValueError(f"Invalid pixel size. Must be one of {list(coefficients.keys())} meters")
            
        C1, C2, C3 = coefficients[pixel_size]
        fx = C1 + C2 * np.log(r + C3)

        # Calculate probability using the formula: exp(-C1·ln(r·1000) + C2)
        prob = np.exp(fx) / (1 + np.exp(fx))

        return prob