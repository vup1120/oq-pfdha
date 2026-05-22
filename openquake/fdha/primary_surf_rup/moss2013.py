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
Module :mod:`openquake.hazardlib.fdha.Moss2013` 
"""

import numpy as np
from openquake.fdha.primary_surf_rup.base import BasePrimarySurfRup


class Moss2013PrimarySR(BasePrimarySurfRup):
    """Principal surface-rupture probability model of Moss et al. (2013).

    Logistic model of the probability of principal surface rupture as a
    function of magnitude, faulting style, and site Vs30.

    References
    ----------
    Moss, R.E.S., et al. (2013). Probabilistic fault displacement hazard
    analysis for reverse faults.
    """

    def get_prob(self, mag: float, style, vs30: float) -> float:
        """
        Model of Moss et al. (2013) for the probability of surface rupture
        based on the rupture mechanism (faulting style) and site conditions 
        (shear-wave velocity).

        This model estimates the likelihood of surface rupture for two styles 
        of faulting: "reverse" and "strike-slip." The probability is calculated 
        as a function of event magnitude (mag) and the time-averaged shear-wave 
        velocity to a depth of 30 meters (Vs30), which characterizes the site as 
        soft or stiff soil.

        :param mag:
            The magnitude of the seismic event (float). Higher magnitudes generally 
            result in greater probabilities of surface rupture.
            
        :param style:
            The faulting style (string). Accepted values are:
            - "reverse": Reverse faulting mechanism.
            - "strike-slip": Strike-slip faulting mechanism.
            
            Default is "reverse." If an unsupported style is provided, a ValueError 
            will be raised.

        :param vs30:
            The time-averaged shear-wave velocity (float) to a depth of 30 meters 
            (Vs30). This parameter distinguishes between stiff soil or rock (Vs30 > 600 m/s) 
            and soft soil (Vs30 ≤ 600 m/s). The rupture probability varies based on 
            this classification.
            
        :raises ValueError:
            If an invalid style is provided, the function raises an error indicating 
            the acceptable faulting styles.

        :return:
            The probability of surface rupture as a float value between 0 and 1.

        :example:
            # Example usage
            model = Moss2013Primary()
            probability = model.get_prob(mag=7.5, style="reverse", vs30=500)
            print(f"Probability of surface rupture: {probability}")
        """

        # Define the accepted style of faultings
        accepted_styles = ["reverse", "strike-slip"]
        
        # Validate the style
        if style not in accepted_styles:
            raise ValueError(
                f"Invalid style '{style}'. Accepted values are: {', '.join(accepted_styles)}"
            )

        m = np.asarray(mag, dtype=float)
        if style == 'reverse':
            if vs30 > 600:
                a, b = 13.9745, 2.1395
                prob = 1.0 / (1.0 + np.exp(a - b * m))
            else:
                a, b = 6.2548, 0.8308
                prob = 1.0 / (1.0 + np.exp(a - b * m))
        else:  # strike-slip
            if vs30 > 600:
                a, b = 11.4071, 1.8465
                prob = 1.0 / (1.0 + np.exp(a - b * m))
            else:
                a, b = 12.2908, 1.9520
                prob = 1.0 / (1.0 + np.exp(a - b * m))
        return prob.item() if prob.shape == () else prob