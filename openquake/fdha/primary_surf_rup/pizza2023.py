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
Module :mod:`openquake.fdha.primary_surf_rup.pizza2023`
"""

import numpy as np
from openquake.fdha.primary_surf_rup.base import BasePrimarySurfRup


class Pizza2023PrimarySR(BasePrimarySurfRup):
    """Principal surface-rupture probability model of Pizza et al. (2023).

    Logistic model of the probability of principal surface rupture as a
    function of magnitude, with coefficients for the ``all``, ``normal``,
    ``reverse``, and ``strike-slip`` faulting styles.

    References
    ----------
    Pizza, M., et al. (2023). Probability of surface rupture for normal and
    other faulting styles.
    """

    def get_prob(self, mag, style="all"):
        """
        Model of Pizza et al., 2023 for the probability of surface rupture
        based on rupture mechanism and earthquake magnitude.

        This model estimates the likelihood of surface rupture for four faulting 
        styles: "all", "normal", "reverse", and "strike-slip." The probability is 
        calculated using a logistic regression formula, with coefficients varying 
        based on the selected faulting style.

        :param mag:
            The magnitude of the seismic event (float). Larger magnitudes generally 
            lead to higher probabilities of surface rupture.

        :param style:
            The faulting style (string). Accepted values are:
            - "all": Represents a combination of all faulting styles.
            - "normal": Represents normal faulting mechanisms.
            - "reverse": Represents reverse faulting mechanisms.
            - "strike-slip": Represents strike-slip faulting mechanisms.
            
            Default is "all." If an unsupported style is provided, a ValueError will 
            be raised.

        :raises ValueError:
            If an invalid style is provided, the function raises an error indicating 
            the acceptable faulting styles.

        :return:
            The probability of surface rupture as a float value between 0 and 1, 
            calculated using a logistic regression model.
        """

        # Define the accepted style of faultings
        accepted_styles = ["all", "normal", "reverse", "strike-slip"]
        
        # Validate the style
        if style not in accepted_styles:
            raise ValueError(
                f"Invalid style '{style}'. Accepted values are: {', '.join(accepted_styles)}"
            )

        m = np.asarray(mag, dtype=float)
        
        if style == "all":
            a = -14.47
            b = 2.177
        elif style == "normal":
            a = -13.5
            b = 2.159
        elif style == "reverse":
            a = -10.75
            b = 1.427
        elif style == "strike-slip":
            a = -28.56
            b = 4.436
        else:
            raise ValueError(f"Invalid style '{style}'.")
        fx = a + b * m
        prob = np.exp(fx) / (1.0 + np.exp(fx))
        
        return prob.item() if prob.shape == () else prob
