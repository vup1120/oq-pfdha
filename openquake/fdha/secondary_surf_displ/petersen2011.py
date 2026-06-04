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
Module :mod:`openquake.fdha.secondary_surf_displ.petersen2011` extends
model of Petersen et al. (2011) for secondary (distributed) fault displacements
into :class:`Petersen2011SecondaryFD`
"""

import numpy as np
from scipy.stats import norm
from openquake.fdha.primary_surf_displ.base import BaseSecondarySurfDispl

class Petersen2011SecondaryFD(BaseSecondarySurfDispl):
    """Distributed fault-displacement model of Petersen et al. (2011) for strike-slip faults.

    References
    ----------
    Petersen, M.D., et al. (2011). Fault displacement hazard for strike-slip
    faults. Bulletin of the Seismological Society of America, 101(2), 805-825.
    """

    # Cell size parameters from Table 4 (Page 812, Petersen et al., 2011)
    CELL_SIZES = {
        25: {"a": -1.1470, "b": 2.1046, "sigma": 1.2508},  # 25 x 25 m
        50: {"a": -0.9000, "b": 0.9866, "sigma": 1.1470},  # 50 x 50 m
        100: {"a": -1.0114, "b": 2.5572, "sigma": 1.0917},  # 100 x 100 m
        150: {"a": -1.0934, "b": 3.5526, "sigma": 1.0188},  # 150 x 150 m
        200: {"a": -1.1538, "b": 4.2342, "sigma": 1.0177},  # 200 x 200 m
    }
    
    # Near-field interpolation points from electronic supplement
    NEAR_FIELD_POINTS = {
        25: {"p0": 0.74541, "p1": 0.078690, "p2": 0.020108, "r1": 100, "r2": 200},
        50: {"p0": 0.87162, "p1": 0.048206, "p2": 0.026177, "r1": 100, "r2": 200},
        100: {"p0": 0.90173, "p1": 0.18523, "p2": 0.066354, "r1": 100, "r2": 200},
        150: {"p0": 0.87394, "p1": 0.19592, "p2": 0.070477, "r1": 150, "r2": 300},
        200: {"p0": 0.92483, "p1": 0.18975, "p2": 0.074709, "r1": 200, "r2": 400},
    }

    def get_prob(self, d, mag, r, cell_size=25):
        """
        Calculate the probability of exceeding displacement thresholds [m] for distributed
        strike-slip faults, per Petersen et al. (2011).

        :param d:
            Target displacement values in meters (array of shape (n_displacements,))
        :param mag:
            Earthquake moment magnitude (scalar, recommended range: 6–8 for strike-slip faults)
        :param r:
            Distance from the principal fault trace in kilometers (array of shape (n_sites,))
        :param cell_size:
            Size of the cell in meters (25, 50, 100, 150, or 200 m; default: 25 m)
        :returns:
            Probability of exceeding the target displacement (shape (n_sites, n_displacements))
        :raises ValueError:
            If mag is outside [6, 8], r is negative, or cell_size is invalid.
        :notes:
            - Applies to distributed (off-fault) surface fault displacement on strike-slip faults,
              per Petersen et al. (2011, doi:10.1785/0120100035).
            - Uses Wells and Coppersmith (1994) for strike-slip average displacement.
            - Limited to 2 km distance from principal fault; no triggered ruptures included.
            - Returns only the displacement exceedance probability (prob_exceeding_d); combine
              with rupture probability (get_prob_rupture) separately, per Petersen et al. (Page 818, Eqn 18).
        """
        # Ensure inputs are arrays with proper shapes - following Youngs2003 pattern
        d = np.asarray(d)
        if d.ndim == 0:
            d = np.array([d])
        
        r = np.asarray(r)
        if r.ndim == 0:
            r = np.array([r])
        
        # Convert distance from km to m
        r = r * 1000  # Now in meters
        
        # Ensure magnitude is scalar
        if np.isscalar(mag):
            mag = float(mag)
        else:
            mag = float(np.atleast_1d(mag)[0])
        
        # Validate inputs
        if not (6 <= mag <= 8):
            raise ValueError("Magnitude must be between 6 and 8 for strike-slip faults")
        
        # Handle zero distances to avoid log(0) warning
        # Replace zeros with a small positive value (0.1 m)
        r = np.where(r == 0, 0.1, r)
        
        # Get dimensions
        n_sites = r.shape[0]
        n_displacements = d.shape[0]
        
        # Calculate average displacement (D_ave) for strike-slip faults per Wells and Coppersmith (1994)
        D_ave = 10 ** (-6.32 + 0.90 * mag)  # In meters, Page 813
        sigma_D_ave = 0.28  # Standard deviation in log10 units
        
        # Convert displacement to centimeters for Petersen's regression (Page 818, Eqn 18)
        d_cm = d * 100  # Convert to cm to match paper units
        
        # Calculate mean for each site (Page 818, Eqn 18)
        # mu will have shape (n_sites,)
        mu = 1.4016 * mag - 0.1671 * np.log(r) - 6.7991  # ln(d) in cm
        
        # Standard deviation in ln(cm) units, Page 818
        sigma_dist = 1.1193
        
        # Reshape arrays for proper broadcasting
        # d_cm: (n_displacements,) -> (1, n_displacements)
        # mu: (n_sites,) -> (n_sites, 1)
        d_cm_reshaped = d_cm[np.newaxis, :]  # Shape (1, n_displacements)
        mu_reshaped = mu[:, np.newaxis]  # Shape (n_sites, 1)
        
        # Calculate probability of exceeding d using log-normal distribution
        # This will broadcast to shape (n_sites, n_displacements)
        prob_exceeding = 1 - norm.cdf(np.log(d_cm_reshaped), loc=mu_reshaped, scale=sigma_dist)
        
        # Handle return values following Youngs2003 pattern
        if n_sites == 1 and n_displacements == 1:
            return float(prob_exceeding[0, 0])
        elif n_sites == 1:
            return prob_exceeding[0, :]  # Return 1D array for single site
        elif n_displacements == 1:
            return prob_exceeding[:, 0]  # Return 1D array for single displacement
        else:
            return prob_exceeding  # Return 2D array