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
Module :mod:`openquake.fdha.secondary_surf_rup.petersen2011` 
"""

import numpy as np
from openquake.fdha.secondary_surf_rup.base import BaseSecondarySurfRup


class Petersen2011SecondarySR(BaseSecondarySurfRup):
    """
    Implementation of the Petersen et al. (2011) model for strike-slip faults
    with different pixel sizes
    """
    
    # Define coefficients for different pixel sizes
    # Cell size parameters from Table 4 (Page 812, Petersen et al., 2011)
    CELL_SIZES = {
        25: {"a": -1.1470, "b": 2.1046, "sigma": 1.2508},  # 25 x 25 m
        50: {"a": -0.9000, "b": 0.9866, "sigma": 1.1470},  # 50 x 50 m
        100: {"a": -1.0114, "b": 2.5572, "sigma": 1.0917},  # 100 x 100 m
        150: {"a": -1.0934, "b": 3.5526, "sigma": 1.0188},  # 150 x 150 m
        200: {"a": -1.1538, "b": 4.2342, "sigma": 1.0177},  # 200 x 200 m
    }

    # Near-field interpolation points from Table 5 (Page 812, Petersen et al., 2011)
    NEAR_FIELD_POINTS = {
        25: {"p0": 0.74541, "p1": 0.078690, "p2": 0.020108, "r1": 100, "r2": 200},  # 25 x 25 m
        50: {"p0": 0.87162, "p1": 0.048206, "p2": 0.026177, "r1": 100, "r2": 200},  # 50 x 50 m
        100: {"p0": 0.90173, "p1": 0.18523, "p2": 0.066354, "r1": 100, "r2": 200},  # 100 x 100 m
        150: {"p0": 0.87394, "p1": 0.19592, "p2": 0.070477, "r1": 150, "r2": 300},  # 150 x 150 m
        200: {"p0": 0.92483, "p1": 0.18975, "p2": 0.074709, "r1": 200, "r2": 400},  # 200 x 200 m
    }

    def get_prob(self, r, cell_size=25, version="default"):
        """
        Calculate the probability of distributed-fault surface rupture as a function of distance,
        cell size, and version, per Petersen et al. (2011).

        :param r:
            Distance from the principal fault trace in kilometers (up to 2 km recommended).
        :param cell_size:
            Size of the cell in meters (25, 50, 100, 150, or 200 m; default: 25 m).
        :param version:
            Model version (case-insensitive). Options: 'default' (power function, Page 812, Table 4),
            'near_field' (interpolated near-field, Page 819, electronic supplement). Default: 'default'.
        :returns:
            Probability of rupture (float or array, 0–1) for the given distance, cell_size, and version.
        :raises ValueError:
            If r is negative or exceeds 2000 m, cell_size is invalid, or version is invalid.
        :notes:
            - Uses power function from Table 4 (Page 812) for 'default' (far-field probabilities).
            - Uses near-field interpolation from Page 819 and electronic supplement for 'near_field' (r < r1).
            - No magnitude dependence, per Petersen et al. (2011, Page 818).
            - Limited to 2 km distance from principal fault; no triggered ruptures included.
            - Near-field points derived from text and electronic supplement, not explicitly numbered as Table 5.
        """
        # Validate inputs
        version = version.lower()
        valid_versions = ["default", "near_field"]
        if version not in valid_versions:
            raise ValueError(f"Invalid version '{version}'. Accepted values are: {', '.join(valid_versions)}")

        # Handle scalar/array input following Youngs2003 pattern
        r = np.asarray(r)
        if r.ndim == 0:
            r = np.array([r])
            was_scalar = True
        else:
            was_scalar = False

        # Convert distance from km to m
        r = r * 1000
        
        # Validate distance range
        #if not (r >= 0).all() or (r > 2000).any():
        #    raise ValueError("Distance r must be non-negative and ≤ 2000 m")
            
        if isinstance(cell_size, str):
            try:
                cell_size = int(cell_size)
            except ValueError:
                raise ValueError(f"Cell size must be convertible to an integer")
        if cell_size not in self.CELL_SIZES:
            raise ValueError(f"Cell size must be one of {list(self.CELL_SIZES.keys())} m")

        # Get cell size parameters
        params = self.CELL_SIZES[cell_size]
        a, b, _ = params["a"], params["b"], params["sigma"]
        near_params = self.NEAR_FIELD_POINTS[cell_size]
        
        if version == "default":
            # Far-field power function (Page 819, Eqn 20, Table 4)
            r_safe = np.where(r == 0, 0.1, r)  # Use 0.1 m as minimum distance
            ln_P = a * np.log(r_safe) + b
            P_rupture = np.exp(ln_P)  # Convert ln(P) to probability
            P_rupture = np.clip(P_rupture, 0, 1)  # Ensure probability is in [0, 1]
        else:  # version == "near_field"
            # Near-field interpolation (Page 819, electronic supplement)
            p0, p1, p2 = near_params["p0"], near_params["p1"], near_params["p2"]
            r1, r2 = near_params["r1"], near_params["r2"]

            # Initialize with default probability (p0)
            P_rupture = np.full_like(r, p0, dtype=float)
            
            # Linear interpolation for r < r1
            mask_near = r < r1
            if np.any(mask_near):
                # Interpolate between p0 and p1 at r1
                P_rupture[mask_near] = p0 + (p1 - p0) * (r[mask_near] / r1)
                
            # Linear interpolation for r1 <= r <= r2
            mask_mid = (r >= r1) & (r <= r2)
            if np.any(mask_mid):
                # Interpolate between p1 and p2
                P_rupture[mask_mid] = p1 + (p2 - p1) * ((r[mask_mid] - r1) / (r2 - r1))
                
            # Use power function for far field (r > r2)
            mask_far = r > r2
            if np.any(mask_far):
                r_safe = np.where(r[mask_far] == 0, 0.1, r[mask_far])
                ln_P_far = a * np.log(r_safe) + b
                P_rupture[mask_far] = np.exp(ln_P_far)
                
            P_rupture = np.clip(P_rupture, 0, 1)  # Ensure probability is in [0, 1]

        # Return scalar if input was scalar, following Youngs2003 pattern
        if was_scalar:
            return float(P_rupture[0])
        else:
            return P_rupture


class Petersen2011SecondarySR_default(Petersen2011SecondarySR):
    """Petersen et al. (2011) distributed-rupture model fixed to the 'default' variant."""

    def get_prob(self, r, cell_size=25, version="default"):
        """Return distributed surface-rupture probability using the 'default' variant."""
        return super().get_prob(r=r, cell_size=cell_size, version="default")
