"""
Hazard-curve interpolation helpers (POE and return-period lookups).

The interpolation itself is OpenQuake's own
:func:`openquake.hazardlib.map_array.compute_hazard_maps` (exp-log
interpolation with the engine's EPSILON cutoff and extrapolate-to-zero
behaviour below the curve minimum); this module only adapts it to the
FDHA calling convention.
"""
import numpy as np

from openquake.hazardlib.map_array import compute_hazard_maps


def get_map_from_curves(
    imls: np.ndarray,
    poes_matrix: np.ndarray,
    pex: float,
) -> np.ndarray:
    """
    Convert hazard curves to displacement map values for a target
    probability (or annual rate) of exceedance.

    Thin wrapper around :func:`openquake.hazardlib.map_array.
    compute_hazard_maps` for the single-target case, so FDHA map values
    interpolate exactly like OpenQuake engine hazard maps.

    Parameters:
    -----------
    imls : array-like
        Intensity measure levels (displacement values)
    poes_matrix : array-like
        Matrix of probabilities/rates of exceedance for each site and
        displacement level, shape (n_sites, n_imls)
    pex : float
        Target probability/rate of exceedance (e.g., 1/return_period)

    Returns:
    --------
    array
        Displacement values corresponding to the target for each site
    """
    hazard_maps = compute_hazard_maps(poes_matrix, imls, [pex])
    # Single target: return the only column.
    return hazard_maps[:, 0]
