"""
Hazard-curve interpolation helpers (POE and return-period lookups),
following OpenQuake engine conventions.
"""
import numpy as np
import warnings
from typing import Union, List

# Use the same EPSILON as OpenQuake Engine
EPSILON = 1e-30

def interpolate_poe(
    imls: np.ndarray,
    poes: np.ndarray,
    target: float
) -> float:
    """
    Interpolate probability of exceedance for a given displacement level.
    
    Parameters:
    -----------
    imls : array-like
        Displacement levels (intensity measure levels)
    poes : array-like
        Probabilities of exceedance corresponding to each displacement level
    target : float
        Target displacement level to interpolate for
        
    Returns:
    --------
    float
        Interpolated probability of exceedance for the target displacement
    """
    if target <= imls[0]:
        return poes[0]
    if target >= imls[-1]:
        return poes[-1]
    return np.exp(np.interp(np.log(target), np.log(imls), np.log(poes)))

def compute_hazard_maps(
    curves: np.ndarray,
    imls: np.ndarray,
    poes: Union[float, List[float], np.ndarray]
) -> np.ndarray:
    """
    OpenQuake Engine's standard implementation for computing hazard maps from curves.
    
    Given a set of hazard curve poes, interpolate hazard maps at the specified poes.
    
    Parameters:
    -----------
    curves : array-like
        Array of floats of shape N x L. Each row represents a curve, where the
        values in the row are the PoEs (Probabilities of Exceedance)
        corresponding to the imls. Each curve corresponds to a geographical location.
    imls : array-like
        Intensity Measure Levels associated with these hazard curves. Type
        should be an array-like of floats.
    poes : array-like
        Value(s) on which to interpolate a hazard map from the input curves.
        
    Returns:
    --------
    array
        An array of shape N x P, where N is the number of curves and P the
        number of poes.
    """
    P = len(poes)
    N, L = curves.shape  # number of levels
    if L != len(imls):
        raise ValueError('The curves have %d levels, %d were passed' %
                         (L, len(imls)))

    log_poes = np.log(poes)
    hmap = np.zeros((N, P))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # avoid RuntimeWarning: divide by zero for zero levels
        imls = np.log(np.array(imls[::-1]))
    for n, curve in enumerate(curves):
        # the hazard curve, having replaced the too small poes with EPSILON
        log_curve = np.log([max(poe, EPSILON) for poe in curve[::-1]])
        for p, log_poe in enumerate(log_poes):
            if log_poe > log_curve[-1]:
                # special case when the interpolation poe is bigger than the
                # maximum, i.e the iml must be smaller than the minimum;
                # extrapolate the iml to zero
                # then the hmap goes automatically to zero
                pass
            else:
                # exp-log interpolation, to reduce numerical errors
                hmap[n, p] = np.exp(np.interp(log_poe, log_curve, imls))
    return hmap

def get_map_from_curves(
    imls: np.ndarray,
    poes_matrix: np.ndarray,
    pex: float
) -> np.ndarray:
    """
    Convert hazard curves to displacement map values for a target probability of exceedance.
    
    This function uses OpenQuake Engine's standard interpolation method to ensure
    consistency with the official OpenQuake implementation.
    
    Parameters:
    -----------
    imls : array-like
        Intensity measure levels (displacement values)
    poes_matrix : array-like
        Matrix of probabilities of exceedance for each site and displacement level
    pex : float
        Target probability of exceedance (e.g., 1/return_period)
    
    Returns:
    --------
    array
        Displacement values corresponding to the target probability for each site
    """
    # Use OpenQuake's standard compute_hazard_maps function
    # This ensures we get the exact same results as OpenQuake Engine
    hazard_maps = compute_hazard_maps(poes_matrix, imls, [pex])
    
    # Return the first (and only) column since we only have one target probability
    return hazard_maps[:, 0]