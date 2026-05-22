# -*- coding: utf-8 -*-
"""
Probability utility functions for Monte Carlo reduction and array normalization.

These functions handle reduction of Monte Carlo samples and standardization
of probability arrays to consistent shapes for hazard calculations.
"""
import numpy as np
from typing import Union, Dict, Any, Optional, Tuple


def _reduce_mc(
    prob: Union[float, np.ndarray, None],
    method: str = "median",
    q: float = 50
) -> Union[float, np.ndarray]:
    """
    Reduce Monte Carlo (MC) samples to a single value per site.

    Accepts:
        - scalar
        - (n_mc,)
        - (n_sites,)
        - (n_sites, n_mc)
    Returns:
        - scalar or (n_sites,) vector

    If prob is None, return 1.0 (as requested).

    Parameters
    ----------
    prob : array-like or scalar or None
        Probability values, possibly with Monte Carlo samples
    method : str, default="median"
        Reduction method: "median", "mean", or "percentile"
    q : float, default=50
        Percentile value (0-100) if method="percentile"

    Returns
    -------
    float or ndarray
        Reduced probability value(s)
    """
    if prob is None:
        return 1.0
    
    arr = np.asarray(prob)
    
    # Fast path: scalar input
    if arr.ndim == 0:
        return float(arr)
    
    # For 1D arrays, reduce over the only axis
    if arr.ndim == 1:
        if method == "percentile":
            qq = float(q)
            if 0.0 <= qq <= 1.0:
                qq *= 100.0
            return np.percentile(arr, qq, axis=0)
        elif method == "mean":
            return np.mean(arr, axis=0)
        else:  # "median" by default
            return np.median(arr, axis=0)
    
    # For 2D arrays, reduce over the last axis (MC dimension)
    if method == "percentile":
        qq = float(q)
        if 0.0 <= qq <= 1.0:
            qq *= 100.0
        return np.percentile(arr, qq, axis=-1)
    elif method == "mean":
        return np.mean(arr, axis=-1)
    else:  # "median" by default
        return np.median(arr, axis=-1)


def _reduce_over_axes(
    arr: np.ndarray,
    axes: Union[int, Tuple[int, ...]],
    method: str = "median",
    q: float = 50
) -> np.ndarray:
    """
    Reduce array over given axes using the chosen statistic.

    Parameters
    ----------
    arr : array-like
        Input array
    axes : int or tuple of int
        Axes over which to reduce
    method : str, default="median"
        Reduction method: "median", "mean", or "percentile"
    q : float, default=50
        Percentile value (0-100) if method="percentile"

    Returns
    -------
    ndarray
        Reduced array
    """
    if method == "percentile":
        qq = float(q)
        if 0.0 <= qq <= 1.0:
            qq *= 100.0
        return np.percentile(arr, qq, axis=axes)
    elif method == "mean":
        return np.mean(arr, axis=axes)
    else:  # "median" by default
        return np.median(arr, axis=axes)


def _to_sites_x_displ(
    arr: Optional[Union[float, np.ndarray]],
    n_sites: int,
    n_displ: int,
    red_cfg: Dict[str, Any]
) -> Optional[np.ndarray]:
    """
    Normalize any probability array to shape ``(n_sites, n_displ)``.

    This utility reduces any present Monte Carlo (MC) dimension by the chosen
    statistic and standardizes array layout for downstream computations.

    Recognized input shapes include:
      - ``(n_sites, n_displ)``: returned as-is
      - ``(n_displ, n_sites)``: transposed
      - ``(n_displ,)``: broadcast across sites
      - ``(n_sites,)``: broadcast across displacements
      - ``(n_displ, n_mc)`` or ``(n_mc, n_displ)``: reduce MC, then broadcast to sites
      - ``(n_sites, n_mc)`` or ``(n_mc, n_sites)``: reduce MC, then broadcast to displacements
      - higher-D arrays with recognizable site/displ axes: move to ``(site, displ, extra...)`` and
        reduce over the remaining axes
      - otherwise: reduce to scalar and broadcast

    Parameters
    ----------
    arr : array-like or None
        Input probability array
    n_sites : int
        Number of sites
    n_displ : int
        Number of displacement levels
    red_cfg : dict
        Reduction configuration with keys:
        - "method": str, reduction method ("median", "mean", "percentile")
        - "q": float, percentile value if method="percentile"

    Returns
    -------
    ndarray or None
        Normalized array of shape (n_sites, n_displ), or None if input was None
    """
    if arr is None:
        return None
    
    # Precompute reduction parameters once
    method = red_cfg.get("method", "median")
    q_val = red_cfg.get("q", 50)
    
    A = np.asarray(arr)
    
    # Fast path: already correct shape
    if A.shape == (n_sites, n_displ):
        return A
    
    # Fast path: scalar input
    if A.ndim == 0:
        return np.full((n_sites, n_displ), float(A), dtype=A.dtype)
    
    # 1D cases
    if A.ndim == 1:
        if A.size == n_displ:
            # (n_displ,) -> broadcast to (n_sites, n_displ)
            return np.broadcast_to(A.reshape(1, n_displ), (n_sites, n_displ))
        if A.size == n_sites:
            # (n_sites,) -> broadcast to (n_sites, n_displ)
            return np.broadcast_to(A.reshape(n_sites, 1), (n_sites, n_displ))
        # treat as MC vector -> reduce to scalar, broadcast
        val = float(_reduce_mc(A, method=method, q=q_val))
        return np.full((n_sites, n_displ), val, dtype=float)
    
    # 2D cases
    if A.ndim == 2:
        # Transpose case
        if A.shape == (n_displ, n_sites):
            return A.T
        
        # Special case: (1, n_displ) should be broadcast directly to (n_sites, n_displ)
        # This happens when model returns a single site's result that needs to be broadcast
        if A.shape == (1, n_displ):
            return np.broadcast_to(A, (n_sites, n_displ))
        
        # Handle cases where one axis matches n_displ or n_sites
        if A.shape[0] == n_displ:
            # (n_displ, n_mc) -> reduce over MC axis, then broadcast to sites
            red = _reduce_over_axes(A, axes=1, method=method, q=q_val)
            return np.broadcast_to(red.reshape(1, n_displ), (n_sites, n_displ))
        
        if A.shape[1] == n_displ:
            # (n_mc, n_displ) -> reduce over MC axis, then broadcast to sites
            # But skip if shape is (1, n_displ) - already handled above
            if A.shape[0] == 1:
                return np.broadcast_to(A, (n_sites, n_displ))
            red = _reduce_over_axes(A, axes=0, method=method, q=q_val)
            return np.broadcast_to(red.reshape(1, n_displ), (n_sites, n_displ))
        
        if A.shape[0] == n_sites:
            # (n_sites, n_mc) -> reduce over MC axis, then broadcast to displ
            red = _reduce_over_axes(A, axes=1, method=method, q=q_val)
            return np.broadcast_to(red.reshape(n_sites, 1), (n_sites, n_displ))
        
        if A.shape[1] == n_sites:
            # (n_mc, n_sites) -> reduce over MC axis, then broadcast to displ
            red = _reduce_over_axes(A, axes=0, method=method, q=q_val)
            return np.broadcast_to(red.reshape(n_sites, 1), (n_sites, n_displ))
    
    # >=3D: find site/displ axes and reduce the rest
    shape = A.shape
    isite = next((i for i, s in enumerate(shape) if s == n_sites), None)
    idispl = next((i for i, s in enumerate(shape) if s == n_displ and i != isite), None)
    
    if isite is not None and idispl is not None:
        # Move site and displ axes to front, then reduce over remaining axes
        A2 = np.moveaxis(A, (isite, idispl), (0, 1))  # -> (site, displ, extra...)
        if A2.ndim > 2:
            A2 = _reduce_over_axes(A2, axes=tuple(range(2, A2.ndim)),
                                   method=method, q=q_val)
        return A2
    
    # Fallback: reduce everything to a scalar, broadcast
    val = float(_reduce_mc(A, method=method, q=q_val))
    return np.full((n_sites, n_displ), val, dtype=float)











