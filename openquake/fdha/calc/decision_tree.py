# -*- coding: utf-8 -*-
"""
Decision tree for Visini et al. (2025) combinations and simple probability combiner.

Case selection logic (per the paper's example narrative):
  - Case 1: faults capable of Rank 1.5 beneath the site cannot be excluded
            -> use combinations A, B, and C
  - Case 2: no Rank 1.5/3 beneath the site, BUT a Rank 1.5 fault exists within 1 km
            -> use combinations A and B
  - Case 3: Rank 1.5/3 faults are excluded beneath the site and within 1 km
            -> use combination A only
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple
import numpy as np
from openquake.hazardlib.site import SiteCollection
from openquake.fdha.calc.utils.rupture_distance import (
    VectorizedRuptureDistanceCalculator,
)


class SiteFaultClassifier:
    """
    Classify the site into Case 1 / 2 / 3 using Rank 1.5 traces and the
    existing rupture-distance utilities.

    Parameters
    ----------
    sitecol : SiteCollection
        SiteCollection for the site (or set of points).
    rank1p5_surfaces : list
        List of surface-like objects (each must support get_fault_trace() and
        get_min_distance(SiteCollection) -> np.ndarray[km]).
    allow_case1_under_site : bool
        If True, directly returns "case1".
    near_radius_m : float
        Distance threshold in meters for Case 2 (default: 1,000 m).

    Rules
    -----
    case1 -> ['A','B','C']
    case2 -> ['A','B'] if any site point is within near_radius_m of any Rank 1.5 trace
    case3 -> ['A'] otherwise
    """

    def __init__(
        self,
        sitecol: SiteCollection,
        rank1p5_surfaces: Optional[List] = None,
        case: Optional[str] = None,
        near_radius_m: float = 1000.0,
    ):
        if sitecol is None:
            raise ValueError("'sitecol' is required.")

        self.sitecol = sitecol
        self.rank1p5_surfaces = rank1p5_surfaces or []
        self.case = case.lower() if case else None
        self.near_radius_m = float(near_radius_m)

    def classify_site(self) -> str:
        """Return ``"case1"``, ``"case2"`` or ``"case3"``."""
        if self.case:
            return self.case
        # !!!! NEED TO CHECK THIS LOGIC
        if not self.rank1p5_surfaces:
            return "case3"

        # Threshold in kilometers for the existing calculators
        near_thresh_km = self.near_radius_m / 1000.0

        # If ANY rank 1.5 surface is within the threshold of ANY site point => case2
        for surf in self.rank1p5_surfaces:
            calc = VectorizedRuptureDistanceCalculator(self.sitecol, surf)
            d_km = calc.calculate_site_to_trace_distances()  # np.ndarray[km]
            if np.any(np.asarray(d_km) <= near_thresh_km):
                return "case2"

        return "case3"


def choose_combinations(case: str) -> List[str]:
    """
    Map a case label to the active combinations.

    case1 -> ['A','B','C']
    case2 -> ['A','B']
    case3 -> ['A']
    """
    table = {"case1": ["A", "B", "C"], "case2": ["A", "B"], "case3": ["A"]}
    c = str(case).lower()
    if c not in table:
        raise ValueError(f"Unknown case: {case}")
    return table[c]


def combine_probabilities(prob_list: List[np.ndarray]) -> np.ndarray:
    """
    Combine multiple probability matrices element-wise:

        P_total = 1 - Π_i (1 - P_i)

    Each matrix must be shape (n_sites, n_displ).
    """
    if not prob_list:
        raise ValueError("No probability matrices provided.")
    if len(prob_list) == 1:
        return prob_list[0]
    one_minus = [1.0 - P for P in prob_list]
    return 1.0 - np.prod(one_minus, axis=0)


def decide_combinations_for_site(
    site: Any = None,
    sitecol: Optional[SiteCollection] = None,
    rank1p5_surfaces: Optional[Any] = None,
    case: Optional[str] = None,
    near_radius_m: float = 1000.0,
    r_distance_m: Optional[float] = None,  # accepted for interface parity (unused)
    style: Optional[str] = None,           # accepted and ignored
    **kwargs,
) -> List[str]:
    """
    Decide which combinations (A, B, C) apply at a site, based on Case 1/2/3 logic.

    Parameters
    ----------
    site : optional
        Deprecated parameter, use 'sitecol' instead.
    sitecol : SiteCollection, optional
        For a point site (or multiple points).
    rank1p5_surfaces : list or dict[name->surface]
        Rank 1.5 surfaces (objects must support get_fault_trace() and get_min_distance()).
    case: 
    near_radius_m : float
        Distance threshold (default 1000 m).
    r_distance_m : float, optional
        Accepted for API compatibility with scripts; not used here.
    style : str, optional
        Accepted and ignored.

    Returns
    -------
    list[str]
        A list like ['A'], ['A','B'], or ['A','B','C'].
    """
    # Normalize rank1p5 surfaces to a list
    if isinstance(rank1p5_surfaces, dict):
        r15_list = list(rank1p5_surfaces.values())
    else:
        r15_list = rank1p5_surfaces or []

    if case is not None:
        return choose_combinations(case)

    if sitecol is None:
        raise ValueError("'sitecol' is required.")

    clf = SiteFaultClassifier(
        sitecol=sitecol,
        rank1p5_surfaces=r15_list,
        near_radius_m=near_radius_m,
    )
    case = clf.classify_site()
    return choose_combinations(case)
