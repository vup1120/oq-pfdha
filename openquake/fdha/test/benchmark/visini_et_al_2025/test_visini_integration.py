# -*- coding: utf-8 -*-
"""
Benchmark test: Visini et al. (2025) integration tests

Integration tests for Visini model components and decision tree logic.
"""

import numpy as np
import pytest

pytestmark = pytest.mark.benchmark

# Models under test
from openquake.fdha.secondary_surf_displ.visini2025 import Visini2025SecondaryFD
from openquake.fdha.secondary_surf_rup.visini2025 import Visini2025SecondarySR

# Decision tree
from openquake.fdha.calc.decision_tree import (
    decide_combinations_for_site,
    choose_combinations,
    SiteFaultClassifier
)

# For Monte Carlo reproducibility:
RNG_SEED = 42


def test_logistic_slice_monotonicity_normal_hw():
    """Test that secondary SR probability increases with magnitude and decreases with distance."""
    model = Visini2025SecondarySR()
    style = 'normal'
    comb = 'A'
    pixel = 100  # 100 m across-strike slice
    rx = np.array([+100.0])  # HW (rx>=0 => fw=0)
    # r in meters (paper Eq. 2 uses meters)
    r_near = np.array([500.0])
    r_far = np.array([3000.0])

    p_m55_near = model.get_prob_slice(mag=5.5, r=r_near, rx=rx, style=style, pixel_size=pixel, combination=comb)
    p_m75_near = model.get_prob_slice(mag=7.5, r=r_near, rx=rx, style=style, pixel_size=pixel, combination=comb)
    p_m65_far = model.get_prob_slice(mag=6.5, r=r_far, rx=rx, style=style, pixel_size=pixel, combination=comb)

    # Higher Mw → higher P; larger r → lower P
    assert p_m75_near > p_m55_near
    assert p_m65_far < model.get_prob_slice(mag=6.5, r=r_near, rx=rx, style=style, pixel_size=pixel, combination=comb)

    # FW indicator lowers P vs HW for same Mw and r
    rx_fw = np.array([-100.0])  # FW
    p_hw = model.get_prob_slice(mag=6.5, r=r_near, rx=rx, style=style, pixel_size=pixel, combination=comb)
    p_fw = model.get_prob_slice(mag=6.5, r=r_near, rx=rx_fw, style=style, pixel_size=pixel, combination=comb)
    assert p_fw < p_hw


def test_along_strike_monte_carlo_trends_normal_hw():
    """Test that wider sites along-strike have higher probability."""
    model = Visini2025SecondarySR()
    np.random.seed(RNG_SEED)

    # Near vs far & width scaling
    fault_length = 40_000.0  # m (40 km PF)
    mechanism = 'normal'
    hw_fw = 'HW'
    near_or_far = 'near'

    # Wider site along-strike → higher probability
    # Note: optimized visini2025.py version doesn't have site_distance parameter
    p50 = model.monte_carlo_rank2_occurrence(fault_length, 50, hw_fw, mechanism, near_or_far, num_simulations=2000)
    p100 = model.monte_carlo_rank2_occurrence(fault_length, 100, hw_fw, mechanism, near_or_far, num_simulations=2000)
    p500 = model.monte_carlo_rank2_occurrence(fault_length, 500, hw_fw, mechanism, near_or_far, num_simulations=2000)
    assert 0 < p50 < p100 < p500 < 1

    # Far-fault should reduce probability relative to near-fault (F-ratio table)
    p100_far = model.monte_carlo_rank2_occurrence(fault_length, 100, hw_fw, mechanism, 'far', num_simulations=2000)
    assert p100_far < p100


def test_tpfm_forward_example():
    """Test TPFM computation from scaling relations."""
    fd = Visini2025SecondaryFD()
    tpfm, _ = fd.compute_tpfm_from_scaling(
        mag=7.0,
        style='normal',
        model='WC1994',
        norm_pos=0.5,
        dip=60.0,
        distance=0.0,  # distance parameter instead of profile
    )
    assert np.isfinite(tpfm)
    # Basic path returns AD*sin(dip); ensure positive and reasonable magnitude
    assert tpfm > 0


def test_choose_combinations():
    """Test that combination selection logic works correctly."""
    # Case 1: all combinations
    assert set(choose_combinations('case1')) == set(['A', 'B', 'C'])
    
    # Case 2: A and B only
    assert set(choose_combinations('case2')) == set(['A', 'B'])
    
    # Case 3: A only
    assert set(choose_combinations('case3')) == set(['A'])


def test_decision_tree_without_rank1p5():
    """Test decision tree for case 3 (no rank 1.5 traces nearby)."""
    from openquake.hazardlib.site import Site, SiteCollection
    from openquake.hazardlib.geo import Point
    
    # Create a simple site collection
    site = Site(Point(16.005, 39.005, 0.0))
    sitecol = SiteCollection([site])
    
    # No rank 1.5 surfaces → case 3 → combination A only
    combos = decide_combinations_for_site(
        sitecol=sitecol,
        rank1p5_surfaces=[],
        near_radius_m=1000.0
    )
    assert set(combos) == set(['A'])


def test_decision_tree_with_explicit_case():
    """Test decision tree with explicitly specified case."""
    from openquake.hazardlib.site import Site, SiteCollection
    from openquake.hazardlib.geo import Point
    
    site = Site(Point(16.005, 39.005, 0.0))
    sitecol = SiteCollection([site])
    
    # Explicit case 1 → all combinations
    combos = decide_combinations_for_site(
        sitecol=sitecol,
        case='case1'
    )
    assert set(combos) == set(['A', 'B', 'C'])
    
    # Explicit case 2 → A and B
    combos = decide_combinations_for_site(
        sitecol=sitecol,
        case='case2'
    )
    assert set(combos) == set(['A', 'B'])
    
    # Explicit case 3 → A only
    combos = decide_combinations_for_site(
        sitecol=sitecol,
        case='case3'
    )
    assert set(combos) == set(['A'])
