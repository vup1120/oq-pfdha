"""
Unit tests: Rupture distance x/L boundary conditions

Tests boundary conditions and edge cases for x/L ratio calculations.
"""

import numpy as np
import math
import pytest

from openquake.fdha.calc.utils.rupture_distance import (
    project_point_onto_trace_km,
    VectorizedRuptureDistanceCalculator,
)

pytestmark = pytest.mark.unit


class _Mesh:
    def __init__(self, lons, lats, depths):
        self.lons = np.asarray(lons)
        self.lats = np.asarray(lats)
        self.depths = np.asarray(depths)


class _DummySurface:
    def __init__(self, trace):
        # trace: array of shape (N, 2) [lon, lat], assume surface at z=0
        trace = np.asarray(trace, dtype=float)
        # create a 1-row mesh along the trace
        lons = trace[:, 0][np.newaxis, :]
        lats = trace[:, 1][np.newaxis, :]
        deps = np.zeros_like(lons)
        self.mesh = _Mesh(lons, lats, deps)


class _Loc:
    def __init__(self, lon, lat):
        self.longitude = float(lon)
        self.latitude = float(lat)


class _Site:
    def __init__(self, lon, lat):
        self.location = _Loc(lon, lat)


def _ratio(site_lon, site_lat, trace):
    x_km, L_km = project_point_onto_trace_km(np.array([site_lon, site_lat]), np.asarray(trace))
    return 0.0 if L_km <= 0.0 else x_km / L_km


def test_project_point_horizontal_line_endpoints():
    # Trace: (0,0) -> (1,0)
    trace = np.array([[0.0, 0.0], [1.0, 0.0]])
    # Exactly at start
    r0 = _ratio(0.0, 0.0, trace)
    assert math.isfinite(r0)
    assert abs(r0 - 0.0) <= 1e-10
    # Exactly at end
    r1 = _ratio(1.0, 0.0, trace)
    assert math.isfinite(r1)
    assert abs(r1 - 1.0) <= 1e-10
    # Middle
    rm = _ratio(0.5, 0.0, trace)
    assert math.isfinite(rm)
    assert abs(rm - 0.5) <= 1e-10


def test_project_point_near_end_robust_to_eps():
    trace = np.array([[0.0, 0.0], [1.0, 0.0]])
    # Slightly off the end in latitude; orthogonal projection parameter clamps to 1
    r = _ratio(1.0, 1e-9, trace)
    assert math.isfinite(r)
    assert abs(r - 1.0) <= 1e-10


def test_vectorized_ratios_two_sites_start_and_end():
    trace = np.array([[0.0, 0.0], [1.0, 0.0]])
    surf = _DummySurface(trace)
    # Two sites: start and end
    sitecol = [_Site(0.0, 0.0), _Site(1.0, 0.0)]
    calc = VectorizedRuptureDistanceCalculator(sitecol, surf)
    ratios, L_km = calc.calculate_x_l_ratios()
    assert ratios.shape == (2,)
    assert L_km > 0
    assert abs(ratios[0] - 0.0) <= 1e-10
    assert abs(ratios[1] - 1.0) <= 1e-10


