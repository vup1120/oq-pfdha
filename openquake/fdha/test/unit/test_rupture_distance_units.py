"""
Unit tests: Rupture distance calculation units

Tests the rupture distance calculation utilities, including:
- Point projection onto trace
- Distance calculations in km
- Vectorized operations
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from openquake.fdha.calc.utils.rupture_distance import project_point_onto_trace_km

pytestmark = pytest.mark.unit


def _make_dummy_surface_from_trace(trace_lonlat):
    class Mesh:
        def __init__(self, trace):
            # Minimal mesh with a single row at depth 0 following the trace
            self.lons = np.array([trace[:, 0]])
            self.lats = np.array([trace[:, 1]])
            self.depths = np.zeros_like(self.lons)

    class Surf:
        def __init__(self, trace):
            self.mesh = Mesh(trace)

        def get_min_distance(self, sitecol):
            # Return a dummy fixed small distance in km per site; not used in these tests
            return np.array([0.1] * len(sitecol))

    return Surf(np.asarray(trace_lonlat, dtype=float))


def _make_site_collection(lon, lat):
    from openquake.hazardlib.site import Site, SiteCollection
    from openquake.hazardlib.geo import Point
    return SiteCollection([Site(Point(float(lon), float(lat), 0.0))])


@pytest.mark.parametrize(
    "trace, sites, expected_x_over_L, expected_L_km, rtol",
    [
        # Straight line on equator: (0,0)->(2,0); L ≈ 222.39 km
        (
            np.array([[0.0, 0.0], [2.0, 0.0]]),
            [(0.5, 0.0), (1.0, 0.0), (2.0, 0.0)],
            [0.25, 0.5, 1.0],
            2.0 * 111.195,
            5e-4,
        ),
        # L-shape: (0,0)->(1,0)->(1,1); L ~ 111.195 + 111.195 ≈ 222.39 km
        (
            np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]),
            [(1.0, 0.5), (1.5, 0.2)],
            [0.75, 0.60],
            2.0 * 111.195,
            1e-2,
        ),
        # High latitude 60°: (-1,60)->(1,60); L ≈ 2° * 111.195 * cos(60°)
        (
            np.array([[-1.0, 60.0], [1.0, 60.0]]),
            [(0.0, 60.0)],
            [0.5],
            2.0 * 111.195 * np.cos(np.deg2rad(60.0)),
            5e-4,
        ),
        # IDL crossing: (179,0)->(-179,0); L ≈ 2° * 111.195
        (
            np.array([[179.0, 0.0], [-179.0, 0.0]]),
            [(180.0, 0.0)],
            [0.5],
            2.0 * 111.195,
            5e-4,
        ),
    ],
)
def test_project_point_onto_trace_km_cases(trace, sites, expected_x_over_L, expected_L_km, rtol):
    for (slon, slat), xL_exp in zip(sites, expected_x_over_L):
        x_km, L_km = project_point_onto_trace_km(np.array([slon, slat]), trace)
        x_over_L = 0.0 if L_km <= 0 else x_km / L_km
        assert_allclose(x_over_L, xL_exp, rtol=rtol, atol=0.0)
        assert_allclose(L_km, expected_L_km, rtol=rtol, atol=0.05)


def test_degenerate_trace_returns_zero():
    trace = np.array([[10.0, 0.0]])  # single point
    x_km, L_km = project_point_onto_trace_km(np.array([10.0, 0.0]), trace)
    assert x_km == 0.0 and L_km == 0.0

