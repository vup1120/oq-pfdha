# -*- coding: utf-8 -*-
"""
Unit tests: Horizontal trace distance (r) and consistency with x/L.

Verifies that ``calculate_site_to_trace_distance(s)`` returns the
horizontal distance from the site to the surface trace polyline (km),
computed via the local equirectangular projection - NOT the OQ
``get_min_distance`` (Rrup) which depends on the 3-D mesh density.

All expected values are hand-calculable from simple plane geometry
(equirectangular projection at small scales) or from known great-circle
distances.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from openquake.fdha.calc.utils.rupture_distance import (
    _min_distance_point_to_polyline_xy,
    _horizontal_distance_to_trace_km,
    RuptureDistanceCalculator,
    VectorizedRuptureDistanceCalculator,
    project_point_onto_trace_km,
)

pytestmark = pytest.mark.unit

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _make_dummy_surface_from_trace(trace_lonlat):
    """Create a minimal surface stub whose mesh top row matches the trace."""
    class Mesh:
        def __init__(self, trace):
            self.lons = np.array([trace[:, 0]])
            self.lats = np.array([trace[:, 1]])
            self.depths = np.zeros_like(self.lons)

    class Surf:
        def __init__(self, trace):
            self.mesh = Mesh(trace)

    return Surf(np.asarray(trace_lonlat, dtype=float))


def _make_site_collection(lon, lat):
    from openquake.hazardlib.site import Site, SiteCollection
    from openquake.hazardlib.geo import Point
    return SiteCollection([Site(Point(float(lon), float(lat), 0.0))])


def _make_multi_site_collection(lonlats):
    from openquake.hazardlib.site import Site, SiteCollection
    from openquake.hazardlib.geo import Point
    return SiteCollection([Site(Point(float(lon), float(lat), 0.0)) for lon, lat in lonlats])


# ---------------------------------------------------------------------------
# Low-level: _min_distance_point_to_polyline_xy
# ---------------------------------------------------------------------------
class TestMinDistancePointToPolylineXY:
    """Tests for the XY-plane polyline distance helper."""

    def test_perpendicular_to_horizontal_segment(self):
        """Point 3 km above the midpoint of a 10 km horizontal segment."""
        poly = np.array([[0.0, 0.0], [10.0, 0.0]])
        pxy = np.array([5.0, 3.0])
        assert_allclose(_min_distance_point_to_polyline_xy(pxy, poly), 3.0, atol=1e-12)

    def test_beyond_segment_start(self):
        """Point off the start end - closest point is the start vertex."""
        poly = np.array([[0.0, 0.0], [10.0, 0.0]])
        pxy = np.array([-3.0, 4.0])  # distance = sqrt(9+16) = 5
        assert_allclose(_min_distance_point_to_polyline_xy(pxy, poly), 5.0, atol=1e-12)

    def test_beyond_segment_end(self):
        """Point off the end - closest point is the end vertex."""
        poly = np.array([[0.0, 0.0], [10.0, 0.0]])
        pxy = np.array([13.0, 4.0])  # distance = sqrt(9+16) = 5
        assert_allclose(_min_distance_point_to_polyline_xy(pxy, poly), 5.0, atol=1e-12)

    def test_l_shaped_polyline(self):
        """Point inside the L-bend - closest to the second segment."""
        poly = np.array([[0.0, 0.0], [5.0, 0.0], [5.0, 5.0]])
        pxy = np.array([6.0, 3.0])  # 1 km from second segment
        assert_allclose(_min_distance_point_to_polyline_xy(pxy, poly), 1.0, atol=1e-12)

    def test_point_on_segment(self):
        """Point exactly on the polyline."""
        poly = np.array([[0.0, 0.0], [10.0, 0.0]])
        pxy = np.array([5.0, 0.0])
        assert_allclose(_min_distance_point_to_polyline_xy(pxy, poly), 0.0, atol=1e-12)


# ---------------------------------------------------------------------------
# Mid-level: _horizontal_distance_to_trace_km (single site, geographic coords)
# ---------------------------------------------------------------------------
class TestHorizontalDistanceToTraceKm:
    """Tests for the geographic-coord wrapper (single point to trace)."""

    def test_straight_trace_on_equator_perpendicular(self):
        """Site 1 degree north of an equatorial trace ≈ 111.195 km away."""
        trace = np.array([[0.0, 0.0], [2.0, 0.0]])
        site = np.array([1.0, 1.0])  # 1 deg N of midpoint
        d = _horizontal_distance_to_trace_km(site, trace)
        # equirectangular: dy = R * dlat_rad
        expected_km = 6371.0088 * np.deg2rad(1.0)
        assert_allclose(d, expected_km, rtol=1e-4)

    def test_site_on_trace_zero_distance(self):
        """Site on the trace itself."""
        trace = np.array([[0.0, 0.0], [2.0, 0.0]])
        site = np.array([1.0, 0.0])
        d = _horizontal_distance_to_trace_km(site, trace)
        assert_allclose(d, 0.0, atol=1e-6)

    def test_degenerate_trace_returns_zero(self):
        trace = np.array([[5.0, 5.0]])
        site = np.array([5.0, 6.0])
        assert _horizontal_distance_to_trace_km(site, trace) == 0.0


# ---------------------------------------------------------------------------
# High-level: Calculator classes (scalar & vectorized)
# ---------------------------------------------------------------------------
class TestCalculatorTraceDistance:
    """Ensure calculator classes return horizontal trace distance, not Rrup."""

    def test_scalar_perpendicular_distance(self):
        """Scalar calculator: site 1 deg N of equatorial trace."""
        trace = np.array([[0.0, 0.0], [2.0, 0.0]])
        surf = _make_dummy_surface_from_trace(trace)
        sc = _make_site_collection(1.0, 1.0)
        calc = RuptureDistanceCalculator(sc, surf)
        d = calc.calculate_site_to_trace_distance()
        expected_km = 6371.0088 * np.deg2rad(1.0)
        assert_allclose(d, expected_km, rtol=1e-4)

    def test_vectorized_perpendicular_distances(self):
        """Vectorized calculator: three sites at different offsets."""
        trace = np.array([[0.0, 0.0], [2.0, 0.0]])
        surf = _make_dummy_surface_from_trace(trace)
        sites = [(1.0, 0.5), (1.0, 1.0), (1.0, 2.0)]
        sc = _make_multi_site_collection(sites)
        calc = VectorizedRuptureDistanceCalculator(sc, surf)
        dists = calc.calculate_site_to_trace_distances()
        # Geodesic truth (hazardlib EARTH_RADIUS). The local frame is
        # hazardlib's OrthographicProjection, whose planar distance is
        # R*sin(delta) vs the geodesic R*delta - a relative shortfall of
        # ~delta^2/6 (2e-4 at the 2 deg / 222 km site here; irrelevant at
        # FDHA's real near-fault ranges), hence rtol=5e-4.
        R = 6371.0
        expected = np.array([
            R * np.deg2rad(0.5),
            R * np.deg2rad(1.0),
            R * np.deg2rad(2.0),
        ])
        assert_allclose(dists, expected, rtol=5e-4)

    def test_scalar_vs_vectorized_distance_consistency(self):
        """Scalar and vectorized distance calculations must agree."""
        trace = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]])
        surf = _make_dummy_surface_from_trace(trace)
        sites_ll = [(0.5, 0.3), (1.2, 0.5), (0.0, -0.2)]
        sc_all = _make_multi_site_collection(sites_ll)
        vcalc = VectorizedRuptureDistanceCalculator(sc_all, surf)
        vec_dists = vcalc.calculate_site_to_trace_distances()

        for i, (lon, lat) in enumerate(sites_ll):
            sc = _make_site_collection(lon, lat)
            scalc = RuptureDistanceCalculator(sc, surf)
            assert_allclose(scalc.calculate_site_to_trace_distance(), vec_dists[i], rtol=1e-9)


# ---------------------------------------------------------------------------
# Consistency: r (distance) and x/L from the SAME projection
# ---------------------------------------------------------------------------
class TestDistanceXLConsistency:
    """Verify r and x/L share the same geometric projection."""

    def test_site_on_trace_has_zero_r_and_known_xL(self):
        """Site sitting exactly on the trace midpoint: r=0, x/L=0.5."""
        trace = np.array([[0.0, 0.0], [2.0, 0.0]])
        surf = _make_dummy_surface_from_trace(trace)
        sc = _make_site_collection(1.0, 0.0)
        calc = RuptureDistanceCalculator(sc, surf)
        d = calc.calculate_site_to_trace_distance()
        xL, L = calc.calculate_x_l_ratio()
        assert_allclose(d, 0.0, atol=1e-6)
        assert_allclose(xL, 0.5, rtol=1e-3)

    def test_vectorized_r_and_xL_share_L(self):
        """Vectorized r and x/L must use the same trace length."""
        trace = np.array([[0.0, 0.0], [2.0, 0.0]])
        surf = _make_dummy_surface_from_trace(trace)
        sites = [(0.5, 0.5), (1.5, -0.5)]
        sc = _make_multi_site_collection(sites)
        vcalc = VectorizedRuptureDistanceCalculator(sc, surf)
        dists = vcalc.calculate_site_to_trace_distances()
        xL_arr, L_km = vcalc.calculate_x_l_ratios()

        # Both sites 0.5 deg from trace => same r
        assert_allclose(dists[0], dists[1], rtol=1e-3)
        # x/L should be 0.25 and 0.75
        assert_allclose(xL_arr, [0.25, 0.75], rtol=1e-3)
        # L_km ≈ 222.39
        assert_allclose(L_km, 2.0 * 111.195, rtol=5e-4)

    def test_r_independent_of_mock_get_min_distance(self):
        """The new r does NOT call get_min_distance; verify by poisoning it."""
        trace = np.array([[0.0, 0.0], [2.0, 0.0]])

        class PoisonSurf:
            class mesh:
                lons = np.array([[0.0, 2.0]])
                lats = np.array([[0.0, 0.0]])
                depths = np.array([[0.0, 0.0]])

            def get_min_distance(self, sitecol):
                raise RuntimeError("Should not be called for trace distance")

        sc = _make_site_collection(1.0, 1.0)
        calc = RuptureDistanceCalculator(sc, PoisonSurf())
        d = calc.calculate_site_to_trace_distance()
        expected = 6371.0088 * np.deg2rad(1.0)
        assert_allclose(d, expected, rtol=1e-4)


# ---------------------------------------------------------------------------
# Real SimpleFaultSurface (small mesh) - optional, requires hazardlib
# ---------------------------------------------------------------------------
class TestWithSimpleFaultSurface:
    """Use a real OQ SimpleFaultSurface to verify trace distance vs manual calc."""

    def test_short_fault_coarse_mesh(self):
        """Short 2-deg fault, 2 km mesh: r must match manual trace distance."""
        from openquake.hazardlib.geo import Line, Point
        from openquake.hazardlib.geo.surface.simple_fault import SimpleFaultSurface

        trace = Line([Point(0.0, 0.0), Point(1.0, 0.0)])
        surface = SimpleFaultSurface.from_fault_data(
            fault_trace=trace,
            upper_seismogenic_depth=0.0,
            lower_seismogenic_depth=15.0,
            dip=45.0,
            mesh_spacing=2.0,
        )

        site_lonlat = np.array([0.5, 0.5])  # 0.5 deg N of trace midpoint
        expected_km = 6371.0088 * np.deg2rad(0.5)

        # Single-site
        sc1 = _make_site_collection(0.5, 0.5)
        calc1 = RuptureDistanceCalculator(sc1, surface)
        d1 = calc1.calculate_site_to_trace_distance()
        assert_allclose(d1, expected_km, rtol=5e-3)

        # Vectorized (same site + one on trace)
        sc2 = _make_multi_site_collection([(0.5, 0.5), (0.5, 0.0)])
        vcalc = VectorizedRuptureDistanceCalculator(sc2, surface)
        dists = vcalc.calculate_site_to_trace_distances()
        assert_allclose(dists[0], expected_km, rtol=5e-3)
        assert_allclose(dists[1], 0.0, atol=0.5)  # on-trace site

    def test_mesh_spacing_does_not_affect_trace_distance(self):
        """Coarse vs fine mesh must yield the same trace distance."""
        from openquake.hazardlib.geo import Line, Point
        from openquake.hazardlib.geo.surface.simple_fault import SimpleFaultSurface

        trace = Line([Point(0.0, 0.0), Point(1.0, 0.0)])
        kw = dict(
            fault_trace=trace,
            upper_seismogenic_depth=0.0,
            lower_seismogenic_depth=20.0,
            dip=30.0,
        )
        surf_coarse = SimpleFaultSurface.from_fault_data(mesh_spacing=5.0, **kw)
        surf_fine = SimpleFaultSurface.from_fault_data(mesh_spacing=1.0, **kw)

        sc = _make_multi_site_collection([(0.5, 0.3), (0.2, -0.1)])
        d_coarse = VectorizedRuptureDistanceCalculator(sc, surf_coarse).calculate_site_to_trace_distances()
        d_fine = VectorizedRuptureDistanceCalculator(sc, surf_fine).calculate_site_to_trace_distances()
        # Trace distance must be nearly identical regardless of mesh spacing.
        # Small differences arise only from the trace polyline vertex count.
        assert_allclose(d_coarse, d_fine, rtol=0.02)


# ---------------------------------------------------------------------------
# Fine-mesh trace extraction regression tests
# ---------------------------------------------------------------------------
class TestFineMeshTraceExtraction:
    """Verify _extract_fault_trace_from_mesh does not collapse fine-mesh traces.

    At rupture_mesh_spacing = 0.02 km and lon ~120°, adjacent mesh nodes
    differ by ~2e-4 degrees.  The old np.allclose (rtol=1e-5) had an
    effective threshold of 1e-5*120 = 1.2e-3 degrees, which treated all
    consecutive nodes as duplicates and collapsed the trace to 1 point.
    """

    def test_fine_spacing_preserves_all_vertices(self):
        """Nodes spaced ~1e-4 degrees apart must NOT be collapsed."""
        from openquake.fdha.calc.utils.rupture_distance import _extract_fault_trace_from_mesh

        n = 100
        lons = np.linspace(120.7, 120.7 + (n - 1) * 1.97e-4, n)
        lats = np.linspace(24.0, 24.0 + (n - 1) * 1.80e-4, n)

        class Mesh:
            pass
        mesh = Mesh()
        mesh.lons = lons.reshape(1, -1)
        mesh.lats = lats.reshape(1, -1)
        mesh.depths = np.zeros_like(mesh.lons)

        class Surf:
            pass
        surf = Surf()
        surf.mesh = mesh

        trace = _extract_fault_trace_from_mesh(surf)
        assert trace.shape[0] == n, (
            f"Expected {n} trace vertices but got {trace.shape[0]} "
            f"(fine-mesh nodes were incorrectly collapsed)"
        )

    def test_exact_duplicates_are_removed(self):
        """Genuine bitwise-identical consecutive nodes must be removed."""
        from openquake.fdha.calc.utils.rupture_distance import _extract_fault_trace_from_mesh

        lons = np.array([[120.7, 120.7, 120.8, 120.8, 120.9]])
        lats = np.array([[24.0, 24.0, 24.1, 24.1, 24.2]])

        class Mesh:
            pass
        mesh = Mesh()
        mesh.lons = lons
        mesh.lats = lats
        mesh.depths = np.zeros_like(lons)

        class Surf:
            pass
        surf = Surf()
        surf.mesh = mesh

        trace = _extract_fault_trace_from_mesh(surf)
        # 5 nodes → 3 unique consecutive nodes
        assert trace.shape[0] == 3
        assert_allclose(trace[:, 0], [120.7, 120.8, 120.9])

    def test_all_identical_raises(self):
        """Surface with all-identical nodes must raise, not silently return r=0."""
        from openquake.fdha.calc.utils.rupture_distance import _extract_fault_trace_from_mesh

        lons = np.array([[120.7, 120.7, 120.7]])
        lats = np.array([[24.0, 24.0, 24.0]])

        class Mesh:
            pass
        mesh = Mesh()
        mesh.lons = lons
        mesh.lats = lats
        mesh.depths = np.zeros_like(lons)

        class Surf:
            pass
        surf = Surf()
        surf.mesh = mesh

        with pytest.raises(ValueError, match="only 1 distinct point"):
            _extract_fault_trace_from_mesh(surf)

    def test_real_surface_fine_mesh(self):
        """Real SimpleFaultSurface with fine mesh preserves many trace vertices."""
        from openquake.hazardlib.geo import Line, Point
        from openquake.hazardlib.geo.surface.simple_fault import SimpleFaultSurface
        from openquake.fdha.calc.utils.rupture_distance import _extract_fault_trace_from_mesh

        trace_line = Line([Point(120.7, 24.0), Point(120.75, 24.05)])
        surface = SimpleFaultSurface.from_fault_data(
            fault_trace=trace_line,
            upper_seismogenic_depth=0.0,
            lower_seismogenic_depth=12.0,
            dip=30.0,
            mesh_spacing=0.05,
        )
        trace = _extract_fault_trace_from_mesh(surface)
        # ~7 km trace at 0.05 km spacing → ~140 vertices
        assert trace.shape[0] >= 50, (
            f"Expected >= 50 trace vertices for 0.05 km mesh but got {trace.shape[0]}"
        )
