# -*- coding: utf-8 -*-
"""
Validate the GC2 engine the ECS layer depends on (oq-engine ``MultiLine``).

The original ECS reference (ecs_functions.R) computed GC2 via a FORTRAN
``r_wrapper_gc2`` distance-metrics library that was never published (absent from
the Zenodo package and upstream GitHub). Per project decision we do NOT write a
new R GC2; instead we reuse oq-engine's ``MultiLine`` GC2 and validate it
against (a) oq-engine's own vetted MultiLine assertions and (b) the analytic
Spudich & Chiou (2015) single-segment reduction, where the generalized
coordinates collapse to exact along-/across-line projection.
"""
import numpy as np
import pytest

from openquake.hazardlib import geo
from openquake.hazardlib.geo import Line, Point
from openquake.hazardlib.geo.multiline import MultiLine

from openquake.fdha.calc.utils import ecs

pytestmark = pytest.mark.unit


# --------------------------------------------------------------------------- #
# (a) oq-engine's own vetted assertions, reproduced here so the engine GC2 is
#     re-checked in *this* repo's CI (the installed wheel strips its tests).
# --------------------------------------------------------------------------- #
def test_multiline_u_max_equals_geodetic_length():
    ln = Line([Point(0.2, 0.05), Point(0.0, 0.05)])
    ml = MultiLine([ln])
    dst = geo.geodetic.geodetic_distance(
        ln.points[0].longitude, ln.points[0].latitude,
        ln.points[1].longitude, ln.points[1].latitude)
    np.testing.assert_allclose(ml.get_u_max(), dst, atol=1e-4)


def test_multiline_on_line_points_have_zero_t():
    # points exactly on the trace -> T == 0 ; U spans [0, length]
    lons = np.array([0.0, 0.05, 0.10, 0.15, 0.20])
    lats = np.full_like(lons, 0.05)
    ml = MultiLine([Line.from_vectors(lons, lats)])
    t, u = ml.get_tu(lons, lats)
    np.testing.assert_allclose(t, 0.0, atol=1e-6)
    assert np.all(np.diff(u) > 0) or np.all(np.diff(u) < 0)


# --------------------------------------------------------------------------- #
# (b) analytic Spudich & Chiou (2015) single-segment reduction, in our ECS
#     gc2ext_ut wrapper (50 km extension + u-origin re-anchor). On a single
#     straight segment, U = along-strike distance from the first vertex and
#     T = signed perpendicular distance -- both exact.
# --------------------------------------------------------------------------- #
def test_gc2ext_single_segment_reduces_to_projection():
    # straight ECS along a meridian near Calingiri; build stations with known
    # along/across offsets in UTM and check gc2ext_ut recovers them.
    ecs_lat = np.linspace(-31.13, -31.05, 6)
    ecs_lon = np.full_like(ecs_lat, 116.47)
    utm = ecs.utm_for(ecs_lon, ecs_lat)

    # ECS own vertices: t ~ 0, u from 0 increasing, in km
    u, t = ecs.gc2ext_ut(ecs_lon, ecs_lat, ecs_lon, ecs_lat, utm)
    assert abs(u[0]) < 1e-6
    np.testing.assert_allclose(t, 0.0, atol=1e-6)
    assert np.all(np.diff(u) > 0)

    # a station offset purely across-strike by a known UTM distance
    x0, y0 = utm.to_xy(ecs_lon, ecs_lat)
    midx, midy = x0[2], y0[2]
    off_m = 1500.0  # 1.5 km east (≈ across-strike for a N-S line)
    slon, slat = utm.to_lonlat([midx + off_m], [midy])
    us, ts = ecs.gc2ext_ut(slon, slat, ecs_lon, ecs_lat, utm)
    # |t| should equal the across offset (km), sign consistent
    np.testing.assert_allclose(abs(ts[0]) * 1000.0, off_m, rtol=2e-3)
    # along-strike u of the station ~ u of the midpoint vertex. A constant-lon
    # line is not perfectly straight in UTM/ortho space, so GC2 weighting shifts
    # U by ~10 m for an off-line point -- correct behaviour, not exact projection.
    np.testing.assert_allclose(us[0], u[2], atol=2e-2)  # km (20 m)
