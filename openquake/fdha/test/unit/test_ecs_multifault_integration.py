# -*- coding: utf-8 -*-
"""
Integration: multi-section ruptures route x/L through the ECS reference line.

A ``multiFaultSource`` rupture arrives as a MultiSurface whose sections each
expose a top-of-rupture line (``.tor``). ``RuptureDistanceCalculator`` detects
this and builds the ECS instead of the (invalid) single-trace extraction.
Single-strand ruptures keep the existing path. These tests use lightweight
surface/site stubs (the real emme26 multiFaultSource is 130 MB).
"""
import numpy as np
import pytest

from openquake.fdha.calc.utils import rupture_distance as rd
from openquake.fdha.calc.utils import ecs

pytestmark = pytest.mark.unit


# --- lightweight stubs mimicking the oq API the calculator touches --------- #
class _Loc:
    def __init__(self, lon, lat):
        self.longitude = float(lon)
        self.latitude = float(lat)


class _Site:
    def __init__(self, lon, lat):
        self.location = _Loc(lon, lat)


class _Tor:
    def __init__(self, coo):
        self.coo = np.asarray(coo, dtype=float)


class _Section:
    def __init__(self, coo):
        self.tor = _Tor(coo)


class _MultiSurface:
    def __init__(self, sections):
        self.surfaces = sections


def _section(lons, lats, depth=0.0):
    return _Section(np.column_stack([lons, lats, np.full(len(lons), depth)]))


def _two_section_fault():
    # two roughly colinear ~N-S segments offset along strike (a 2-section fault)
    s1 = _section(np.full(8, 116.47), np.linspace(-31.13, -31.09, 8))
    s2 = _section(np.full(8, 116.472), np.linspace(-31.085, -31.05, 8))
    return _MultiSurface([s1, s2])


# --------------------------------------------------------------------------- #
def test_section_traces_detects_multisection():
    surf = _two_section_fault()
    traces = rd._section_traces(surf)
    assert traces is not None and len(traces) == 2
    # single-section surface -> None (single-strand path)
    assert rd._section_traces(_MultiSurface([surf.surfaces[0]])) is None
    # non-multisurface object -> None
    assert rd._section_traces(object()) is None


def test_multisection_routes_xl_through_ecs():
    surf = _two_section_fault()
    sites = [_Site(116.46, -31.12), _Site(116.475, -31.10), _Site(116.48, -31.06)]
    calc = rd.VectorizedRuptureDistanceCalculator(sites, surf)
    # the ECS was built and adopted as the reference trace
    assert calc._ecs is not None
    assert isinstance(calc._ecs, ecs.EcsResult)
    xl, L_km = calc.calculate_x_l_ratios()
    assert len(xl) == len(sites)
    assert np.all(xl >= 0.0) and np.all(xl <= 1.0)
    assert L_km > 0.0
    # site-to-trace distance also works (uses the ECS nodes as the trace)
    dists = calc.calculate_site_to_trace_distances()
    assert len(dists) == len(sites) and np.all(dists >= 0.0)


def test_ecs_from_traces_direct():
    surf = _two_section_fault()
    traces = rd._section_traces(surf)
    res = ecs.ecs_from_traces(traces)
    assert isinstance(res, ecs.EcsResult)
    assert res.u.max() > res.u.min()
    xl, L = res.x_l(np.array([116.47]), np.array([-31.10]))
    assert 0.0 <= xl[0] <= 1.0 and L > 0
