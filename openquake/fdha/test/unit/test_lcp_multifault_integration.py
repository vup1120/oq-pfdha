# -*- coding: utf-8 -*-
"""
Integration: the ``reference_line_method`` switch routes multi-section x/L
through the LCP reference line ('lcp') instead of the ECS default ('ecs').

Uses the same lightweight surface/site stubs as
test_ecs_multifault_integration.py.
"""
import numpy as np
import pytest

from openquake.fdha.calc.utils import rupture_distance as rd
from openquake.fdha.calc.utils import ecs, lcp

pytestmark = pytest.mark.unit


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
    s1 = _section(np.full(8, 116.47), np.linspace(-31.13, -31.09, 8))
    s2 = _section(np.full(8, 116.472), np.linspace(-31.085, -31.05, 8))
    return _MultiSurface([s1, s2])


def _sites():
    return [_Site(116.46, -31.12), _Site(116.475, -31.10), _Site(116.48, -31.06)]


def test_multisection_routes_xl_through_lcp():
    calc = rd.VectorizedRuptureDistanceCalculator(
        _sites(), _two_section_fault(), reference_line_method="lcp")
    assert isinstance(calc._refline, lcp.LcpResult)
    xl, L_km = calc.calculate_x_l_ratios()
    assert len(xl) == 3
    assert np.all((xl >= 0.0) & (xl <= 1.0))
    # both builders agree on the rupture scale (~8.9 km incl. stepover)
    assert 8.0 < L_km < 11.0
    # southern site maps near the start, northern near the end (x/L is
    # oriented along the mean rupture strike, like the ECS)
    assert xl[0] < xl[2]
    dists = calc.calculate_site_to_trace_distances()
    assert len(dists) == 3 and np.all(dists >= 0.0)


def test_default_method_is_still_ecs():
    calc = rd.VectorizedRuptureDistanceCalculator(_sites(), _two_section_fault())
    assert isinstance(calc._refline, ecs.EcsResult)
    assert calc._ecs is calc._refline      # compatibility alias


def test_lcp_and_ecs_xl_are_consistent():
    # the two reference lines are built differently but describe the same
    # rupture: x/L for well-inside sites should agree to first order
    surf = _two_section_fault()
    sites = _sites()
    xl_e, L_e = rd.VectorizedRuptureDistanceCalculator(
        sites, surf, reference_line_method="ecs").calculate_x_l_ratios()
    xl_l, L_l = rd.VectorizedRuptureDistanceCalculator(
        sites, surf, reference_line_method="lcp").calculate_x_l_ratios()
    assert L_l == pytest.approx(L_e, rel=0.2)
    np.testing.assert_allclose(xl_l, xl_e, atol=0.1)


def test_invalid_method_rejected():
    with pytest.raises(ValueError, match="reference_line_method"):
        rd.VectorizedRuptureDistanceCalculator(
            _sites(), _two_section_fault(), reference_line_method="spline")


def test_contextmaker_builds_declared_reference_lines():
    """Per-model routing: the context maker receives the union of the
    configured models' MULTIFAULT_REFERENCE_LINE declarations (collected by
    the calculator, hazardlib REQUIRES_DISTANCES-style) and builds one
    distance calculator per required method."""
    from openquake.fdha.calc.contexts import FDHAContextMaker

    class _SiteCol(list):
        lons = np.array([116.46, 116.475, 116.48])
        lats = np.array([-31.12, -31.10, -31.06])
        vs30 = np.full(3, 760.0)
        sids = np.arange(3, dtype=np.uint32)

    sitecol = _SiteCol(_sites())
    cmaker = FDHAContextMaker(
        sitecol, {'multifault_reference_lines': ('lcp',)},
        maximum_distance=50.0)
    assert cmaker.multifault_reference_lines == ('lcp',)

    surf = _two_section_fault()
    calc = cmaker._get_distance_calculator(surf, 'lcp')
    assert isinstance(calc._refline, lcp.LcpResult)
    # per-method caching: same surface, different method -> different calc
    calc_ecs = cmaker._get_distance_calculator(surf, 'ecs')
    assert isinstance(calc_ecs._refline, ecs.EcsResult)
    assert cmaker._get_distance_calculator(surf, 'lcp') is calc
