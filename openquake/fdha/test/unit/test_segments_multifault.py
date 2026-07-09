# -*- coding: utf-8 -*-
"""
Segmentation-direct distances for multi-section ruptures
(``reference_line_method='segments'``): no smoothed representative line —
r is the distance to the nearest actual section trace (gaps NOT bridged)
and x/L comes from raw MultiLine GC2. This is the distance treatment
required by segmentation-calibrated distributed-FD models (Visini et al.
2025), where a smoothed ECS/LCP trace would wrongly assign r ~ 0 inside an
inter-section gap.
"""
import numpy as np
import pytest

from openquake.hazardlib.geo.line import Line
from openquake.hazardlib.geo.multiline import MultiLine

from openquake.fdha.calc.utils import rupture_distance as rd
from openquake.fdha.calc.utils import segments as seg

pytestmark = pytest.mark.unit


# --- stubs shared with the ECS/LCP integration tests ----------------------- #
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


def _gapped_fault():
    # two E-W collinear sections at lat 42.0 with a ~0.1 deg (~8.3 km) gap:
    # [13.00 .. 13.10]  GAP  [13.20 .. 13.30]
    s1 = _section(np.linspace(13.00, 13.10, 5), np.full(5, 42.0))
    s2 = _section(np.linspace(13.20, 13.30, 5), np.full(5, 42.0))
    return _MultiSurface([s1, s2])


KM_PER_DEG_LON_42 = 111.19 * np.cos(np.radians(42.0))   # ~82.6 km


# --------------------------------------------------------------------------- #
# core property: gaps are NOT bridged
# --------------------------------------------------------------------------- #
def test_gap_site_gets_true_segment_distance():
    surf = _gapped_fault()
    gap_mid = _Site(13.15, 42.0)                # centre of the gap, ON strike
    calc = rd.VectorizedRuptureDistanceCalculator(
        [gap_mid], surf, reference_line_method="segments")
    assert isinstance(calc._refline, seg.SegmentsResult)
    (r,) = calc.calculate_site_to_trace_distances()
    # true distance to the nearest section end = 0.05 deg of lon at lat 42
    expected = 0.05 * KM_PER_DEG_LON_42
    assert r == pytest.approx(expected, rel=0.02)

    # contrast: the ECS representative line passes through the gap -> r ~ 0
    calc_ecs = rd.VectorizedRuptureDistanceCalculator(
        [gap_mid], surf, reference_line_method="ecs")
    (r_ecs,) = calc_ecs.calculate_site_to_trace_distances()
    assert r_ecs < 0.2 * expected


def test_on_segment_sites_have_zero_distance():
    surf = _gapped_fault()
    sites = [_Site(13.05, 42.0), _Site(13.25, 42.0)]     # on each section
    calc = rd.VectorizedRuptureDistanceCalculator(
        [sites[0], sites[1]], surf, reference_line_method="segments")
    r = calc.calculate_site_to_trace_distances()
    assert np.all(r < 0.01)


def test_off_strike_distance_is_perpendicular():
    surf = _gapped_fault()
    site = _Site(13.05, 42.05)                  # 0.05 deg N of section 1
    calc = rd.VectorizedRuptureDistanceCalculator(
        [site], surf, reference_line_method="segments")
    (r,) = calc.calculate_site_to_trace_distances()
    assert r == pytest.approx(0.05 * 111.19, rel=0.02)


def test_buried_sections_do_not_attract_r():
    # rupture = one surface section + the top edge of a deep (8 km) section
    # 10 km to the south: surface-displacement r must ignore the buried one
    s_surf = _section(np.linspace(13.00, 13.10, 5), np.full(5, 42.0), depth=0.0)
    s_deep = _section(np.linspace(13.00, 13.10, 5), np.full(5, 41.91), depth=8.2)
    surf = _MultiSurface([s_surf, s_deep])
    site = _Site(13.05, 41.93)                  # 2.2 km from the buried top,
    calc = rd.VectorizedRuptureDistanceCalculator(  # ~7.8 km from the trace
        [site], surf, reference_line_method="segments")
    (r,) = calc.calculate_site_to_trace_distances()
    assert r == pytest.approx(0.07 * 111.19, rel=0.03)   # to the SURFACE trace
    # ECS (unchanged behaviour) smooths through both edges -> much smaller r
    calc_ecs = rd.VectorizedRuptureDistanceCalculator(
        [site], surf, reference_line_method="ecs")
    (r_ecs,) = calc_ecs.calculate_site_to_trace_distances()
    assert r_ecs < r


# --------------------------------------------------------------------------- #
# x/L and L from raw MultiLine GC2 (no smoothing)
# --------------------------------------------------------------------------- #
def test_x_l_matches_raw_multiline_gc2():
    surf = _gapped_fault()
    traces = rd._section_traces(surf)
    res = seg.segments_from_traces(traces)

    lons = np.array([13.02, 13.15, 13.28])
    lats = np.array([42.01, 42.0, 41.99])
    xl, L_m = res.x_l(lons, lats)

    # independent reference computation straight on oq MultiLine
    ml = MultiLine([Line.from_vectors(np.array(lo), np.array(la))
                    for lo, la in traces])
    vlon = np.concatenate([lo for lo, _ in traces])
    vlat = np.concatenate([la for _, la in traces])
    _, u_v = ml.get_tu(vlon.copy(), vlat.copy())
    _, u_q = ml.get_tu(lons.copy(), lats.copy())
    L_ref = u_v.max() - u_v.min()
    xl_ref = np.clip((u_q - u_v.min()) / L_ref, 0, 1)

    np.testing.assert_allclose(xl, xl_ref, atol=1e-12)
    assert L_m == pytest.approx(L_ref * 1000.0)
    # L spans the whole system incl. the gap (~0.30 deg of lon at lat 42)
    assert L_m / 1000.0 == pytest.approx(0.30 * KM_PER_DEG_LON_42, rel=0.03)
    # sites near the W end / gap centre / E end map to low / mid / high x/L
    order = xl if xl[0] < xl[-1] else xl[::-1]
    assert order[0] < 0.2 and 0.3 < order[1] < 0.7 and order[2] > 0.8


def test_x_l_ratios_through_calculator():
    surf = _gapped_fault()
    sites = [_Site(13.01, 42.0), _Site(13.29, 42.0)]
    calc = rd.VectorizedRuptureDistanceCalculator(
        sites, surf, reference_line_method="segments")
    xl, L_km = calc.calculate_x_l_ratios()
    assert np.all((xl >= 0.0) & (xl <= 1.0))
    assert abs(xl[0] - xl[1]) > 0.8              # opposite ends of the system
    assert L_km == pytest.approx(0.30 * KM_PER_DEG_LON_42, rel=0.03)


# --------------------------------------------------------------------------- #
# per-model wiring (MULTIFAULT_REFERENCE_LINE declarations)
# --------------------------------------------------------------------------- #
def test_contextmaker_always_computes_segments_mask_r():
    """The canonical ctx.r on a multi-section rupture is ALWAYS the
    segments (nearest surface-reaching section) distance — the
    principal/distributed mask and the Visini near/far label must not see
    a bridged inter-section gap — even when no configured model declared
    'segments' (here the union is just the 'lcp' default)."""
    from openquake.fdha.calc.contexts import FDHAContextMaker

    class _SiteCol(list):
        lons = np.array([13.05, 13.15])
        lats = np.array([42.0, 42.0])
        vs30 = np.full(2, 760.0)
        sids = np.arange(2, dtype=np.uint32)

    sitecol = _SiteCol([_Site(13.05, 42.0), _Site(13.15, 42.0)])
    cmaker = FDHAContextMaker(
        sitecol, {'multifault_reference_lines': ('lcp',)},
        maximum_distance=50.0)
    assert cmaker.multifault_reference_lines == ('lcp',)

    calc = cmaker._get_distance_calculator(_gapped_fault(), 'segments')
    assert isinstance(calc._refline, seg.SegmentsResult)
    r = calc.calculate_site_to_trace_distances()
    # on-section site ~0; gap-centre site keeps its true segment distance
    assert r[0] < 0.01
    assert r[1] == pytest.approx(0.05 * KM_PER_DEG_LON_42, rel=0.02)


def test_visini_models_declare_segments():
    """The Visini et al. (2025) models declare the segmentation-direct
    treatment; the generic bases default to the smoothed LCP line."""
    from openquake.fdha.secondary_surf_rup.visini2025 import (
        Visini2025SecondarySR)
    from openquake.fdha.secondary_surf_displ.visini2025 import (
        Visini2025SecondaryFD)
    from openquake.fdha.secondary_surf_rup.petersen2011 import (
        Petersen2011SecondarySR)

    assert Visini2025SecondarySR.MULTIFAULT_REFERENCE_LINE == 'segments'
    assert Visini2025SecondaryFD.MULTIFAULT_REFERENCE_LINE == 'segments'
    assert Petersen2011SecondarySR.MULTIFAULT_REFERENCE_LINE == 'lcp'


def test_invalid_method_message_mentions_segments():
    with pytest.raises(ValueError, match="segments"):
        rd.VectorizedRuptureDistanceCalculator(
            [_Site(13.0, 42.0)], _gapped_fault(),
            reference_line_method="spline")


def test_single_trace_input_validation():
    with pytest.raises(ValueError, match=">=1 trace"):
        seg.SegmentsResult([(np.array([13.0]), np.array([42.0]))])
