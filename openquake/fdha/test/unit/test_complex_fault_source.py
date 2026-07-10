# -*- coding: utf-8 -*-
"""
End-to-end tests for complex-fault support.

Geometry under test: a straight N-S fault at lon 0.0, lat 0.0 -> 0.3
(L = 0.3 deg x 111.195 km/deg = 33.36 km), dipping 30 deg east from 0 to
15 km depth. The down-dip horizontal extent is 15/tan(30) = 25.98 km, so
the bottom edge sits 0.23365 deg east of the top edge. The same plane is
written three ways:

- a standalone ``complexFaultSource`` (floating ruptures over the
  interpolated edge mesh; requires ``complex_fault_mesh_spacing``);
- a ``simpleFaultSource`` twin of identical trace/dip/depths and MFD, used
  as the reference in the equivalence test;
- a ``characteristicFaultSource`` carrying the same ``complexFaultGeometry``
  (single full-surface rupture with the exact NRML top edge retained as
  ``original_trace``).

The straight trace makes the FDHA site metrics analytic: a site 3 km east
of the trace midpoint must see r = rx = 3.0 km (hanging wall), x/L = 0.5,
L = 33.36 km, ztor = 0.
"""

import numpy as np
import pytest

from openquake.fdha.calc.utils.parsing import parse_source_model_faults
from openquake.fdha.calc.contexts import FDHAContextMaker
from openquake.hazardlib.site import Site, SiteCollection
from openquake.hazardlib.geo import Point

pytestmark = pytest.mark.unit


KM_PER_DEG = 111.195          # hazardlib spherical-earth degree length
TRACE_LEN_KM = 0.3 * KM_PER_DEG   # 33.36 km
DIP_DEG = 30.0
BOTTOM_LON = 0.23365          # 15/tan(30 deg) = 25.98 km east, in degrees

# Aki & Richards: points ordered S->N so the east-dipping plane lies to the
# right of the strike direction.
COMPLEX_SOURCE_XML = """<?xml version='1.0' encoding='utf-8'?>
<nrml xmlns:gml="http://www.opengis.net/gml"
      xmlns="http://openquake.org/xmlns/nrml/0.5">
  <sourceModel name="standalone complex fault">
    <sourceGroup name="g1" tectonicRegion="Active Shallow Crust">
      <complexFaultSource id="cfs1" name="Straight complex fault"
                          tectonicRegion="Active Shallow Crust">
        <complexFaultGeometry>
          <faultTopEdge>
            <gml:LineString><gml:posList>
              0.0 0.0 0.0
              0.0 0.3 0.0
            </gml:posList></gml:LineString>
          </faultTopEdge>
          <faultBottomEdge>
            <gml:LineString><gml:posList>
              0.23365 0.0 15.0
              0.23365 0.3 15.0
            </gml:posList></gml:LineString>
          </faultBottomEdge>
        </complexFaultGeometry>
        <magScaleRel>WC1994</magScaleRel>
        <ruptAspectRatio>2.0</ruptAspectRatio>
        <truncGutenbergRichterMFD aValue="4.0" bValue="1.0"
                                  minMag="6.5" maxMag="7.0"/>
        <rake>-90.0</rake>
      </complexFaultSource>
    </sourceGroup>
  </sourceModel>
</nrml>
"""

SIMPLE_TWIN_XML = """<?xml version='1.0' encoding='utf-8'?>
<nrml xmlns:gml="http://www.opengis.net/gml"
      xmlns="http://openquake.org/xmlns/nrml/0.5">
  <sourceModel name="simple twin of the complex fault">
    <sourceGroup name="g1" tectonicRegion="Active Shallow Crust">
      <simpleFaultSource id="sfs1" name="Straight simple fault"
                         tectonicRegion="Active Shallow Crust">
        <simpleFaultGeometry>
          <gml:LineString><gml:posList>
            0.0 0.0
            0.0 0.3
          </gml:posList></gml:LineString>
          <dip>30.0</dip>
          <upperSeismoDepth>0.0</upperSeismoDepth>
          <lowerSeismoDepth>15.0</lowerSeismoDepth>
        </simpleFaultGeometry>
        <magScaleRel>WC1994</magScaleRel>
        <ruptAspectRatio>2.0</ruptAspectRatio>
        <truncGutenbergRichterMFD aValue="4.0" bValue="1.0"
                                  minMag="6.5" maxMag="7.0"/>
        <rake>-90.0</rake>
      </simpleFaultSource>
    </sourceGroup>
  </sourceModel>
</nrml>
"""

CHAR_COMPLEX_XML = """<?xml version='1.0' encoding='utf-8'?>
<nrml xmlns:gml="http://www.opengis.net/gml"
      xmlns="http://openquake.org/xmlns/nrml/0.5">
  <sourceModel name="characteristic source with complex geometry">
    <sourceGroup name="g1" tectonicRegion="Active Shallow Crust">
      <characteristicFaultSource id="chc1" name="Char complex fault"
                                 tectonicRegion="Active Shallow Crust">
        <incrementalMFD binWidth="0.1" minMag="7.0">
          <occurRates>0.001</occurRates>
        </incrementalMFD>
        <rake>-90.0</rake>
        <surface>
          <complexFaultGeometry>
            <faultTopEdge>
              <gml:LineString><gml:posList>
                0.0 0.0 0.0
                0.0 0.3 0.0
              </gml:posList></gml:LineString>
            </faultTopEdge>
            <faultBottomEdge>
              <gml:LineString><gml:posList>
                0.23365 0.0 15.0
                0.23365 0.3 15.0
              </gml:posList></gml:LineString>
            </faultBottomEdge>
          </complexFaultGeometry>
        </surface>
      </characteristicFaultSource>
    </sourceGroup>
  </sourceModel>
</nrml>
"""

# 3 km east of the trace midpoint: hanging wall of the east-dipping plane.
HW_SITE = (3.0 / KM_PER_DEG, 0.15)
# 3 km west: footwall.
FW_SITE = (-3.0 / KM_PER_DEG, 0.15)


def _sitecol(*lonlats):
    return SiteCollection([
        Site(location=Point(lon, lat), vs30=760.0) for lon, lat in lonlats])


@pytest.fixture(scope="module")
def complex_source(tmp_path_factory):
    base = tmp_path_factory.mktemp("cfs")
    xml = base / "complex_source.xml"
    xml.write_text(COMPLEX_SOURCE_XML)
    srcs = parse_source_model_faults(
        [str(xml)], rupture_mesh_spacing=1.0, complex_fault_mesh_spacing=1.0,
        width_of_mfd_bin=0.1)
    return srcs


@pytest.fixture(scope="module")
def char_complex_source(tmp_path_factory):
    base = tmp_path_factory.mktemp("chc")
    xml = base / "char_complex.xml"
    xml.write_text(CHAR_COMPLEX_XML)
    srcs = parse_source_model_faults(
        [str(xml)], rupture_mesh_spacing=1.0, complex_fault_mesh_spacing=1.0,
        width_of_mfd_bin=0.1)
    return srcs


class TestParsing:
    def test_requires_complex_fault_mesh_spacing(self, tmp_path):
        """The OQ converter refuses complexFaultGeometry without the
        dedicated mesh-spacing parameter (rupture_mesh_spacing is not a
        fallback)."""
        xml = tmp_path / "complex_source.xml"
        xml.write_text(COMPLEX_SOURCE_XML)
        with pytest.raises(Exception, match="complex_fault_mesh_spacing"):
            parse_source_model_faults([str(xml)], rupture_mesh_spacing=1.0)

    def test_parse_standalone(self, complex_source):
        (src,) = complex_source.values()
        assert src.__class__.__name__ == "ComplexFaultSource"
        rups = list(src.iter_ruptures())
        assert len(rups) > 0
        # parametric floating ruptures with annual rates, GR total rate
        total = sum(r.occurrence_rate for r in rups)
        expected_total = 10 ** (4.0 - 1.0 * 6.5) - 10 ** (4.0 - 1.0 * 7.0)
        assert total == pytest.approx(expected_total, rel=1e-6)
        # the fault reaches the surface
        assert min(r.surface.get_top_edge_depth() for r in rups) == (
            pytest.approx(0.0, abs=0.01))

    def test_parse_characteristic_with_complex_geometry(
            self, char_complex_source):
        (src,) = char_complex_source.values()
        assert src.__class__.__name__ == "CharacteristicFaultSource"
        rups = list(src.iter_ruptures())
        assert len(rups) == 1
        assert rups[0].mag == pytest.approx(7.0)
        assert rups[0].occurrence_rate == pytest.approx(0.001)

    def test_original_trace_attached_from_top_edge(self, char_complex_source):
        """characteristicFaultSource + complexFaultGeometry gets the exact
        NRML faultTopEdge as original_trace (depth column dropped)."""
        (src,) = char_complex_source.values()
        trace = getattr(src.surface, "original_trace", None)
        assert trace is not None
        np.testing.assert_allclose(
            trace, [[0.0, 0.0], [0.0, 0.3]], rtol=0.0, atol=1e-12)


class TestSiteMetrics:
    """FDHA distance metrics against the analytic straight-fault geometry."""

    def test_char_complex_hanging_wall_metrics(self, char_complex_source):
        (src,) = char_complex_source.values()
        (rup,) = src.iter_ruptures()
        cm = FDHAContextMaker(_sitecol(HW_SITE), {}, maximum_distance=50.0)
        ctx = cm.get_ctx(rup)
        assert ctx is not None
        assert ctx.r[0] == pytest.approx(3.0, abs=0.05)
        assert ctx.rx[0] == pytest.approx(3.0, abs=0.05)   # hanging wall: +
        assert ctx.x_L[0] == pytest.approx(0.5, abs=0.01)
        assert ctx.L[0] == pytest.approx(TRACE_LEN_KM, abs=0.2)
        assert ctx.ztor[0] == pytest.approx(0.0, abs=0.01)

    def test_char_complex_footwall_metrics(self, char_complex_source):
        (src,) = char_complex_source.values()
        (rup,) = src.iter_ruptures()
        cm = FDHAContextMaker(_sitecol(FW_SITE), {}, maximum_distance=50.0)
        ctx = cm.get_ctx(rup)
        assert ctx is not None
        assert ctx.r[0] == pytest.approx(3.0, abs=0.05)
        assert ctx.rx[0] == pytest.approx(-3.0, abs=0.05)  # footwall: -
        assert ctx.x_L[0] == pytest.approx(0.5, abs=0.01)

    def test_standalone_full_length_rupture_metrics(self, complex_source):
        """The largest floating ruptures span the whole 33.4 km trace, so
        the mid-trace hanging-wall site sees the same analytic metrics
        (within the resampled-mesh discretisation)."""
        (src,) = complex_source.values()
        rups = [r for r in src.iter_ruptures()
                if r.surface.get_top_edge_depth() < 0.01]
        assert rups, "no surface-reaching ruptures"
        big = max(rups, key=lambda r: r.mag)
        assert big.surface.get_area() > 500.0  # near-full plane at Mw~6.95
        cm = FDHAContextMaker(_sitecol(HW_SITE), {}, maximum_distance=50.0)
        ctx = cm.get_ctx(big)
        assert ctx is not None
        assert ctx.r[0] == pytest.approx(3.0, abs=0.15)
        assert ctx.rx[0] > 0.0
        assert ctx.x_L[0] == pytest.approx(0.5, abs=0.05)
        assert ctx.L[0] == pytest.approx(TRACE_LEN_KM, abs=1.5)


class TestBuriedRuptureGate:
    """Ruptures whose top edge lies deeper than 0.5 km must be skipped:
    they cannot produce surface fault displacement. The gate is
    ``FDHAContextMaker.is_surface_rupturing`` applied per rupture in
    ``calculate_fdha_hazard`` (curve and map paths alike)."""

    def test_tolerance_is_half_km(self):
        assert FDHAContextMaker.SURFACE_DEPTH_TOLERANCE_KM == 0.5

    def test_buried_complex_ruptures_flagged(self, complex_source):
        (src,) = complex_source.values()
        rups = list(src.iter_ruptures())
        tops = [float(np.nanmin(r.surface.mesh.depths)) for r in rups]
        buried = [r for r, t in zip(rups, tops) if t > 0.5]
        surface = [r for r, t in zip(rups, tops) if t <= 0.5]
        # the floating-rupture scenario must exercise both sides of the gate
        assert buried, "fixture has no buried ruptures - gate not exercised"
        assert surface, "fixture has no surface ruptures"
        cm = FDHAContextMaker(_sitecol(HW_SITE), {}, maximum_distance=50.0)
        assert all(not cm.is_surface_rupturing(r) for r in buried)
        assert all(cm.is_surface_rupturing(r) for r in surface)

    def test_buried_only_source_contributes_zero(self, tmp_path):
        """With iter_ruptures restricted to top edge > 0.5 km, the hazard
        must be exactly zero everywhere: every rupture is skipped before
        any probability model is evaluated."""
        calc, src = _make_calc(tmp_path, COMPLEX_SOURCE_XML, "buried")
        all_rups = list(src.iter_ruptures())
        buried = [r for r in all_rups
                  if float(np.nanmin(r.surface.mesh.depths)) > 0.5]
        assert buried
        src.iter_ruptures = lambda **kw: iter(buried)
        results = calc.run()
        assert np.asarray(results['poes']).max() == 0.0
        assert np.asarray(results['rate_principal']).max() == 0.0
        assert np.asarray(results['rate_distributed']).max() == 0.0

    def test_surface_only_source_reproduces_full_hazard(self, tmp_path):
        """Dropping the buried ruptures from the source changes nothing:
        the full-source hazard already excludes them."""
        calc_full, _ = _make_calc(tmp_path, COMPLEX_SOURCE_XML, "full")
        poes_full = np.asarray(calc_full.run()['poes'])

        calc_surf, src = _make_calc(tmp_path, COMPLEX_SOURCE_XML, "surfonly")
        surface = [r for r in src.iter_ruptures()
                   if float(np.nanmin(r.surface.mesh.depths)) <= 0.5]
        src.iter_ruptures = lambda **kw: iter(surface)
        poes_surf = np.asarray(calc_surf.run()['poes'])

        np.testing.assert_array_equal(poes_full, poes_surf)


INI_TEMPLATE = (
    "[general]\n"
    "description = complex fault end-to-end test\n"
    "calculation_mode = fdha_classical\n\n"
    "[geometry]\n"
    "sites = {sites}\n\n"
    "[site_params]\n"
    "reference_vs30_value = 760.0\n\n"
    "[erf]\n"
    "rupture_mesh_spacing = 1.0\n"
    "complex_fault_mesh_spacing = 1.0\n"
    "width_of_mfd_bin = 0.1\n\n"
    "[calculation]\n"
    "investigation_time = 1.0\n"
    'displacement_measure_levels = {{"FD": [0.001, 0.01, 0.1, 1.0]}}\n'
    "r_threshold_km = 1.0\n\n"
    "[parameters]\n"
    "target_displacement = [0.001, 0.01, 0.1, 1.0]\n\n"
    "[models.primary_surf_rup]\n"
    "type = Youngs2003PrimarySR\n\n"
    "[models.primary_surf_rup.parameters]\n"
    "style = normal\n\n"
    "[models.primary_surf_displ]\n"
    "type = Youngs2003PrimaryFD\n\n"
    "[models.primary_surf_displ.parameters]\n"
    "style = normal\n"
    "norm_disp_type = AD\n\n"
    "[models.secondary_surf_rup]\n"
    "type = Youngs2003SecondarySR\n\n"
    "[models.secondary_surf_displ]\n"
    "type = Youngs2003SecondaryFD\n"
)


def _make_calc(tmp_path, source_xml_text, name, sites="0.0 0.15",
               rupture_mesh_spacing=1.0, complex_fault_mesh_spacing=1.0,
               width_of_mfd_bin=0.1, **converterparams):
    """Build a ready-to-run calculator; returns (calculator, first source)."""
    from openquake.fdha.calc.calculators import (
        FaultRuptureProbabilityCalculator)
    src_xml = tmp_path / f"{name}.xml"
    src_xml.write_text(source_xml_text)
    branch_dir = tmp_path / "branch_configs"
    branch_dir.mkdir(parents=True, exist_ok=True)
    ini = branch_dir / f"branch_{name}.ini"
    ini.write_text(INI_TEMPLATE.format(sites=sites))
    calc = FaultRuptureProbabilityCalculator(
        str(ini), [str(src_xml)],
        rupture_mesh_spacing=rupture_mesh_spacing,
        complex_fault_mesh_spacing=complex_fault_mesh_spacing,
        width_of_mfd_bin=width_of_mfd_bin, **converterparams)
    (src,) = calc.fault_sources.values()
    return calc, src


def _run_calc(tmp_path, source_xml_text, name, sites="0.0 0.15",
              **converterparams):
    calc, _src = _make_calc(
        tmp_path, source_xml_text, name, sites=sites, **converterparams)
    return calc.run()


class TestHazardEndToEnd:
    def test_hazard_curve_standalone_complex(self, tmp_path):
        """On-trace site: principal hazard must accumulate, rates bounded by
        the GR total rate, and the low-displacement plateau must sit near
        the reference values of this fixed scenario."""
        results = _run_calc(
            tmp_path, COMPLEX_SOURCE_XML, "complex",
            rupture_mesh_spacing=1.0, complex_fault_mesh_spacing=1.0,
            width_of_mfd_bin=0.1)
        poes = np.asarray(results['poes'])
        assert poes.shape == (1, 4)
        assert np.isfinite(poes).all() and (poes >= 0.0).all()
        assert np.asarray(results['rate_principal']).max() > 0.0
        total_rate = 10 ** (4.0 - 1.0 * 6.5) - 10 ** (4.0 - 1.0 * 7.0)
        assert poes.max() <= total_rate * 1.0001
        # regression anchor: frozen from a verified run of this scenario
        # (plateau = 11% of the 2.1623e-3/yr GR total rate, and within 5%
        # of the simpleFaultSource twin — see the equivalence test below)
        expected = np.array(
            [2.38460639e-04, 2.38002779e-04, 2.22840039e-04, 8.79603903e-05])
        np.testing.assert_allclose(poes[0], expected, rtol=1e-6)

    def test_complex_matches_simple_twin(self, tmp_path):
        """The complexFaultSource written as top/bottom edges of the exact
        plane of the simpleFaultSource twin must produce the same hazard
        curve (identical MFD; only the surface construction differs)."""
        res_c = _run_calc(
            tmp_path, COMPLEX_SOURCE_XML, "complex",
            rupture_mesh_spacing=1.0, complex_fault_mesh_spacing=1.0,
            width_of_mfd_bin=0.1)
        res_s = _run_calc(
            tmp_path, SIMPLE_TWIN_XML, "simple",
            rupture_mesh_spacing=1.0, width_of_mfd_bin=0.1)
        poes_c = np.asarray(res_c['poes'])[0]
        poes_s = np.asarray(res_s['poes'])[0]
        assert (poes_s > 0.0).all()
        np.testing.assert_allclose(poes_c, poes_s, rtol=0.05)

    def test_hazard_curve_characteristic_complex(self, tmp_path):
        """Single characteristic rupture at 1e-3/yr: the on-trace plateau is
        bounded by (and close to) rate x P_sr(M7)."""
        results = _run_calc(
            tmp_path, CHAR_COMPLEX_XML, "charcomplex",
            rupture_mesh_spacing=1.0, complex_fault_mesh_spacing=1.0,
            width_of_mfd_bin=0.1)
        poes = np.asarray(results['poes'])
        assert poes.shape == (1, 4)
        assert np.asarray(results['rate_principal']).max() > 0.0
        # Youngs2003 normal-style P(SR|M7) gate
        from openquake.fdha.primary_surf_rup import Youngs2003PrimarySR
        p_sr = float(Youngs2003PrimarySR().get_prob(7.0, style="normal"))
        assert poes.max() <= 0.001 * p_sr * 1.0001
        assert poes.max() > 0.5 * 0.001 * p_sr
