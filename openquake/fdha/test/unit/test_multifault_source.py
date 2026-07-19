# -*- coding: utf-8 -*-
"""
End-to-end tests for multiFaultSource support.

A synthetic 2-section (kite) multiFaultSource is written as a NRML source +
geometryModel pair and driven through the whole chain:

- ``parse_source_model_faults``: investigation_time sniffing from the XML
  header and automatic auxiliary hdf5 provisioning;
- ``MultiFaultSource.iter_ruptures``: non-parametric ruptures with PMF
  occurrence probabilities on single- and multi-section MultiSurfaces;
- ``FDHAContextMaker.get_ctx``: PMF -> Poisson-equivalent annual rate
  conversion, coarse far-rupture pre-filter, section-id surface hashing,
  reference-line (ECS/LCP) x/L for multi-section ruptures and the plain
  top-trace path for single-section ruptures;
- ``calculate_fdha_hazard`` via ``FaultRuptureProbabilityCalculator``.
"""

import numpy as np
import pytest

from openquake.fdha.calc.utils.parsing import parse_source_model_faults
from openquake.fdha.calc.contexts import FDHAContextMaker
from openquake.hazardlib.site import Site, SiteCollection
from openquake.hazardlib.geo import Point

pytestmark = pytest.mark.unit


# Two near-vertical kite sections striking ~E-W at lat 42, top at 0 km depth,
# arranged en echelon with a ~2.5 km gap / ~1.1 km stepover between them.
SECTIONS_XML = """<?xml version='1.0' encoding='utf-8'?>
<nrml xmlns:gml="http://www.opengis.net/gml"
      xmlns="http://openquake.org/xmlns/nrml/0.5">
  <geometryModel name="test sections" investigation_time="1.0">
    <section name="S1" id="s1">
      <kiteSurface>
        <profile>
          <gml:LineString><gml:posList>
            13.00 42.00 0.0 13.00 42.02 10.0
          </gml:posList></gml:LineString>
        </profile>
        <profile>
          <gml:LineString><gml:posList>
            13.12 42.00 0.0 13.12 42.02 10.0
          </gml:posList></gml:LineString>
        </profile>
      </kiteSurface>
    </section>
    <section name="S2" id="s2">
      <kiteSurface>
        <profile>
          <gml:LineString><gml:posList>
            13.15 42.01 0.0 13.15 42.03 10.0
          </gml:posList></gml:LineString>
        </profile>
        <profile>
          <gml:LineString><gml:posList>
            13.27 42.01 0.0 13.27 42.03 10.0
          </gml:posList></gml:LineString>
        </profile>
      </kiteSurface>
    </section>
  </geometryModel>
</nrml>
"""

# Three non-parametric ruptures: each single section, and both together.
SOURCES_XML = """<?xml version='1.0' encoding='utf-8'?>
<nrml xmlns:gml="http://www.opengis.net/gml"
      xmlns="http://openquake.org/xmlns/nrml/0.5">
  <sourceModel name="test model" investigation_time="1.0">
    <sourceGroup name="g1" rup_interdep="indep" src_interdep="indep"
                 tectonicRegion="Active Shallow Crust">
      <multiFaultSource id="mf1" name="Test multifault">
        <multiPlanesRupture probs_occur="0.99 0.01">
          <magnitude>6.0</magnitude>
          <sectionIndexes indexes="s1"/>
          <rake>-90</rake>
        </multiPlanesRupture>
        <multiPlanesRupture probs_occur="0.98 0.02">
          <magnitude>6.1</magnitude>
          <sectionIndexes indexes="s2"/>
          <rake>-90</rake>
        </multiPlanesRupture>
        <multiPlanesRupture probs_occur="0.995 0.005">
          <magnitude>6.5</magnitude>
          <sectionIndexes indexes="s1,s2"/>
          <rake>-90</rake>
        </multiPlanesRupture>
      </multiFaultSource>
    </sourceGroup>
  </sourceModel>
</nrml>
"""


@pytest.fixture(scope="module")
def multifault_source(tmp_path_factory):
    base = tmp_path_factory.mktemp("mfsrc")
    src_xml = base / "sources.xml"
    sec_xml = base / "sections.xml"
    src_xml.write_text(SOURCES_XML)
    sec_xml.write_text(SECTIONS_XML)
    # No hdf5path and no investigation_time given on purpose: both must be
    # provided automatically (temp file / sniffed from the XML header).
    srcs = parse_source_model_faults(
        [str(src_xml), str(sec_xml)], rupture_mesh_spacing=1.0)
    return srcs


@pytest.fixture(scope="module")
def ruptures(multifault_source):
    (src,) = multifault_source.values()
    return src, list(src.iter_ruptures())


def _sitecol(*lonlats):
    return SiteCollection([
        Site(location=Point(lon, lat), vs30=760.0) for lon, lat in lonlats])


class TestParsing:
    def test_parse_and_iter(self, ruptures):
        src, rups = ruptures
        assert src.investigation_time == 1.0  # sniffed from the XML header
        assert len(rups) == 3
        # engine yields non-parametric ruptures carrying a PMF, no rate
        for rup in rups:
            assert not hasattr(rup, 'occurrence_rate')
            assert rup.probs_occur.sum() == pytest.approx(1.0)
        assert [len(r.surface.surfaces) for r in rups] == [1, 1, 2]


class TestGetCtx:
    def test_pmf_to_annual_rate(self, ruptures):
        src, rups = ruptures
        sitecol = _sitecol((13.06, 42.005))
        cm = FDHAContextMaker(sitecol, {}, maximum_distance=50.0)
        for rup in rups:
            ctx = cm.get_ctx(rup, investigation_time=src.investigation_time)
            assert ctx is not None
            expected = -np.log(rup.probs_occur[0]) / src.investigation_time
            assert ctx.occurrence_rate[0] == pytest.approx(expected)

    def test_investigation_time_scales_rate(self, ruptures):
        _src, rups = ruptures
        sitecol = _sitecol((13.06, 42.005))
        cm = FDHAContextMaker(sitecol, {}, maximum_distance=50.0)
        r1 = cm.get_ctx(rups[0], investigation_time=1.0).occurrence_rate[0]
        r50 = cm.get_ctx(rups[0], investigation_time=50.0).occurrence_rate[0]
        assert r50 == pytest.approx(r1 / 50.0)

    def test_single_section_rupture_context(self, ruptures):
        """Single-section MultiSurface: trace = section top edge, no ECS."""
        src, rups = ruptures
        # ~0.5 km north of the middle of s1
        sitecol = _sitecol((13.06, 42.005))
        cm = FDHAContextMaker(sitecol, {}, maximum_distance=50.0)
        ctx = cm.get_ctx(rups[0], investigation_time=src.investigation_time)
        assert ctx is not None
        assert 0.3 < ctx.r[0] < 1.0
        # site is mid-section: x/L well inside the trace
        assert 0.3 < ctx.x_L[0] < 0.7
        assert 8.0 < ctx.L[0] < 11.0
        assert ctx.ztor[0] == pytest.approx(0.0, abs=0.01)

    def test_multi_section_rupture_context(self, ruptures):
        src, rups = ruptures
        sitecol = _sitecol((13.06, 42.005),   # near s1
                           (13.135, 42.006))  # inside the stepover gap
        cm = FDHAContextMaker(sitecol, {}, maximum_distance=50.0)
        ctx = cm.get_ctx(rups[2], investigation_time=src.investigation_time)
        assert ctx is not None and len(ctx) == 2
        # reference line spans both sections plus the gap
        assert 18.0 < ctx.L[0] < 26.0
        assert np.all((ctx.x_L >= 0.0) & (ctx.x_L <= 1.0))
        # the gap site is close to the reference line bridge
        assert ctx.r[1] < 2.0
        assert np.isfinite(ctx.rx).all()

    def test_per_model_reference_line_metrics(self, ruptures):
        """The context carries one metric set per reference-line method in
        the union declared by the configured models (plus 'segments' for the
        canonical mask r), and metrics_for() routes models to their own set:
        e.g. Chiou2025 ('ecs') and Visini2025 ('segments') in the same job."""
        src, rups = ruptures
        sitecol = _sitecol((13.06, 42.005),   # near s1
                           (13.135, 42.006))  # inside the stepover gap
        cm = FDHAContextMaker(
            sitecol, {'multifault_reference_lines': ('lcp', 'ecs')},
            maximum_distance=50.0)
        ctx = cm.get_ctx(rups[2], investigation_time=src.investigation_time)
        assert ctx is not None
        assert set(ctx.ref_metrics) == {'lcp', 'ecs', 'segments'}
        # canonical arrays are the smoothed default ('lcp') for x/L ...
        r_lcp, xl_lcp, L_lcp = ctx.metrics_for('lcp')
        np.testing.assert_array_equal(ctx.x_L, xl_lcp)
        np.testing.assert_array_equal(ctx.L, L_lcp)
        # ... but ctx.r (mask / near-far) is the segments distance: the gap
        # site keeps its true ~1.2 km to the nearest section, while the
        # smoothed lcp line bridges the gap (r ~ 0)
        r_seg, _, _ = ctx.metrics_for('segments')
        np.testing.assert_array_equal(ctx.r, r_seg)
        assert r_seg[1] > 3 * r_lcp[1]
        assert 0.8 < r_seg[1] < 2.0
        # every method spans the two-section system
        for method in ('lcp', 'ecs'):
            _, xl_m, L_m = ctx.metrics_for(method)
            assert 18.0 < L_m[0] < 26.0
            assert np.all((xl_m >= 0.0) & (xl_m <= 1.0))
        # an undeclared method falls back to the canonical arrays
        r_fb, xl_fb, _ = ctx.metrics_for('unknown-method')
        np.testing.assert_array_equal(r_fb, ctx.r)
        np.testing.assert_array_equal(xl_fb, ctx.x_L)

    def test_far_rupture_prefilter(self, ruptures):
        """A far site must be screened out before any reference-line build."""
        src, rups = ruptures
        cm = FDHAContextMaker(_sitecol((20.0, 50.0)), {}, maximum_distance=50.0)
        assert cm._all_sites_far(rups[2].surface)
        assert cm.get_ctx(rups[2], investigation_time=1.0) is None
        # screened before _get_distance_calculator: no cache activity at all
        assert cm.get_cache_stats()['misses'] == 0

    def test_surface_hash_by_section_ids(self, ruptures):
        _src, rups = ruptures
        cm = FDHAContextMaker(_sitecol((13.06, 42.005)), {})
        h0 = cm._get_surface_hash(rups[0].surface)
        h1 = cm._get_surface_hash(rups[1].surface)
        h2 = cm._get_surface_hash(rups[2].surface)
        assert len({h0, h1, h2}) == 3

    def test_is_surface_rupturing(self, ruptures):
        _src, rups = ruptures
        cm = FDHAContextMaker(_sitecol((13.06, 42.005)), {})
        for rup in rups:
            assert cm.is_surface_rupturing(rup)  # sections reach 0 km


class TestHazardEndToEnd:
    def _write_ini(self, base_dir, sites_line):
        branch_dir = base_dir / "branch_configs"
        branch_dir.mkdir(parents=True, exist_ok=True)
        ini = branch_dir / "branch_multifault.ini"
        ini.write_text(
            "[general]\n"
            "description = multifault end-to-end test\n"
            "calculation_mode = fdha_classical\n\n"
            "[geometry]\n"
            f"sites = {sites_line}\n\n"
            "[site_params]\n"
            "reference_vs30_value = 760.0\n\n"
            "[erf]\n"
            "rupture_mesh_spacing = 1.0\n"
            "width_of_mfd_bin = 0.1\n\n"
            "[calculation]\n"
            "investigation_time = 1.0\n"
            'displacement_measure_levels = {"FD": [0.001, 0.01, 0.1, 1.0]}\n'
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
            "[models.secondary_surf_rup.parameters]\n"
            "version = 3\n\n"
            "[models.secondary_surf_displ]\n"
            "type = Youngs2003SecondaryFD\n\n"
            "[models.secondary_surf_displ.parameters]\n"
            "percentile = 85\n"
        )
        return ini

    def test_hazard_curve_multifault(self, tmp_path):
        from openquake.fdha.calc.calculators import (
            FaultRuptureProbabilityCalculator)

        src_xml = tmp_path / "sources.xml"
        sec_xml = tmp_path / "sections.xml"
        src_xml.write_text(SOURCES_XML)
        sec_xml.write_text(SECTIONS_XML)
        # site sits on the s1 top trace -> principal zone for s1 ruptures
        ini = self._write_ini(tmp_path, "13.06 42.0")

        calc = FaultRuptureProbabilityCalculator(
            str(ini), [str(src_xml), str(sec_xml)])
        # [erf] rupture_mesh_spacing from the ini reached the converter:
        # section tor length ~9-10 km (5 km sampling would truncate it)
        (src,) = calc.fault_sources.values()
        results = calc.run()

        rates = np.asarray(results['rates'])
        assert rates.shape == (1, 4)
        assert np.isfinite(rates).all()
        assert (rates >= 0.0).all()
        # the on-trace site must accumulate a nonzero principal hazard rate
        assert np.asarray(results['rate_principal']).max() > 0.0
        # rates are bounded by the total annual occurrence rate of the model
        total_rate = sum(-np.log(p[0]) for p in src.probs_occur)
        assert rates.max() <= total_rate * 1.0001
