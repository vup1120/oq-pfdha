# -*- coding: utf-8 -*-
"""
Regression tests for the PMF investigation_time resolution.

The PMF (``probs_occur``) of a non-parametric source is defined over the
time span declared in its NRML header, which the engine converter stores on
the source object (``src.investigation_time``). The hazard loop must use
that per-source value - reading only the INI (default 1.0) silently rescales
every non-parametric rate (an XML span of 50 yr read as 1 yr inflates all
rates 50x). A conflicting INI value must be rejected loudly.
"""
import numpy as np
import pytest
from types import SimpleNamespace

from openquake.fdha.calc.hazard import _resolve_investigation_time
from openquake.fdha.test.unit.test_multifault_source import (
    SECTIONS_XML, SOURCES_XML)

pytestmark = pytest.mark.unit


# --------------------------------------------------------------- unit level
class TestResolveInvestigationTime:
    def test_source_value_is_authoritative(self):
        src = SimpleNamespace(investigation_time=2.0, source_id="mf1")
        assert _resolve_investigation_time(src, None) == 2.0

    def test_matching_ini_value_is_accepted(self):
        src = SimpleNamespace(investigation_time=2.0, source_id="mf1")
        assert _resolve_investigation_time(src, 2.0) == 2.0

    def test_conflicting_ini_value_raises(self):
        src = SimpleNamespace(investigation_time=50.0, source_id="mf1")
        with pytest.raises(ValueError, match="investigation_time"):
            _resolve_investigation_time(src, 1.0)

    def test_parametric_source_falls_back_to_ini(self):
        src = SimpleNamespace(source_id="cf1")  # no attribute
        assert _resolve_investigation_time(src, 5.0) == 5.0

    def test_parametric_source_defaults_to_one_year(self):
        src = SimpleNamespace(source_id="cf1")
        assert _resolve_investigation_time(src, None) == 1.0


# --------------------------------------------------------------- end-to-end
def _write_ini(base_dir, investigation_time_line=""):
    # [models.*] sections are only legal in materialised branch INIs, which
    # the loader recognises by the branch_configs/ parent directory.
    branch_dir = base_dir / "branch_configs"
    branch_dir.mkdir(parents=True, exist_ok=True)
    ini = branch_dir / "job.ini"
    ini.write_text(
        "[general]\n"
        "description = investigation_time resolution test\n"
        "calculation_mode = fdha_classical\n\n"
        "[geometry]\n"
        "sites = 13.06 42.0\n\n"
        "[site_params]\n"
        "reference_vs30_value = 760.0\n\n"
        "[erf]\n"
        "rupture_mesh_spacing = 1.0\n"
        "width_of_mfd_bin = 0.1\n\n"
        "[calculation]\n"
        f"{investigation_time_line}"
        'displacement_measure_levels = {"FD": [0.001, 0.01, 0.1, 1.0]}\n'
        "r_threshold_km = 1.0\n\n"
        "[models.primary_surf_rup]\n"
        "type = Youngs2003PrimarySR\n\n"
        "[models.primary_surf_rup.parameters]\n"
        "style = normal\n\n"
        "[models.primary_surf_displ]\n"
        "type = Youngs2003PrimaryFD\n\n"
        "[models.primary_surf_displ.parameters]\n"
        "style = normal\n"
        "norm_disp_type = AD\n"
    )
    return ini


def _run_multifault(tmp_path, xml_time: str, ini_time_line: str):
    from openquake.fdha.calc.calculators import (
        FaultRuptureProbabilityCalculator)

    workdir = tmp_path / f"t{xml_time}_{abs(hash(ini_time_line))}"
    workdir.mkdir()
    src_xml = workdir / "sources.xml"
    sec_xml = workdir / "sections.xml"
    src_xml.write_text(
        SOURCES_XML.replace('investigation_time="1.0"',
                            f'investigation_time="{xml_time}"'))
    sec_xml.write_text(
        SECTIONS_XML.replace('investigation_time="1.0"',
                             f'investigation_time="{xml_time}"'))
    ini = _write_ini(workdir, ini_time_line)
    calc = FaultRuptureProbabilityCalculator(
        str(ini), [str(src_xml), str(sec_xml)])
    return calc


def test_xml_time_span_scales_rates_without_ini_key(tmp_path):
    """XML investigation_time=2.0 with no INI key must halve the rates
    relative to the 1.0-yr model (rate = -ln(P0) / t)."""
    calc1 = _run_multifault(tmp_path, "1.0", "")
    calc2 = _run_multifault(tmp_path, "2.0", "")
    rates1 = np.asarray(calc1.run()["poes"], dtype=float)
    rates2 = np.asarray(calc2.run()["poes"], dtype=float)
    assert rates1.max() > 0.0, "baseline rates are zero - vacuous test"
    np.testing.assert_allclose(rates2, rates1 * 0.5, rtol=1e-10)


def test_conflicting_ini_time_raises(tmp_path):
    """XML declares 2.0 yr; an INI stating 1.0 yr must be rejected, not
    silently double every rate."""
    calc = _run_multifault(tmp_path, "2.0", "investigation_time = 1.0\n")
    with pytest.raises(ValueError, match="investigation_time"):
        calc.run()


def test_matching_ini_time_is_accepted(tmp_path):
    calc_match = _run_multifault(tmp_path, "2.0", "investigation_time = 2.0\n")
    calc_bare = _run_multifault(tmp_path, "2.0", "")
    r_match = np.asarray(calc_match.run()["poes"], dtype=float)
    r_bare = np.asarray(calc_bare.run()["poes"], dtype=float)
    np.testing.assert_allclose(r_match, r_bare, rtol=1e-12)
