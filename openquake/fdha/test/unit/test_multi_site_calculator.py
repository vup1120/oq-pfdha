# -*- coding: utf-8 -*-
"""
Tests for multi-site hazard-curve support.

Two groups, marked separately:
- TestInitializeSite (unit): FaultRuptureProbabilityCalculator._initialize_site()
  building a multi-point SiteCollection from site_location.sites_list, plus the
  single-site backward-compatible path.
- TestMultiSiteHazardCurve (integration / regression): end-to-end runs verifying
  a multi-site job reproduces, per site, the standalone single-site result.
"""

import os
import numpy as np
import pytest
from numpy.testing import assert_allclose

from openquake.fdha.calc.calculators import FaultRuptureProbabilityCalculator


# Norcia benchmark fault + source model (real geometry, frozen science).
# Resolved relative to this test file so the test is cwd-independent.
NORCIA_SOURCE_MODEL = os.path.normpath(os.path.join(
    os.path.dirname(__file__),
    "..", "benchmark", "norcia_sensitivity_youngs2003", "source_model_norcia.xml",
))


def _make_uninitialised_calculator(config):
    """Return a FaultRuptureProbabilityCalculator with only .config set.

    Bypasses __init__ so _initialize_site() can be exercised in isolation,
    without parsing source models or instantiating displacement models.
    """
    calc = FaultRuptureProbabilityCalculator.__new__(FaultRuptureProbabilityCalculator)
    calc.config = config
    calc.sitecol = None
    return calc


class TestInitializeSite:
    """FaultRuptureProbabilityCalculator._initialize_site() - site construction."""

    pytestmark = pytest.mark.unit

    def test_initialize_site_multi_site(self):
        calc = _make_uninitialised_calculator({
            'site_location': {
                'sites_list': [
                    {'longitude': 13.0, 'latitude': 42.0},
                    {'longitude': 13.5, 'latitude': 42.5},
                    {'longitude': 14.0, 'latitude': 43.0},
                ]
            }
        })
        calc._initialize_site()
        assert len(calc.sitecol) == 3
        assert_allclose(calc.sitecol.lons, [13.0, 13.5, 14.0])
        assert_allclose(calc.sitecol.lats, [42.0, 42.5, 43.0])
        # Site ids are 0..N-1 - required by the np.add.at accumulation path
        assert_allclose(calc.sitecol.sids, [0, 1, 2])

    def test_initialize_site_single_site_compat(self):
        """A config with scalar lon/lat (no sites_list) still yields 1 site."""
        calc = _make_uninitialised_calculator({
            'site_location': {'longitude': 13.278, 'latitude': 42.767}
        })
        calc._initialize_site()
        assert len(calc.sitecol) == 1
        assert_allclose(calc.sitecol.lons, [13.278])
        assert_allclose(calc.sitecol.lats, [42.767])

    def test_initialize_site_multi_site_global_vs30(self):
        """Global site_location.vs30 is applied to every site."""
        calc = _make_uninitialised_calculator({
            'site_location': {
                'vs30': 800.0,
                'sites_list': [
                    {'longitude': 13.0, 'latitude': 42.0},
                    {'longitude': 13.5, 'latitude': 42.5},
                ],
            }
        })
        calc._initialize_site()
        assert len(calc.sitecol) == 2
        assert_allclose(calc.sitecol.vs30, [800.0, 800.0])

    def test_initialize_site_missing_coords_raises(self):
        calc = _make_uninitialised_calculator({'site_location': {}})
        with pytest.raises(ValueError, match="latitude and longitude"):
            calc._initialize_site()


def _write_multi_site_branch_ini(base_dir, sites_line):
    """Write a materialised branch INI (under branch_configs/) with a custom
    [geometry].sites line, based on the Norcia Youngs2003 AD-85 branch.

    The branch_configs/ path component lets the loader accept the legacy
    [models.*] sections that FaultRuptureProbabilityCalculator consumes.
    ``base_dir`` may be a not-yet-existing subdirectory of tmp_path, so that
    several independent jobs can be written under one test.
    """
    branch_dir = base_dir / "branch_configs"
    branch_dir.mkdir(parents=True, exist_ok=True)
    ini = branch_dir / "branch_multisite.ini"
    ini.write_text(
        "[general]\n"
        "description = multi-site test (Norcia Youngs2003 AD 85)\n"
        "calculation_mode = fdha_classical\n\n"
        "[geometry]\n"
        f"{sites_line}\n\n"
        "[site_params]\n"
        "reference_vs30_value = 760.0\n\n"
        "[erf]\n"
        "rupture_mesh_spacing = 0.1\n"
        "width_of_mfd_bin = 0.1\n\n"
        "[calculation]\n"
        "investigation_time = 1.0\n"
        'displacement_measure_levels = {"FD": [0.001, 0.01, 0.05, 0.1, 0.5, 1.0]}\n'
        "r_threshold_km = 0.1\n\n"
        "[parameters]\n"
        "target_displacement = [0.001, 0.01, 0.05, 0.1, 0.5, 1.0]\n\n"
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


# Two Norcia test sites:
#   site 0 - on a fault-trace vertex (13.2264, 42.8846), principal zone
#   site 1 - same latitude, ~0.27 deg (~22 km) east, distributed only
_SITE_0 = "13.2264 42.8846"
_SITE_1 = "13.5000 42.8846"


def _run_curve(base_dir, sites_line):
    """Run an end-to-end hazard curve for the given [geometry].sites line."""
    assert os.path.exists(NORCIA_SOURCE_MODEL), (
        f"Norcia source model missing: {NORCIA_SOURCE_MODEL}"
    )
    ini = _write_multi_site_branch_ini(base_dir, f"sites = {sites_line}")
    calc = FaultRuptureProbabilityCalculator(str(ini), [NORCIA_SOURCE_MODEL])
    return calc, calc.run()


class TestMultiSiteHazardCurve:
    """End-to-end multi-site hazard curve via FaultRuptureProbabilityCalculator."""

    pytestmark = [
        pytest.mark.integration,
        pytest.mark.regression,
        pytest.mark.youngs2003,
    ]

    def test_multi_site_matches_separate_single_site_runs(self, tmp_path):
        """A multi-site run must reproduce, per site, the standalone run of
        each site.

        Site 0 and site 1 are each calculated alone end-to-end, then together
        in one multi-site job. Row ``i`` of the multi-site ``rates`` must match
        the standalone run of site ``i`` to floating-point tolerance - adding
        more sites to a job must not perturb any site's numerics.
        """
        # Standalone single-site runs
        _, result_0 = _run_curve(tmp_path / "s0", _SITE_0)
        _, result_1 = _run_curve(tmp_path / "s1", _SITE_1)
        rates_0 = np.asarray(result_0["rates"])
        rates_1 = np.asarray(result_1["rates"])

        # Combined multi-site run
        calc_multi, result_multi = _run_curve(
            tmp_path / "multi", f"{_SITE_0}, {_SITE_1}"
        )
        rates_multi = np.asarray(result_multi["rates"])

        assert rates_0.shape == (1, 6)
        assert rates_1.shape == (1, 6)
        assert rates_multi.shape == (2, 6)
        assert len(calc_multi.sitecol) == 2

        # Per-site numerics must be identical whether run alone or together
        assert_allclose(rates_multi[0], rates_0[0], rtol=1e-10, atol=1e-14)
        assert_allclose(rates_multi[1], rates_1[0], rtol=1e-10, atol=1e-14)

        # Sanity: the near-fault site produced a non-trivial curve
        assert np.any(rates_multi[0] > 0.0), "near-fault site curve is all zero"

    def test_multi_site_result_metadata(self, tmp_path):
        """Multi-site result carries per-site coordinates aligned with rates rows."""
        _, result = _run_curve(tmp_path, f"{_SITE_0}, {_SITE_1}")

        assert result["n_sites"] == 2
        assert_allclose(result["site_lons"], [13.2264, 13.5000])
        assert_allclose(result["site_lats"], [42.8846, 42.8846])


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
