"""
Unit tests: Configuration loader and validation

Tests the config_loader module functionality, including:
- Configuration loading (INI parsing, multi-site geometry)
- Legacy entry-point rejection
- Coordinate validation
"""

import pytest

pytestmark = pytest.mark.unit
from openquake.fdha.calc.config_loader import (
    load_config,
    ConfigurationError,
    ConfigValidationError
)


class TestLegacyTomlConfiguration:
    """TOML is an explicit hard stop."""

    def test_load_config_rejects_toml(self, tmp_path):
        toml_file = tmp_path / "config.toml"
        toml_file.write_text("[calculation]\n")

        with pytest.raises(ConfigurationError, match="TOML configuration files are no longer supported"):
            load_config(toml_file)


def _write_ini(tmp_path, geometry_body):
    """Write a minimal valid INI with a custom [geometry] body and return its path."""
    ini = tmp_path / "job.ini"
    ini.write_text(
        "[general]\n"
        "description = multi_site_test\n"
        "calculation_mode = fdha_classical\n\n"
        "[geometry]\n"
        f"{geometry_body}\n\n"
        "[erf]\n"
        "rupture_mesh_spacing = 2.0\n"
        "width_of_mfd_bin = 0.1\n\n"
        "[calculation]\n"
        "investigation_time = 1.0\n"
        'displacement_measure_levels = {"FD": [0.001, 0.01, 0.1]}\n'
        "r_threshold_km = 0.1\n"
        "fdha_logic_tree_file = lt.xml\n"
        "source_model_logic_tree_file = smlt.xml\n"
    )
    return ini


class TestMultiSiteParsing:
    """Multi-site parsing of [geometry].sites and [geometry].sites_csv."""

    def test_single_site_inline_unchanged(self, tmp_path):
        """Single-site 'sites = lon lat' still sets scalar lon/lat (regression)."""
        ini = _write_ini(tmp_path, "sites = 16.16573727 39.64704451")
        config = load_config(ini)
        site_loc = config['site_location']
        assert site_loc['longitude'] == pytest.approx(16.16573727)
        assert site_loc['latitude'] == pytest.approx(39.64704451)

    def test_multi_site_inline_two_sites(self, tmp_path):
        ini = _write_ini(tmp_path, "sites = 10.0 43.0, 10.1 43.1")
        config = load_config(ini)
        sites_list = config['site_location']['sites_list']
        assert len(sites_list) == 2
        assert sites_list[0] == {'longitude': 10.0, 'latitude': 43.0}
        assert sites_list[1] == {'longitude': 10.1, 'latitude': 43.1}

    def test_multi_site_inline_depth(self, tmp_path):
        ini = _write_ini(tmp_path, "sites = 10.0 43.0 0.5, 10.1 43.1 1.0")
        config = load_config(ini)
        sites_list = config['site_location']['sites_list']
        assert len(sites_list) == 2
        assert sites_list[0]['depth'] == pytest.approx(0.5)
        assert sites_list[1]['depth'] == pytest.approx(1.0)

    def test_single_site_still_sets_scalars(self, tmp_path):
        """sites_list of length 1 also exposes scalar lon/lat for backward compat."""
        ini = _write_ini(tmp_path, "sites = 10.0 43.0")
        config = load_config(ini)
        site_loc = config['site_location']
        assert len(site_loc['sites_list']) == 1
        assert site_loc['longitude'] == pytest.approx(10.0)
        assert site_loc['latitude'] == pytest.approx(43.0)

    def test_multi_site_sets_first_as_scalar(self, tmp_path):
        """sites_list of length > 1 sets scalar lon/lat to the first site."""
        ini = _write_ini(tmp_path, "sites = 10.0 43.0, 11.0 44.0")
        config = load_config(ini)
        site_loc = config['site_location']
        assert site_loc['longitude'] == pytest.approx(10.0)
        assert site_loc['latitude'] == pytest.approx(43.0)

    def test_sites_csv_two_sites(self, tmp_path):
        csv = tmp_path / "sites.csv"
        csv.write_text("site_id,lon,lat\n0,10.0,43.0\n1,10.1,43.1\n")
        ini = _write_ini(tmp_path, "sites_csv = sites.csv")
        config = load_config(ini)
        sites_list = config['site_location']['sites_list']
        assert len(sites_list) == 2
        assert sites_list[0] == {'longitude': 10.0, 'latitude': 43.0}
        assert sites_list[1] == {'longitude': 10.1, 'latitude': 43.1}
        assert config['site_location']['longitude'] == pytest.approx(10.0)
        assert config['site_location']['latitude'] == pytest.approx(43.0)

    def test_sites_csv_with_depth(self, tmp_path):
        csv = tmp_path / "sites.csv"
        csv.write_text("site_id,lon,lat,depth\n0,10.0,43.0,0.5\n1,10.1,43.1,1.0\n")
        ini = _write_ini(tmp_path, "sites_csv = sites.csv")
        config = load_config(ini)
        sites_list = config['site_location']['sites_list']
        assert sites_list[0]['depth'] == pytest.approx(0.5)
        assert sites_list[1]['depth'] == pytest.approx(1.0)

    def test_sites_csv_missing_column_raises(self, tmp_path):
        csv = tmp_path / "sites.csv"
        csv.write_text("site_id,longitude,lat\n0,10.0,43.0\n")
        ini = _write_ini(tmp_path, "sites_csv = sites.csv")
        with pytest.raises(ConfigValidationError, match="lon.*lat"):
            load_config(ini)

    def test_multi_site_lon_out_of_range(self, tmp_path):
        ini = _write_ini(tmp_path, "sites = 400.0 43.0, 10.1 43.1")
        with pytest.raises(ConfigValidationError, match="Longitude"):
            load_config(ini)

    def test_multi_site_lat_out_of_range(self, tmp_path):
        ini = _write_ini(tmp_path, "sites = 10.0 95.0, 10.1 43.1")
        with pytest.raises(ConfigValidationError, match="Latitude"):
            load_config(ini)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
