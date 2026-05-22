"""Legacy public INI forms are hard errors after v5 canonicalisation."""
from __future__ import annotations

from pathlib import Path

import pytest


def _repo_root() -> Path:
    for p in [Path(__file__).resolve(), *Path(__file__).resolve().parents]:
        if (p / "examples").is_dir() and (p / "openquake").is_dir():
            return p
    raise RuntimeError("Could not locate repository root.")


def test_legacy_single_model_job_raises_configuration_error(tmp_path, monkeypatch):
    legacy_ini = tmp_path / "legacy.ini"
    legacy_ini.write_text(
        """[general]
description = legacy
[geometry]
sites = 0 0
[erf]
rupture_mesh_spacing = 1.0
width_of_mfd_bin = 0.1
[calculation]
displacement_measure_levels = {"FD": [0.01]}
[models.primary_surf_rup]
type = Youngs2003PrimarySR
"""
    )

    from openquake.fdha.calc.config_loader import ConfigurationError
    from openquake.fdha.main import run_calculation

    with pytest.raises(ConfigurationError, match="models"):
        run_calculation(str(legacy_ini))


def test_missing_both_logic_tree_and_models_raises_configuration_error(tmp_path):
    ini = tmp_path / "empty.ini"
    ini.write_text(
        "[general]\n"
        "description = empty\n"
        "[geometry]\n"
        "sites = 0 0\n"
        "[erf]\n"
        "rupture_mesh_spacing = 1.0\n"
        "width_of_mfd_bin = 0.1\n"
        "[calculation]\n"
        "fdha_logic_tree_file = lt.xml\n"
    )

    from openquake.fdha.calc.config_loader import ConfigurationError
    from openquake.fdha.main import run_calculation

    with pytest.raises(ConfigurationError):
        run_calculation(str(ini))
