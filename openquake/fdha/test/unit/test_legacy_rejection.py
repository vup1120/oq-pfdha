from __future__ import annotations

import pytest

from openquake.fdha.calc.config_loader import ConfigurationError, load_config


def _write_ini(tmp_path, body: str):
    path = tmp_path / "job.ini"
    path.write_text(body)
    return path


def test_source_model_file_is_rejected(tmp_path):
    ini = _write_ini(
        tmp_path,
        """[calculation]
source_model_file = source_model.xml
fdha_logic_tree_file = fdha.xml
""",
    )

    with pytest.raises(ConfigurationError) as exc:
        load_config(ini)

    msg = str(exc.value)
    assert "Legacy key `source_model_file` is no longer supported" in msg
    assert "source_model_logic_tree_file" in msg


def test_fdha_logic_tree_files_plural_is_rejected(tmp_path):
    ini = _write_ini(
        tmp_path,
        """[calculation]
source_model_logic_tree_file = source_model_logic_tree.xml
fdha_logic_tree_files = a.xml, b.xml
""",
    )

    with pytest.raises(ConfigurationError) as exc:
        load_config(ini)

    msg = str(exc.value)
    assert "Legacy key `fdha_logic_tree_files` (plural) is no longer supported" in msg
    assert "fdha_logic_tree_file" in msg


def test_models_section_is_rejected(tmp_path):
    ini = _write_ini(
        tmp_path,
        """[calculation]
source_model_logic_tree_file = source_model_logic_tree.xml
fdha_logic_tree_file = fdha.xml

[models.primary_surf_rup]
type = Youngs2003PrimarySR
""",
    )

    with pytest.raises(ConfigurationError) as exc:
        load_config(ini)

    msg = str(exc.value)
    assert "Legacy `[models]` section is no longer supported" in msg
    assert "uncertaintyModel" in msg


def test_target_displacement_is_rejected(tmp_path):
    ini = _write_ini(
        tmp_path,
        """[calculation]
source_model_logic_tree_file = source_model_logic_tree.xml
fdha_logic_tree_file = fdha.xml

[parameters]
target_displacement = [0.01]
""",
    )

    with pytest.raises(ConfigurationError) as exc:
        load_config(ini)

    msg = str(exc.value)
    assert "Legacy key `target_displacement` is no longer supported" in msg
    assert "displacement_measure_levels" in msg
