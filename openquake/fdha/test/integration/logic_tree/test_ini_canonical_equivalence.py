"""Phase G tests 9–10 after legacy removal.

The pre-removal equivalence checks are now hard-stop checks for the legacy
keys plus a smoke run of the canonical replacement.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

from openquake.fdha.calc.config_loader import ConfigurationError
from openquake.fdha.logic_tree.driver import FdhaLogicTree


DATA = Path(__file__).resolve().parent / "data"


def _repo_root() -> Path:
    for p in [Path(__file__).resolve(), *Path(__file__).resolve().parents]:
        if (p / "examples").is_dir() and (p / "openquake").is_dir():
            return p
    raise RuntimeError("Could not locate repository root.")


@pytest.fixture(scope="session")
def validation_demo_root() -> Path:
    """Shared validation demo assets."""
    root = (
        _repo_root()
        / "openquake/fdha/test/fixtures/examples_archive/logic_tree_validation"
    )
    return root


def _rates_from_lt(ini_path: Path) -> np.ndarray:
    out = ini_path.parent / "out_ini_equiv_test_tmp"
    if out.exists():
        shutil.rmtree(out)
    lt = FdhaLogicTree.from_ini(str(ini_path))
    result = lt.run(outdir=out)
    rates = np.asarray(result.mean_rates, dtype=float)
    return rates


@pytest.mark.integration
def test_equivalence_source_model_file_vs_logic_tree_wrapper(
    validation_demo_root: Path, tmp_path: Path
):
    shutil.copy(
        validation_demo_root / "source_model.xml",
        tmp_path / "source_model.xml",
    )
    shutil.copy(
        validation_demo_root / "source_model_logic_tree.xml",
        tmp_path / "source_model_logic_tree.xml",
    )
    shutil.copy(
        validation_demo_root / "single_bilinear" / "fdha_logic_tree.xml",
        tmp_path / "fdha_logic_tree.xml",
    )

    lv = "[0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]"
    common_geom = """[general]
calculation_mode = fdha_classical
[geometry]
sites = 16.16878333 39.66247618
[site_params]
reference_vs30_value = 760.0
[logic_tree]
number_of_logic_tree_samples = 0
[erf]
rupture_mesh_spacing = 2.0
width_of_mfd_bin = 0.1
"""

    ini_legacy = tmp_path / "legacy_source_model_file.ini"
    ini_legacy.write_text(
        common_geom
        + f"""[calculation]
source_model_file = source_model.xml
investigation_time = 1.0
displacement_measure_levels = {{"FD": {lv}}}
fdha_logic_tree_file = fdha_logic_tree.xml
"""
    )

    ini_canon = tmp_path / "canonical_smlt.ini"
    ini_canon.write_text(
        common_geom
        + f"""[calculation]
source_model_logic_tree_file = source_model_logic_tree.xml
investigation_time = 1.0
displacement_measure_levels = {{"FD": {lv}}}
fdha_logic_tree_file = fdha_logic_tree.xml
"""
    )

    with pytest.raises(ConfigurationError, match="source_model_logic_tree_file"):
        _rates_from_lt(ini_legacy)
    r_canon = _rates_from_lt(ini_canon)

    assert np.any(r_canon >= 0)


@pytest.mark.integration
def test_equivalence_fdha_logic_tree_files_plural_vs_singular_merged(tmp_path):
    shutil.copy(DATA / "lt_equiv_frag1.xml", tmp_path / "frag1.xml")
    shutil.copy(DATA / "lt_equiv_frag2.xml", tmp_path / "frag2.xml")
    shutil.copy(DATA / "lt_equiv_full.xml", tmp_path / "full.xml")
    shutil.copy(
        _repo_root()
        / "openquake/fdha/test/fixtures/examples_archive/logic_tree_validation"
        / "source_model.xml",
        tmp_path / "source_model.xml",
    )
    smlt = """<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns:gml="http://www.opengis.net/gml"
      xmlns="http://openquake.org/xmlns/nrml/0.4">
    <logicTree logicTreeID="lt">
            <logicTreeBranchSet uncertaintyType="sourceModel" branchSetID="bs1">
                <logicTreeBranch branchID="b1">
                    <uncertaintyModel>source_model.xml</uncertaintyModel>
                    <uncertaintyWeight>1.0</uncertaintyWeight>
                </logicTreeBranch>
            </logicTreeBranchSet>
    </logicTree>
</nrml>
"""
    (tmp_path / "source_model_logic_tree.xml").write_text(smlt)

    lv = "[0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]"
    common = """[general]
calculation_mode = fdha_classical
[geometry]
sites = 16.16878333 39.66247618
[site_params]
reference_vs30_value = 760.0
[logic_tree]
number_of_logic_tree_samples = 0
[erf]
rupture_mesh_spacing = 2.0
width_of_mfd_bin = 0.1
"""

    ini_plural = tmp_path / "plural_files.ini"
    ini_plural.write_text(
        common
        + f"""[calculation]
source_model_logic_tree_file = source_model_logic_tree.xml
investigation_time = 1.0
displacement_measure_levels = {{"FD": {lv}}}
fdha_logic_tree_files = frag1.xml, frag2.xml
"""
    )

    ini_singular = tmp_path / "singular_file.ini"
    ini_singular.write_text(
        common
        + f"""[calculation]
source_model_logic_tree_file = source_model_logic_tree.xml
investigation_time = 1.0
displacement_measure_levels = {{"FD": {lv}}}
fdha_logic_tree_file = full.xml
"""
    )

    with pytest.raises(ConfigurationError, match="fdha_logic_tree_file"):
        _rates_from_lt(ini_plural)
    r_sg = _rates_from_lt(ini_singular)

    assert np.any(r_sg >= 0)
