from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np

from openquake.fdha.logic_tree.driver import FdhaLogicTree


def _write_lt_xml(path: Path, primary_sr: str, primary_fd: str, secondary_sr: str, secondary_fd: str) -> None:
    path.write_text(
        f"""<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="pfdha">
    <logicTreeBranchingLevel branchingLevelID="bl1">
      <logicTreeBranchSet branchSetID="bs1" uncertaintyType="fdhaPrimarySRModel">
        <logicTreeBranch branchID="B1">
          <uncertaintyModel>{primary_sr}</uncertaintyModel>
          <uncertaintyWeight>1.0</uncertaintyWeight>
        </logicTreeBranch>
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
    <logicTreeBranchingLevel branchingLevelID="bl2">
      <logicTreeBranchSet branchSetID="bs2" uncertaintyType="fdhaPrimaryFDModel" applyToBranches="B1">
        <logicTreeBranch branchID="B2">
          <uncertaintyModel>{primary_fd}</uncertaintyModel>
          <uncertaintyWeight>1.0</uncertaintyWeight>
        </logicTreeBranch>
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
    <logicTreeBranchingLevel branchingLevelID="bl3">
      <logicTreeBranchSet branchSetID="bs3" uncertaintyType="fdhaSecondarySRModel" applyToBranches="B1">
        <logicTreeBranch branchID="B3">
          <uncertaintyModel>{secondary_sr}</uncertaintyModel>
          <uncertaintyWeight>1.0</uncertaintyWeight>
        </logicTreeBranch>
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
    <logicTreeBranchingLevel branchingLevelID="bl4">
      <logicTreeBranchSet branchSetID="bs4" uncertaintyType="fdhaSecondaryFDModel" applyToBranches="B3">
        <logicTreeBranch branchID="B4">
          <uncertaintyModel>{secondary_fd}</uncertaintyModel>
          <uncertaintyWeight>1.0</uncertaintyWeight>
        </logicTreeBranch>
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
  </logicTree>
</nrml>
"""
    )


def _write_smlt_xml(path: Path, source_model_file: str = "source_model.xml") -> None:
    path.write_text(
        f"""<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="single_source_model">
    <logicTreeBranchSet branchSetID="bs_source_model"
                        uncertaintyType="sourceModel">
      <logicTreeBranch branchID="SM0">
        <uncertaintyModel>{source_model_file}</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
  </logicTree>
</nrml>
"""
    )


def test_single_branch_equivalence(tmp_path):
    demo_dir = Path(__file__).resolve().parents[3] / "demo" / "hazard_curve"
    source_model = tmp_path / "source_model.xml"
    shutil.copyfile(demo_dir / "source_model.xml", source_model)

    lt_xml = tmp_path / "lt.xml"
    _write_lt_xml(
        lt_xml,
        primary_sr="[Youngs2003PrimarySR]\nstyle = 'all'",
        primary_fd="[Youngs2003PrimaryFD]\nstyle = 'normal'\nscaling_model = 'WC1994'\nnorm_disp_type = 'AD'",
        secondary_sr="[Youngs2003SecondarySR]\nstyle = 'all'",
        secondary_fd="[Youngs2003SecondaryFD]\nstyle = 'all'",
    )
    smlt_xml = tmp_path / "source_model_logic_tree.xml"
    _write_smlt_xml(smlt_xml, source_model.name)

    ini = tmp_path / "job.ini"
    ini.write_text(
        f"""[general]\ndescription = lt_single_branch\n\n[geometry]\nsites = 16.16573727 39.64704451\n\n[site_params]\nreference_vs30_value = 760.0\n\n[erf]\nrupture_mesh_spacing = 1.0\nwidth_of_mfd_bin = 0.1\n\n[calculation]\nsource_model_logic_tree_file = {smlt_xml}\ndisplacement_measure_levels = {{\"FD\": [0.0001, 0.001, 0.01]}}\nfdha_logic_tree_file = {lt_xml}\n"""
    )

    # Logic-tree path
    lt = FdhaLogicTree.from_ini(str(ini))
    res_lt = lt.run(outdir=tmp_path / "out")
    lt_rates = np.array(res_lt.mean_rates, dtype=float)

    branch_rates = np.loadtxt(
        tmp_path / "out" / "hazard_curves" / "branch_0000.csv",
        delimiter=",",
        skiprows=1,
        usecols=1,
    )
    assert np.allclose(lt_rates[0], branch_rates, atol=1e-12, rtol=0)


def test_visini_single_branch_full_logic_tree_equivalence(tmp_path, monkeypatch):
    """Visini2025 single-branch canonical LT has one branch and a stable mean."""
    demo_dir = Path(__file__).resolve().parents[3] / "demo" / "hazard_curve"
    source_model = tmp_path / "source_model.xml"
    shutil.copyfile(demo_dir / "source_model.xml", source_model)

    # Visini's along-strike probability uses Monte Carlo.  Fix the generator
    # so the LT and legacy paths see identical stochastic samples.
    monkeypatch.setattr(
        np.random,
        "default_rng",
        lambda *args, **kwargs: np.random.Generator(np.random.PCG64(12345)),
    )

    fdha_lt_xml = tmp_path / "fdha_lt.xml"
    _write_lt_xml(
        fdha_lt_xml,
        primary_sr="[Youngs2003PrimarySR]\nstyle = 'all'",
        primary_fd=(
            "[Youngs2003PrimaryFD]\n"
            "style = 'all'\n"
            "scaling_model = 'WC1994'\n"
            "norm_disp_type = 'AD'"
        ),
        secondary_sr=(
            "[Visini2025SecondarySR]\n"
            "style = 'reverse'\n"
            "pixel_size = 100\n"
            "along_strike_width = 500\n"
            "distribution_type = 'uniform'\n"
            "segment_sampling = 'truncated'"
        ),
        secondary_fd=(
            "[Visini2025SecondaryFD]\n"
            "style = 'reverse'\n"
            "scaling_model = 'WC1994'"
        ),
    )
    smlt_xml = tmp_path / "source_model_logic_tree.xml"
    _write_smlt_xml(smlt_xml)

    common = f"""[general]
description = visini_single_branch_equivalence

[geometry]
sites = 16.185 39.64

[site_params]
reference_vs30_value = 600.0

[erf]
rupture_mesh_spacing = 1.0
width_of_mfd_bin = 0.1

[calculation]
source_model_logic_tree_file = {smlt_xml.name}
displacement_measure_levels = {{"FD": [0.0001, 0.001, 0.01, 0.1]}}
case = case3
near_far_threshold_km = 0.2
r_threshold_km = 0.1
"""

    lt_ini = tmp_path / "lt.ini"
    lt_ini.write_text(
        common
        + f"""fdha_logic_tree_file = {fdha_lt_xml.name}
"""
    )

    lt = FdhaLogicTree.from_ini(str(lt_ini))
    res_lt = lt.run(outdir=tmp_path / "out_lt")
    lt_rates = np.array(res_lt.mean_rates, dtype=float)

    manifest = json.loads((tmp_path / "out_lt" / "manifest.json").read_text())
    branch_rates = np.loadtxt(
        tmp_path / "out_lt" / manifest["branches"][0]["curve_file"],
        delimiter=",",
        skiprows=1,
        usecols=1,
    )

    assert np.any(lt_rates > 0), "Visini LT result is zero; test is vacuous"
    np.testing.assert_allclose(lt_rates[0], branch_rates, atol=1e-12, rtol=0)

