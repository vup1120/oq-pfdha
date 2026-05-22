from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np

from openquake.fdha.logic_tree.driver import FdhaLogicTree


def test_two_branch_weighted_sum(tmp_path):
    demo_dir = Path(__file__).resolve().parents[3] / "demo" / "hazard_curve"
    source_model = tmp_path / "source_model.xml"
    shutil.copyfile(demo_dir / "source_model.xml", source_model)
    smlt_xml = tmp_path / "source_model_logic_tree.xml"
    smlt_xml.write_text(
        f"""<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="smlt">
    <logicTreeBranchSet branchSetID="bs1" uncertaintyType="sourceModel">
      <logicTreeBranch branchID="sm0">
        <uncertaintyModel>{source_model.name}</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
  </logicTree>
</nrml>
"""
    )

    w = 0.3

    lt_xml = tmp_path / "lt.xml"
    lt_xml.write_text(
        f"""<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="pfdha">
    <logicTreeBranchingLevel branchingLevelID="bl1">
      <logicTreeBranchSet branchSetID="bs1" uncertaintyType="fdhaPrimarySRModel">
        <logicTreeBranch branchID="SR0">
          <uncertaintyModel>[FixedPrimarySR]\nvalue = 0.0</uncertaintyModel>
          <uncertaintyWeight>{w}</uncertaintyWeight>
        </logicTreeBranch>
        <logicTreeBranch branchID="SR1">
          <uncertaintyModel>[FixedPrimarySR]\nvalue = 1.0</uncertaintyModel>
          <uncertaintyWeight>{1-w}</uncertaintyWeight>
        </logicTreeBranch>
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
    <logicTreeBranchingLevel branchingLevelID="bl2">
      <logicTreeBranchSet branchSetID="bs2" uncertaintyType="fdhaPrimaryFDModel" applyToBranches="SR0 SR1">
        <logicTreeBranch branchID="FD">
          <uncertaintyModel>[Youngs2003PrimaryFD]\nstyle = 'normal'\nscaling_model = 'WC1994'\nnorm_disp_type = 'AD'</uncertaintyModel>
          <uncertaintyWeight>1.0</uncertaintyWeight>
        </logicTreeBranch>
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
    <logicTreeBranchingLevel branchingLevelID="bl3">
      <logicTreeBranchSet branchSetID="bs3" uncertaintyType="fdhaSecondarySRModel" applyToBranches="SR0 SR1">
        <logicTreeBranch branchID="SSR">
          <uncertaintyModel>[FixedSecondarySR]\nvalue = 0.0</uncertaintyModel>
          <uncertaintyWeight>1.0</uncertaintyWeight>
        </logicTreeBranch>
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
    <logicTreeBranchingLevel branchingLevelID="bl4">
      <logicTreeBranchSet branchSetID="bs4" uncertaintyType="fdhaSecondaryFDModel" applyToBranches="SSR">
        <logicTreeBranch branchID="SFD">
          <uncertaintyModel>[Youngs2003SecondaryFD]\nstyle = 'all'</uncertaintyModel>
          <uncertaintyWeight>1.0</uncertaintyWeight>
        </logicTreeBranch>
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
  </logicTree>
</nrml>
"""
    )

    ini = tmp_path / "job.ini"
    ini.write_text(
        f"""[general]\ndescription = lt_two_branch\n\n[geometry]\nsites = 16.16573727 39.64704451\n\n[site_params]\nreference_vs30_value = 760.0\n\n[erf]\nrupture_mesh_spacing = 1.0\nwidth_of_mfd_bin = 0.1\n\n[calculation]\nsource_model_logic_tree_file = {smlt_xml}\ndisplacement_measure_levels = {{\"FD\": [0.01, 0.1]}}\nfdha_logic_tree_file = {lt_xml}\n"""
    )

    lt = FdhaLogicTree.from_ini(str(ini))
    res_lt = lt.run(outdir=tmp_path / "out")
    mean_lt = np.array(res_lt.mean_rates, dtype=float)

    r0 = np.loadtxt(tmp_path / "out" / "hazard_curves" / "branch_0000.csv", delimiter=",", skiprows=1, usecols=1)
    r1 = np.loadtxt(tmp_path / "out" / "hazard_curves" / "branch_0001.csv", delimiter=",", skiprows=1, usecols=1)

    expected = w * r0 + (1 - w) * r1
    assert np.allclose(mean_lt, expected, atol=1e-12, rtol=0)

