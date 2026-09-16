from __future__ import annotations

import json
import shutil
from pathlib import Path

from openquake.fdha.logic_tree.driver import FdhaLogicTree


def test_manifest_completeness(tmp_path):
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

    lt_xml = tmp_path / "lt.xml"
    lt_xml.write_text(
        """<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="pfdha">
    <logicTreeBranchSet branchSetID="bs1" uncertaintyType="fdhaPrimarySRModel">
      <logicTreeBranch branchID="B1">
        <uncertaintyModel>Pizza2023PrimarySR</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
    <logicTreeBranchSet branchSetID="bs2" uncertaintyType="fdhaPrimaryFDModel" applyToBranches="B1">
      <logicTreeBranch branchID="B2">
        <uncertaintyModel>[Youngs2003PrimaryFD]
style = normal
scaling_model = WC1994
norm_disp_type = AD</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
    <logicTreeBranchSet branchSetID="bs3" uncertaintyType="fdhaSecondarySRModel" applyToBranches="B1">
      <logicTreeBranch branchID="B3">
        <uncertaintyModel>[Youngs2003SecondarySR]
style = all</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
    <logicTreeBranchSet branchSetID="bs4" uncertaintyType="fdhaSecondaryFDModel" applyToBranches="B3">
      <logicTreeBranch branchID="B4">
        <uncertaintyModel>[Youngs2003SecondaryFD]
style = all</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
  </logicTree>
</nrml>
"""
    )

    ini = tmp_path / "job.ini"
    ini.write_text(
        f"""[general]\ndescription = manifest\n\n[geometry]\nsites = 16.16573727 39.64704451\n\n[site_params]\nreference_vs30_value = 760.0\n\n[erf]\nrupture_mesh_spacing = 1.0\nwidth_of_mfd_bin = 0.1\n\n[calculation]\nsource_model_logic_tree_file = {smlt_xml}\ndisplacement_measure_levels = {{\"FD\": [0.01]}}\nfdha_logic_tree_file = {lt_xml}\n"""
    )

    outdir = tmp_path / "out"
    FdhaLogicTree.from_ini(str(ini)).run(outdir=outdir)

    manifest = json.loads((outdir / "manifest.json").read_text())
    assert "branches" in manifest and manifest["branches"]
    b0 = manifest["branches"][0]
    assert "combined_branch_weight" in b0
    assert "source_model_logic_tree_file" in manifest
    assert "fdha_branch_id" in b0
    assert "source_model_branch_id" in b0
    assert "curve_file" in b0

