# -*- coding: utf-8 -*-
"""
C4 model contract, driver level:

- an aggregate-definition primary FD chain carrying secondary slots must
  HALT at FDLT-013 before any calculator runs;
- the same chain without secondary slots runs single-bucket: distributed
  rates are exactly zero and the total equals the principal column.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

from openquake.fdha.logic_tree.driver import FdhaLogicTree
from openquake.fdha.logic_tree.types import LogicTreeValidationError

pytestmark = pytest.mark.integration

_SECONDARY_LEVELS = """
    <logicTreeBranchSet branchSetID="bs3" uncertaintyType="fdhaSecondarySRModel" applyToBranches="FD">
      <logicTreeBranch branchID="SSR">
        <uncertaintyModel>[FixedSecondarySR]\nvalue = 0.0</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
    <logicTreeBranchSet branchSetID="bs4" uncertaintyType="fdhaSecondaryFDModel" applyToBranches="SSR">
      <logicTreeBranch branchID="SFD">
        <uncertaintyModel>[Youngs2003SecondaryFD]\nstyle = 'all'</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
"""


def _write_job(tmp_path: Path, with_secondary: bool) -> Path:
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
        f"""<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="pfdha">
    <logicTreeBranchSet branchSetID="bs1" uncertaintyType="fdhaPrimarySRModel">
      <logicTreeBranch branchID="SR1">
        <uncertaintyModel>[FixedPrimarySR]\nvalue = 1.0</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
    <logicTreeBranchSet branchSetID="bs2" uncertaintyType="fdhaPrimaryFDModel" applyToBranches="SR1">
      <logicTreeBranch branchID="FD">
        <uncertaintyModel>[Lavrentiadis2023PrimaryFD_aggregate]</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
{_SECONDARY_LEVELS if with_secondary else ''}
  </logicTree>
</nrml>
"""
    )

    ini = tmp_path / "job.ini"
    ini.write_text(
        # On-trace site (fault-trace vertex, as in the demo job).
        f"""[general]\ndescription = lt_aggregate_contract\n\n[geometry]\nsites = 16.1455213236 39.6196231258\n\n[site_params]\nreference_vs30_value = 760.0\n\n[erf]\nrupture_mesh_spacing = 1.0\nwidth_of_mfd_bin = 0.1\n\n[calculation]\nsource_model_logic_tree_file = {smlt_xml}\ndisplacement_measure_levels = {{\"FD\": [0.01, 0.1]}}\nfdha_logic_tree_file = {lt_xml}\n"""
    )
    return ini


def test_aggregate_chain_with_secondary_halts_at_fdlt013(tmp_path, monkeypatch):
    ini = _write_job(tmp_path, with_secondary=True)

    # No calculator may run before the halt.
    import openquake.fdha.logic_tree.driver as driver_mod

    def _must_not_run(*_a, **_k):  # pragma: no cover - must never fire
        raise AssertionError("Calculator ran before FDLT-013 halted the job")

    monkeypatch.setattr(driver_mod, "_run_single", _must_not_run)

    with pytest.raises(LogicTreeValidationError) as excinfo:
        FdhaLogicTree.from_ini(str(ini)).run(outdir=tmp_path / "out")
    msg = str(excinfo.value)
    assert "FDLT-013" in msg
    assert "Lavrentiadis2023PrimaryFD_aggregate" in msg


def test_aggregate_chain_runs_single_bucket(tmp_path):
    ini = _write_job(tmp_path, with_secondary=False)
    res = FdhaLogicTree.from_ini(str(ini)).run(outdir=tmp_path / "out")

    total = np.array(res.mean_rates, dtype=float)
    principal = np.array(res.mean_rates_principal, dtype=float)
    distributed = np.array(res.mean_rates_distributed, dtype=float)

    assert np.any(total > 0), "aggregate chain produced zero hazard (vacuous)"
    # Single bucket: everything in the principal column, distributed == 0.
    assert np.all(distributed == 0.0)
    np.testing.assert_array_equal(total, principal)
