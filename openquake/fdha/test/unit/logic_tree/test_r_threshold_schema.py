"""Schema/parser/validator tests for the ``fdhaCalcRThreshold`` branch set.

Phase 2 of the r_threshold-as-epistemic-uncertainty task: XML representation
only (enumeration/execution is Phase 3, and the enumerator must fail loudly
until then). Weight-sum validation reuses the existing FDLT-001 tolerance.
"""
from __future__ import annotations

import pytest

from openquake.fdha.calc.config_loader import ConfigurationError
from openquake.fdha.logic_tree.enumerator import SourceInfo, enumerate_end_branches
from openquake.fdha.logic_tree.nrml_reader import parse
from openquake.fdha.logic_tree.param_parser import parse_r_threshold_model
from openquake.fdha.logic_tree.validators import (
    check_r_threshold_conflict,
    validate_spec,
)


def _xml_with_branches(branches: str) -> str:
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="lt_rt">
    <logicTreeBranchingLevel branchingLevelID="bl_rt">
      <logicTreeBranchSet branchSetID="bs_rt" uncertaintyType="fdhaCalcRThreshold">
{branches}
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
  </logicTree>
</nrml>
"""


def _branch(bid: str, model: str, weight: str) -> str:
    return f"""        <logicTreeBranch branchID="{bid}">
          <uncertaintyModel>{model}</uncertaintyModel>
          <uncertaintyWeight>{weight}</uncertaintyWeight>
        </logicTreeBranch>
"""


def _write_spec(tmp_path, xml_text: str):
    path = tmp_path / "rt_lt.xml"
    path.write_text(xml_text)
    return parse(path)


TWO_BRANCH_XML = _xml_with_branches(
    _branch("RT_A", "0.05", "0.4") + _branch("RT_B", "0.2", "0.6")
)


# ---------------------------------------------------------------- grammar
# OpenQuake-engine convention for scalar uncertainty types: the
# <uncertaintyModel> text is a single bare float (here: km, positive).

def test_parse_r_threshold_model_valid():
    assert parse_r_threshold_model("0.05") == 0.05
    assert parse_r_threshold_model("  2.5  ") == 2.5
    assert parse_r_threshold_model("\n  0.2\n") == 0.2


@pytest.mark.parametrize(
    "text,fragment",
    [
        ("", "expected single positive float"),
        ("0.1 0.2", "expected single positive float"),
        ("r_threshold_km = 0.1", "expected single positive float"),
        ("abc", "expected single positive float"),
        ("-0.1", "positive finite float"),
        ("0", "positive finite float"),
        ("inf", "positive finite float"),
        ("nan", "positive finite float"),
    ],
)
def test_parse_r_threshold_model_invalid(text, fragment):
    with pytest.raises(ValueError, match=fragment):
        parse_r_threshold_model(text)


# --------------------------------------------------------------- validator

def test_valid_branch_set_parses_to_expected_branch_list(tmp_path):
    spec = _write_spec(tmp_path, TWO_BRANCH_XML)
    report = validate_spec(spec)
    assert report.errors == ()

    (bs,) = [bs for lvl in spec.branching_levels for bs in lvl.branch_sets]
    assert bs.uncertainty_type == "fdhaCalcRThreshold"
    parsed = [
        (br.branch_id, parse_r_threshold_model(br.uncertainty_model),
         float(br.uncertainty_weight))
        for br in bs.branches
    ]
    assert parsed == [("RT_A", 0.05, 0.4), ("RT_B", 0.2, 0.6)]

    # Registered type: neither FDLT-005 (unknown type) nor FDLT-006
    # (model-class check) may fire for a calc-param branch set.
    assert not any(i.code in ("FDLT-005", "FDLT-006") for i in report.issues)


def test_applicability_metadata_unavailable_warning_emitted_once(tmp_path):
    spec = _write_spec(tmp_path, TWO_BRANCH_XML)
    report = validate_spec(spec)
    msgs = [i.message for i in report.warnings
            if "applicability-range metadata unavailable" in i.message]
    assert len(msgs) == 1


def test_weights_must_sum_to_one_reuses_fdlt001(tmp_path):
    xml = _xml_with_branches(
        _branch("RT_A", "0.05", "0.4") + _branch("RT_B", "0.2", "0.4")
    )
    spec = _write_spec(tmp_path, xml)
    report = validate_spec(spec)
    (issue,) = [i for i in report.errors if i.code == "FDLT-001"]
    assert "sum to 0.8" in issue.message
    assert "bs_rt" in issue.message


@pytest.mark.parametrize(
    "model,fragment",
    [
        ("-1", "positive finite float"),
        ("abc", "expected single positive float"),
        ("r_threshold_km = 0.1", "expected single positive float"),
    ],
)
def test_malformed_value_raises_fdlt010(tmp_path, model, fragment):
    xml = _xml_with_branches(_branch("RT_BAD", model, "1.0"))
    spec = _write_spec(tmp_path, xml)
    report = validate_spec(spec)
    (issue,) = [i for i in report.errors if i.code == "FDLT-010"]
    assert fragment in issue.message
    assert "RT_BAD" in issue.message


def test_duplicate_values_raise_fdlt011(tmp_path):
    xml = _xml_with_branches(
        _branch("RT_A", "0.2", "0.5") + _branch("RT_B", "0.20", "0.5")
    )
    spec = _write_spec(tmp_path, xml)
    report = validate_spec(spec)
    (issue,) = [i for i in report.errors if i.code == "FDLT-011"]
    assert "Duplicate r_threshold_km value 0.2" in issue.message


def test_multiple_branch_sets_raise_fdlt012(tmp_path):
    xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="lt_rt">
    <logicTreeBranchingLevel branchingLevelID="bl_rt1">
      <logicTreeBranchSet branchSetID="bs_rt1" uncertaintyType="fdhaCalcRThreshold">
{_branch("RT_A", "0.05", "1.0")}
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
    <logicTreeBranchingLevel branchingLevelID="bl_rt2">
      <logicTreeBranchSet branchSetID="bs_rt2" uncertaintyType="fdhaCalcRThreshold">
{_branch("RT_B", "0.2", "1.0")}
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
  </logicTree>
</nrml>
"""
    spec = _write_spec(tmp_path, xml)
    report = validate_spec(spec)
    (issue,) = [i for i in report.errors if i.code == "FDLT-012"]
    assert "More than one fdhaCalcRThreshold branch set" in issue.message
    assert "bs_rt1" in issue.message and "bs_rt2" in issue.message


# ------------------------------------------------------------ conflict rule

def test_conflict_rule_scalar_plus_branch_set(tmp_path):
    spec = _write_spec(tmp_path, TWO_BRANCH_XML)
    base_config = {"calculation": {"r_threshold_km": 0.1}}
    with pytest.raises(ConfigurationError) as exc:
        check_r_threshold_conflict(
            spec, base_config, "/path/to/job.ini", ["rt_lt.xml"]
        )
    msg = str(exc.value)
    assert "[calculation]" in msg
    assert "/path/to/job.ini" in msg
    assert "bs_rt" in msg
    assert "rt_lt.xml" in msg
    assert "mutually exclusive" in msg


def test_conflict_rule_detects_parameters_section_scalar(tmp_path):
    spec = _write_spec(tmp_path, TWO_BRANCH_XML)
    base_config = {"parameters": {"r_threshold_km": 0.1}}
    with pytest.raises(ConfigurationError, match=r"\[parameters\]"):
        check_r_threshold_conflict(spec, base_config, "job.ini", ["rt_lt.xml"])


def test_no_conflict_without_branch_set(tmp_path):
    xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="lt">
    <logicTreeBranchingLevel branchingLevelID="bl1">
      <logicTreeBranchSet branchSetID="bs1" uncertaintyType="fdhaPrimarySRModel">
{_branch("B1", "Youngs2003PrimarySR", "1.0")}
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
  </logicTree>
</nrml>
"""
    spec = _write_spec(tmp_path, xml)
    base_config = {"calculation": {"r_threshold_km": 0.1}}
    check_r_threshold_conflict(spec, base_config, "job.ini", ["lt.xml"])  # no raise


def test_no_conflict_without_scalar(tmp_path):
    spec = _write_spec(tmp_path, TWO_BRANCH_XML)
    check_r_threshold_conflict(spec, {"calculation": {}}, "job.ini", ["rt_lt.xml"])


# --------------------------------------------------------------- enumerator

def test_enumerator_produces_calc_param_end_branches(tmp_path):
    """A fdhaCalcRThreshold branch set enters the Cartesian product: one end
    branch per value, carried on the ``calc_r_threshold`` pseudo-slot with
    the branch weight (never silently dropped)."""
    spec = _write_spec(tmp_path, TWO_BRANCH_XML)
    ebs = enumerate_end_branches(spec, [SourceInfo(source_id="s1", rake=-90.0)])
    assert len(ebs) == 2
    got = sorted(
        (
            eb.selections["calc_r_threshold"].branch_id,
            eb.selections["calc_r_threshold"].params["r_threshold_km"],
            eb.weight,
        )
        for eb in ebs
    )
    assert got == [("RT_A", 0.05, 0.4), ("RT_B", 0.2, 0.6)]
