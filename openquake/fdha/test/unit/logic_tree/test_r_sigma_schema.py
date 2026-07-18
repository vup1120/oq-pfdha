"""Schema/parser/validator tests for the ``fdhaCalcRSigma`` branch set.

XML representation, grammar (zero allowed!), the per-source coverage rule
(FDLT-012), the unit-slip advisory (FDLT-104), the scalar-XOR-branch-set
conflict rule, and the enumerator pseudo-slot. Weight-sum validation reuses
the existing FDLT-001 tolerance.
"""
from __future__ import annotations

import pytest

from openquake.fdha.calc.config_loader import ConfigurationError
from openquake.fdha.logic_tree.enumerator import SourceInfo, enumerate_end_branches
from openquake.fdha.logic_tree.nrml_reader import parse
from openquake.fdha.logic_tree.param_parser import parse_r_sigma_model
from openquake.fdha.logic_tree.validators import (
    check_r_sigma_conflict,
    validate_spec,
)


def _xml_with_branches(branches: str) -> str:
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="lt_rs">
    <logicTreeBranchingLevel branchingLevelID="bl_rs">
      <logicTreeBranchSet branchSetID="bs_rs" uncertaintyType="fdhaCalcRSigma">
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


def _two_set_xml(scope1: str | None, scope2: str | None) -> str:
    """Two fdhaCalcRSigma branch sets with optional applyToSources scopes."""
    a1 = f' applyToSources="{scope1}"' if scope1 is not None else ""
    a2 = f' applyToSources="{scope2}"' if scope2 is not None else ""
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="lt_rs">
    <logicTreeBranchingLevel branchingLevelID="bl_rs1">
      <logicTreeBranchSet branchSetID="bs_rs1" uncertaintyType="fdhaCalcRSigma"{a1}>
{_branch("RS_A", "0.027", "1.0")}
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
    <logicTreeBranchingLevel branchingLevelID="bl_rs2">
      <logicTreeBranchSet branchSetID="bs_rs2" uncertaintyType="fdhaCalcRSigma"{a2}>
{_branch("RS_B", "0.116", "1.0")}
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
  </logicTree>
</nrml>
"""


def _write_spec(tmp_path, xml_text: str):
    path = tmp_path / "rs_lt.xml"
    path.write_text(xml_text)
    return parse(path)


TWO_BRANCH_XML = _xml_with_branches(
    _branch("RS_A", "0.0", "0.4") + _branch("RS_B", "0.065", "0.6")
)


# ---------------------------------------------------------------- grammar
# OpenQuake-engine convention for scalar uncertainty types: the
# <uncertaintyModel> text is a single bare float (here: km, >= 0).

def test_parse_r_sigma_model_valid():
    assert parse_r_sigma_model("0.027") == 0.027
    assert parse_r_sigma_model("  0.116  ") == 0.116
    assert parse_r_sigma_model("\n  0.05\n") == 0.05


def test_parse_r_sigma_model_zero_is_legal():
    """0 selects the boxcar W_p path (perfectly located trace) and must
    parse - a tree can weigh 'trust the trace' against Gaussian classes."""
    assert parse_r_sigma_model("0") == 0.0
    assert parse_r_sigma_model("0.0") == 0.0


@pytest.mark.parametrize(
    "text,fragment",
    [
        ("", "expected single non-negative float"),
        ("0.1 0.2", "expected single non-negative float"),
        ("r_sigma_km = 0.1", "expected single non-negative float"),
        ("abc", "expected single non-negative float"),
        ("-0.1", "non-negative finite float"),
        ("inf", "non-negative finite float"),
        ("nan", "non-negative finite float"),
    ],
)
def test_parse_r_sigma_model_invalid(text, fragment):
    with pytest.raises(ValueError, match=fragment):
        parse_r_sigma_model(text)


# --------------------------------------------------------------- validator

def test_valid_branch_set_parses_to_expected_branch_list(tmp_path):
    spec = _write_spec(tmp_path, TWO_BRANCH_XML)
    report = validate_spec(spec)
    assert report.errors == ()

    (bs,) = [bs for lvl in spec.branching_levels for bs in lvl.branch_sets]
    assert bs.uncertainty_type == "fdhaCalcRSigma"
    parsed = [
        (br.branch_id, parse_r_sigma_model(br.uncertainty_model),
         float(br.uncertainty_weight))
        for br in bs.branches
    ]
    assert parsed == [("RS_A", 0.0, 0.4), ("RS_B", 0.065, 0.6)]

    # Registered type: neither FDLT-005 (unknown type) nor FDLT-006
    # (model-class check) may fire for a calc-param branch set.
    assert not any(i.code in ("FDLT-005", "FDLT-006") for i in report.issues)


def test_weights_must_sum_to_one_reuses_fdlt001(tmp_path):
    xml = _xml_with_branches(
        _branch("RS_A", "0.0", "0.4") + _branch("RS_B", "0.065", "0.4")
    )
    spec = _write_spec(tmp_path, xml)
    report = validate_spec(spec)
    (issue,) = [i for i in report.errors if i.code == "FDLT-001"]
    assert "sum to 0.8" in issue.message
    assert "bs_rs" in issue.message


@pytest.mark.parametrize(
    "model,fragment",
    [
        ("-1", "non-negative finite float"),
        ("abc", "expected single non-negative float"),
        ("r_sigma_km = 0.1", "expected single non-negative float"),
    ],
)
def test_malformed_value_raises_fdlt010(tmp_path, model, fragment):
    xml = _xml_with_branches(_branch("RS_BAD", model, "1.0"))
    spec = _write_spec(tmp_path, xml)
    report = validate_spec(spec)
    (issue,) = [i for i in report.errors if i.code == "FDLT-010"]
    assert fragment in issue.message
    assert "RS_BAD" in issue.message


def test_duplicate_values_raise_fdlt011(tmp_path):
    xml = _xml_with_branches(
        _branch("RS_A", "0.065", "0.5") + _branch("RS_B", "0.0650", "0.5")
    )
    spec = _write_spec(tmp_path, xml)
    report = validate_spec(spec)
    (issue,) = [i for i in report.errors if i.code == "FDLT-011"]
    assert "Duplicate r_sigma_km value 0.065" in issue.message


def test_suspiciously_large_sigma_warns_fdlt104(tmp_path):
    """sigma = 65 km is a metres-typed-as-km slip: warning, not error."""
    xml = _xml_with_branches(_branch("RS_BIG", "65", "1.0"))
    spec = _write_spec(tmp_path, xml)
    report = validate_spec(spec)
    assert report.errors == ()
    (issue,) = [i for i in report.warnings if i.code == "FDLT-104"]
    assert "check the unit" in issue.message
    assert "RS_BIG" in issue.message


def test_petersen_class_values_do_not_warn(tmp_path):
    """All five published two-sided classes must validate silently."""
    xml = _xml_with_branches(
        _branch("RS_1", "0.0269", "0.2") + _branch("RS_2", "0.0438", "0.2")
        + _branch("RS_3", "0.0655", "0.2") + _branch("RS_4", "0.0727", "0.2")
        + _branch("RS_5", "0.116", "0.2")
    )
    spec = _write_spec(tmp_path, xml)
    report = validate_spec(spec)
    assert report.errors == ()
    assert not [i for i in report.warnings if i.code == "FDLT-104"]


# ---------------------------------------- FDLT-012: per-source coverage rule

def test_two_sets_disjoint_sources_are_valid(tmp_path):
    """Correlation groups: disjoint applyToSources scopes may coexist."""
    spec = _write_spec(tmp_path, _two_set_xml("srcA srcB", "srcC"))
    report = validate_spec(spec, source_ids={"srcA", "srcB", "srcC"})
    assert not [i for i in report.errors if i.code == "FDLT-012"]


def test_two_sets_overlapping_sources_raise_fdlt012(tmp_path):
    spec = _write_spec(tmp_path, _two_set_xml("srcA srcB", "srcB srcC"))
    report = validate_spec(spec, source_ids={"srcA", "srcB", "srcC"})
    (issue,) = [i for i in report.errors if i.code == "FDLT-012"]
    assert "srcB" in issue.message
    assert "bs_rs1" in issue.message and "bs_rs2" in issue.message


def test_unscoped_second_set_raises_fdlt012(tmp_path):
    """A set without applyToSources covers every source, so it cannot
    coexist with any other sigma set."""
    spec = _write_spec(tmp_path, _two_set_xml(None, "srcC"))
    report = validate_spec(spec, source_ids={"srcC"})
    (issue,) = [i for i in report.errors if i.code == "FDLT-012"]
    assert "no applyToSources" in issue.message
    assert "at most one sigma branch set" in issue.message


# ------------------------------------------------------------ conflict rule

def test_conflict_rule_scalar_plus_branch_set(tmp_path):
    spec = _write_spec(tmp_path, TWO_BRANCH_XML)
    base_config = {"calculation": {"r_sigma_km": 0.1}}
    with pytest.raises(ConfigurationError) as exc:
        check_r_sigma_conflict(
            spec, base_config, "/path/to/job.ini", ["rs_lt.xml"]
        )
    msg = str(exc.value)
    assert "[calculation]" in msg
    assert "/path/to/job.ini" in msg
    assert "bs_rs" in msg
    assert "rs_lt.xml" in msg
    assert "mutually exclusive" in msg


def test_conflict_rule_detects_parameters_section_scalar(tmp_path):
    spec = _write_spec(tmp_path, TWO_BRANCH_XML)
    base_config = {"parameters": {"r_sigma_km": 0.1}}
    with pytest.raises(ConfigurationError, match=r"\[parameters\]"):
        check_r_sigma_conflict(spec, base_config, "job.ini", ["rs_lt.xml"])


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
    base_config = {"calculation": {"r_sigma_km": 0.1}}
    check_r_sigma_conflict(spec, base_config, "job.ini", ["lt.xml"])  # no raise


def test_no_conflict_without_scalar(tmp_path):
    spec = _write_spec(tmp_path, TWO_BRANCH_XML)
    check_r_sigma_conflict(spec, {"calculation": {}}, "job.ini", ["rs_lt.xml"])


def test_scalar_r_threshold_does_not_conflict_with_sigma_set(tmp_path):
    """r_threshold_km is a fixed calculation parameter (the sigma=0 boxcar
    half-width), NOT the epistemic quantity: it may coexist with a
    fdhaCalcRSigma branch set."""
    spec = _write_spec(tmp_path, TWO_BRANCH_XML)
    base_config = {"calculation": {"r_threshold_km": 0.1}}
    check_r_sigma_conflict(spec, base_config, "job.ini", ["rs_lt.xml"])  # no raise


# --------------------------------------------------------------- enumerator

def test_enumerator_produces_calc_param_end_branches(tmp_path):
    """A fdhaCalcRSigma branch set enters the Cartesian product: one end
    branch per value, carried on the ``calc_r_sigma`` pseudo-slot with
    the branch weight (never silently dropped)."""
    spec = _write_spec(tmp_path, TWO_BRANCH_XML)
    ebs = enumerate_end_branches(spec, [SourceInfo(source_id="s1", rake=-90.0)])
    assert len(ebs) == 2
    got = sorted(
        (
            eb.selections["calc_r_sigma"].branch_id,
            eb.selections["calc_r_sigma"].params["r_sigma_km"],
            eb.weight,
        )
        for eb in ebs
    )
    assert got == [("RS_A", 0.0, 0.4), ("RS_B", 0.065, 0.6)]


def test_enumerator_scopes_sigma_sets_per_source(tmp_path):
    """Disjoint applyToSources scopes route each source through its own
    sigma set only (correlation groups, D5)."""
    spec = _write_spec(tmp_path, _two_set_xml("s1", "s2"))
    ebs = enumerate_end_branches(
        spec,
        [SourceInfo(source_id="s1", rake=-90.0),
         SourceInfo(source_id="s2", rake=-90.0)],
    )
    by_source = {}
    for eb in ebs:
        by_source.setdefault(eb.source_id, []).append(
            eb.selections["calc_r_sigma"].params["r_sigma_km"]
        )
    assert by_source == {"s1": [0.027], "s2": [0.116]}
