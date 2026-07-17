# -*- coding: utf-8 -*-
"""
C4 model-contract logic-tree guards:

- FDLT-013 (error): aggregate-definition primary FD model + non-empty
  secondary slot in the same end-branch chain (double counting);
- FDLT-014 (error): mixed DISPLACEMENT_DEFINITION within one FD branch set
  (no cross-definition conversion, Sarmiento et al. 2025 Table 1);
- FDLT-105 (warning): mixed DISPLACEMENT_COMPONENT within one FD branch set;
- FDLT-015 (error): wrong-class output_type -- the class choice IS the
  definition (Lavrentiadis2023PrimaryFD is aggregate-only; the _principal
  variant class pins output_type and accepts no explicit value).
"""
import pytest

from openquake.fdha.logic_tree.types import (
    Branch,
    BranchingLevel,
    BranchSet,
    EndBranch,
    LogicTreeSpec,
    ModelChoice,
)
from openquake.fdha.logic_tree.validators import (
    validate_end_branch_chains,
    validate_spec,
)

pytestmark = pytest.mark.unit


# --------------------------------------------------------------- helpers
def _spec(*branch_sets):
    levels = tuple(
        BranchingLevel(branching_level_id=f"bl_{i}", branch_sets=(bs,))
        for i, bs in enumerate(branch_sets)
    )
    return LogicTreeSpec(
        logic_tree_id="lt_test", branching_levels=levels, basepath=".")


def _fd_set(utype, *models, bs_id="bs_fd"):
    n = len(models)
    weight = repr(1.0 / n)
    branches = tuple(
        Branch(branch_id=f"B{i}", uncertainty_model=m,
               uncertainty_weight=weight)
        for i, m in enumerate(models)
    )
    return BranchSet(branch_set_id=bs_id, uncertainty_type=utype,
                     branches=branches)


def _codes(report):
    return [i.code for i in report.issues]


def _choice(class_name, params=None, branch_id="Bx"):
    return ModelChoice(class_name=class_name, params=params or {},
                       branch_id=branch_id, weight=1.0)


def _end_branch(**selections):
    return EndBranch(source_id="src1", style="normal", weight=1.0,
                     selections=selections)


# --------------------------------------------------------------- FDLT-014
def test_mixed_definition_in_primary_fd_set_is_an_error():
    spec = _spec(_fd_set(
        "fdhaPrimaryFDModel",
        "Youngs2003PrimaryFD",       # principal
        "Kuehn2024PrimaryFD",        # aggregate
    ))
    report = validate_spec(spec)
    errors = [i for i in report.issues if i.code == "FDLT-014"]
    assert len(errors) == 1
    assert errors[0].level == "error"
    assert "principal" in errors[0].message
    assert "aggregate" in errors[0].message


def test_same_definition_set_passes_fdlt014():
    spec = _spec(_fd_set(
        "fdhaPrimaryFDModel",
        "Petersen2011PrimaryFD_bilinear",
        "Petersen2011PrimaryFD_elliptical",
        "Petersen2011PrimaryFD_quadratic",
    ))
    report = validate_spec(spec)
    assert "FDLT-014" not in _codes(report)
    assert "FDLT-105" not in _codes(report)  # all lateral


def test_lavrentiadis_principal_class_groups_with_sum_of_principal():
    """The class choice IS the definition: the _principal variant class
    shares the sum-of-principal definition with Chiou2025, while the
    aggregate parent class does not."""
    ok = _spec(_fd_set(
        "fdhaPrimaryFDModel",
        "[Lavrentiadis2023PrimaryFD_principal]\ninclude_zero_slip = True",
        "Chiou2025PrimaryFD",
    ))
    report = validate_spec(ok)
    assert "FDLT-014" not in _codes(report)
    assert "FDLT-015" not in _codes(report)

    bad = _spec(_fd_set(
        "fdhaPrimaryFDModel",
        "Lavrentiadis2023PrimaryFD",  # aggregate
        "Chiou2025PrimaryFD",         # sum-of-principal
    ))
    assert "FDLT-014" in _codes(validate_spec(bad))


def test_mixed_definition_checked_in_secondary_fd_sets_too():
    """All registered secondary FD models are distributed-definition, so a
    homogeneous secondary set passes."""
    spec = _spec(_fd_set(
        "fdhaSecondaryFDModel",
        "Petersen2011SecondaryFD",
        "Youngs2003SecondaryFD",
    ))
    assert "FDLT-014" not in _codes(validate_spec(spec))


# --------------------------------------------------------------- FDLT-105
def test_mixed_component_in_fd_set_is_a_warning():
    spec = _spec(_fd_set(
        "fdhaPrimaryFDModel",
        "Youngs2003PrimaryFD",   # vertical (principal)
        "Takao2013PrimaryFD",    # net (principal)
    ))
    report = validate_spec(spec)
    assert "FDLT-014" not in _codes(report)  # both principal
    warns = [i for i in report.issues if i.code == "FDLT-105"]
    assert len(warns) == 1
    assert warns[0].level == "warning"
    assert "vertical" in warns[0].message
    assert "net" in warns[0].message


def test_component_check_spans_secondary_sets():
    spec = _spec(_fd_set(
        "fdhaSecondaryFDModel",
        "Petersen2011SecondaryFD",  # lateral
        "Youngs2003SecondaryFD",    # vertical
    ))
    assert "FDLT-105" in _codes(validate_spec(spec))


def test_component_check_ignores_non_fd_sets():
    spec = _spec(_fd_set(
        "fdhaSecondarySRModel",
        "[FixedSecondarySR]\nvalue = 0.0",
        "Petersen2011SecondarySR_default",
    ))
    report = validate_spec(spec)
    assert "FDLT-014" not in _codes(report)
    assert "FDLT-105" not in _codes(report)


# --------------------------------------------------------------- FDLT-013
def test_aggregate_primary_with_secondary_slot_is_an_error():
    eb = _end_branch(
        primary_surf_displ=_choice("Kuehn2024PrimaryFD", branch_id="B2"),
        secondary_surf_rup=_choice("FixedSecondarySR", {"value": 0.0},
                                   branch_id="B3"),
        secondary_surf_displ=_choice("Youngs2003SecondaryFD",
                                     branch_id="B4"),
    )
    report = validate_end_branch_chains([eb])
    errors = [i for i in report.issues if i.code == "FDLT-013"]
    assert len(errors) == 1
    assert errors[0].level == "error"
    assert "Kuehn2024PrimaryFD" in errors[0].message
    assert "double counts" in errors[0].message


def test_aggregate_primary_with_only_secondary_sr_is_still_an_error():
    eb = _end_branch(
        primary_surf_displ=_choice("Kuehn2024PrimaryFD"),
        secondary_surf_rup=_choice("FixedSecondarySR", {"value": 0.0}),
    )
    assert "FDLT-013" in _codes(validate_end_branch_chains([eb]))


def test_aggregate_primary_without_secondary_passes():
    eb = _end_branch(
        primary_surf_rup=_choice("FixedPrimarySR", {"value": 1.0}),
        primary_surf_displ=_choice("Kuehn2024PrimaryFD"),
    )
    assert _codes(validate_end_branch_chains([eb])) == []


def test_principal_primary_with_secondary_passes():
    eb = _end_branch(
        primary_surf_displ=_choice("Youngs2003PrimaryFD"),
        secondary_surf_rup=_choice("Youngs2003SecondarySR"),
        secondary_surf_displ=_choice("Youngs2003SecondaryFD"),
    )
    assert _codes(validate_end_branch_chains([eb])) == []


def test_lavrentiadis_principal_class_with_secondary_passes():
    """The IAEA L23 chains (now the Lavrentiadis2023PrimaryFD_principal
    class) are legitimate sum-of-principal chains and may carry secondary
    slots."""
    eb = _end_branch(
        primary_surf_displ=_choice(
            "Lavrentiadis2023PrimaryFD_principal",
            {"include_zero_slip": True}),
        secondary_surf_rup=_choice("FixedSecondarySR", {"value": 0.0}),
        secondary_surf_displ=_choice("Youngs2003SecondaryFD"),
    )
    assert _codes(validate_end_branch_chains([eb])) == []


def test_lavrentiadis_aggregate_variant_with_secondary_fails():
    eb = _end_branch(
        primary_surf_displ=_choice("Lavrentiadis2023PrimaryFD"),
        secondary_surf_displ=_choice("Youngs2003SecondaryFD"),
    )
    assert "FDLT-013" in _codes(validate_end_branch_chains([eb]))


def test_identical_chains_are_reported_once():
    eb1 = _end_branch(
        primary_surf_displ=_choice("Kuehn2024PrimaryFD", branch_id="B2"),
        secondary_surf_displ=_choice("Youngs2003SecondaryFD",
                                     branch_id="B4"),
    )
    eb2 = EndBranch(source_id="src2", style="reverse", weight=1.0,
                    selections=eb1.selections)
    report = validate_end_branch_chains([eb1, eb2])
    assert _codes(report).count("FDLT-013") == 1


# --------------------------------------------------------------- FDLT-015
def test_prnc_output_type_on_aggregate_class_is_an_error():
    """The class choice IS the definition: requesting disp_prnc_prime on
    the aggregate class must fail at validation time, pointing to the
    _principal variant class."""
    spec = _spec(_fd_set(
        "fdhaPrimaryFDModel",
        "[Lavrentiadis2023PrimaryFD]\noutput_type = disp_prnc_prime\n"
        "include_zero_slip = True",
    ))
    report = validate_spec(spec)
    errors = [i for i in report.issues if i.code == "FDLT-015"]
    assert len(errors) == 1
    assert errors[0].level == "error"
    assert "Lavrentiadis2023PrimaryFD_principal" in errors[0].message


@pytest.mark.parametrize("output_type", [
    "disp_prnc_prime", "disp_agg_prime"])
def test_explicit_output_type_on_principal_class_is_an_error(output_type):
    """The _principal class pins output_type; any explicit value (even the
    redundant disp_prnc_prime) is rejected to keep configurations
    canonical."""
    spec = _spec(_fd_set(
        "fdhaPrimaryFDModel",
        f"[Lavrentiadis2023PrimaryFD_principal]\noutput_type = {output_type}",
    ))
    report = validate_spec(spec)
    errors = [i for i in report.issues if i.code == "FDLT-015"]
    assert len(errors) == 1
    assert "fixed by the class choice" in errors[0].message


def test_aggregate_output_types_on_aggregate_class_pass_fdlt015():
    spec = _spec(_fd_set(
        "fdhaPrimaryFDModel",
        "[Lavrentiadis2023PrimaryFD]\noutput_type = disp_agg_seg",
        "Lavrentiadis2023PrimaryFD",
    ))
    assert "FDLT-015" not in _codes(validate_spec(spec))


def test_clean_principal_class_passes_fdlt015():
    spec = _spec(_fd_set(
        "fdhaPrimaryFDModel",
        "[Lavrentiadis2023PrimaryFD_principal]\ninclude_zero_slip = True",
    ))
    assert "FDLT-015" not in _codes(validate_spec(spec))
