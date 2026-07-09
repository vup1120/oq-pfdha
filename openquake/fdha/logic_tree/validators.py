from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Optional

from openquake.fdha.logic_tree.param_parser import (
    parse_r_threshold_model,
    parse_uncertainty_model,
)
from openquake.fdha.logic_tree.types import (
    ALLOWED_STYLES,
    FDHA_CALC_PARAM_UTYPES,
    FDHA_UNCERTAINTY_TYPES,
    LogicTreeSpec,
    ValidatorIssue,
    ValidatorReport,
)


_METADATA_YAML = Path(__file__).with_name("model_metadata.yaml")


@lru_cache(maxsize=1)
def _load_model_metadata() -> dict:
    """Load `model_metadata.yaml`.

    Entries MUST carry a primary-source citation. Classes absent from the file
    are treated as UNKNOWN (validators skip them silently).
    """
    try:
        import yaml  # type: ignore
    except Exception:
        return {}
    if not _METADATA_YAML.is_file():
        return {}
    with _METADATA_YAML.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}
    return data if isinstance(data, dict) else {}


def _primary_fd_metadata() -> dict[str, dict]:
    md = _load_model_metadata()
    out = md.get("primary_fd_models") or {}
    return {str(k): (v or {}) for k, v in out.items() if isinstance(v, dict)}


def _secondary_distance_applicability() -> dict[str, dict]:
    """Machine-readable distance-applicability metadata, if any.

    Read from ``model_metadata.yaml`` key
    ``secondary_model_distance_applicability`` (entries must carry a
    primary-source citation, like every entry in that file). As of this
    writing no such metadata exists, so callers fall back to the single
    "applicability-range metadata unavailable" advisory warning. Ranges are
    never hardcoded here.
    """
    md = _load_model_metadata()
    out = md.get("secondary_model_distance_applicability") or {}
    return {str(k): (v or {}) for k, v in out.items() if isinstance(v, dict)}


def validate_spec(spec: LogicTreeSpec, source_ids: Optional[set[str]] = None) -> ValidatorReport:
    issues: list[ValidatorIssue] = []

    declared_branch_ids: set[str] = set()
    pfd_meta = _primary_fd_metadata()
    # fdhaCalcRThreshold bookkeeping: (branch_set_id, parsed branch values)
    r_threshold_branch_sets: list[str] = []
    r_threshold_values: list[float] = []

    for level in spec.branching_levels:
        for bs in level.branch_sets:
            is_calc_param = bs.uncertainty_type in FDHA_CALC_PARAM_UTYPES
            if bs.uncertainty_type == "fdhaCalcRThreshold":
                r_threshold_branch_sets.append(bs.branch_set_id)
            # FDLT-005 uncertaintyType
            if bs.uncertainty_type not in FDHA_UNCERTAINTY_TYPES:
                issues.append(
                    ValidatorIssue(
                        code="FDLT-005",
                        level="error",
                        message=f"Invalid uncertaintyType '{bs.uncertainty_type}' in {bs.branch_set_id}",
                    )
                )

            # FDLT-004 applyToStyle
            if bs.apply_to_style is not None and bs.apply_to_style not in ALLOWED_STYLES:
                issues.append(
                    ValidatorIssue(
                        code="FDLT-004",
                        level="error",
                        message=f"Invalid applyToStyle '{bs.apply_to_style}' in {bs.branch_set_id}",
                    )
                )

            # FDLT-003 applyToSources resolves
            if bs.apply_to_sources is not None and source_ids is not None:
                unresolved = [sid for sid in bs.apply_to_sources.split() if sid not in source_ids]
                if unresolved:
                    issues.append(
                        ValidatorIssue(
                            code="FDLT-003",
                            level="error",
                            message=f"applyToSources has unknown source ids {unresolved} in {bs.branch_set_id}",
                        )
                    )

            # Per-branch checks
            weights: list[float] = []
            # (displacement_definition, source_citation) tuples seen in this branch set
            disp_def_entries: list[tuple[str, str]] = []
            # parsed r_threshold_km values seen in this branch set
            bs_r_threshold_values: list[float] = []

            for br in bs.branches:
                declared_branch_ids.add(br.branch_id)

                # FDLT-007 PLACEHOLDER_WEIGHT rejected
                if br.uncertainty_weight.strip() == "PLACEHOLDER_WEIGHT":
                    issues.append(
                        ValidatorIssue(
                            code="FDLT-007",
                            level="error",
                            message=f"PLACEHOLDER_WEIGHT is not allowed (branch {br.branch_id})",
                        )
                    )
                    continue
                try:
                    w = float(br.uncertainty_weight)
                    weights.append(w)
                except Exception:
                    issues.append(
                        ValidatorIssue(
                            code="FDLT-001",
                            level="error",
                            message=f"Non-numeric weight '{br.uncertainty_weight}' (branch {br.branch_id})",
                        )
                    )

                if is_calc_param:
                    # FDLT-010: type-specific grammar (one line
                    # `r_threshold_km = <positive float>`); the model-class
                    # check (FDLT-006) does not apply to parameter branches.
                    try:
                        value = parse_r_threshold_model(br.uncertainty_model)
                    except ValueError as exc:
                        issues.append(
                            ValidatorIssue(
                                code="FDLT-010",
                                level="error",
                                message=f"{exc} (branch {br.branch_id} in {bs.branch_set_id})",
                            )
                        )
                    else:
                        # FDLT-011: duplicate values within one branch set
                        if value in bs_r_threshold_values:
                            issues.append(
                                ValidatorIssue(
                                    code="FDLT-011",
                                    level="error",
                                    message=(
                                        f"Duplicate r_threshold_km value {value} in "
                                        f"{bs.branch_set_id} (branch {br.branch_id})"
                                    ),
                                )
                            )
                        bs_r_threshold_values.append(value)
                    continue

                # FDLT-006 class resolves
                class_name, _ = parse_uncertainty_model(br.uncertainty_model)
                if not _class_is_registered(class_name):
                    issues.append(
                        ValidatorIssue(
                            code="FDLT-006",
                            level="error",
                            message=f"Unregistered FDHA model class '{class_name}' (branch {br.branch_id})",
                        )
                    )

                # Collect citation-traced metadata for advisory checks
                if bs.uncertainty_type == "fdhaPrimaryFDModel" and class_name in pfd_meta:
                    entry = pfd_meta[class_name]
                    dd = entry.get("displacement_definition")
                    src = entry.get("source", "")
                    if dd:
                        disp_def_entries.append((str(dd), str(src)))

            # FDLT-001 weights sum
            if weights:
                s = sum(weights)
                if abs(s - 1.0) > 1e-6:
                    issues.append(
                        ValidatorIssue(
                            code="FDLT-001",
                            level="error",
                            message=f"Weights in {bs.branch_set_id} sum to {s} (expected 1.0 ± 1e-6)",
                        )
                    )

            r_threshold_values.extend(bs_r_threshold_values)

            # FDLT-101 (advisory): mixed displacement definitions within one Primary FD branch set.
            # Driven by model_metadata.yaml; entries missing from the file are UNKNOWN and skipped.
            unique_defs = {dd for dd, _ in disp_def_entries}
            if bs.uncertainty_type == "fdhaPrimaryFDModel" and len(unique_defs) > 1:
                # Build "DEF (source)" fragments; dedupe while preserving first-seen order.
                seen: dict[str, str] = {}
                for dd, src in disp_def_entries:
                    if dd not in seen:
                        seen[dd] = src
                fragments = [f"{dd} [{seen[dd]}]" for dd in sorted(seen.keys())]
                issues.append(
                    ValidatorIssue(
                        code="FDLT-101",
                        level="warning",
                        message=(
                            f"Mixed displacement definitions in {bs.branch_set_id}: "
                            + "; ".join(fragments)
                        ),
                    )
                )

    # FDLT-012: at most one fdhaCalcRThreshold branch set per (merged) tree
    if len(r_threshold_branch_sets) > 1:
        issues.append(
            ValidatorIssue(
                code="FDLT-012",
                level="error",
                message=(
                    "More than one fdhaCalcRThreshold branch set is not allowed; "
                    f"found {r_threshold_branch_sets}"
                ),
            )
        )

    # FDLT-102/FDLT-103 (advisory): cross-check branch values against
    # machine-readable distance-applicability metadata of the registered
    # secondary/distributed models. No such metadata is currently published
    # in model_metadata.yaml, so the FDLT-102 fallback fires; ranges are
    # never hardcoded here.
    if r_threshold_branch_sets:
        applicability = _secondary_distance_applicability()
        if not applicability:
            issues.append(
                ValidatorIssue(
                    code="FDLT-102",
                    level="warning",
                    message=(
                        "applicability-range metadata unavailable: cannot "
                        "cross-check fdhaCalcRThreshold branch values against "
                        "secondary-model distance-applicability ranges"
                    ),
                )
            )
        else:
            for model, entry in sorted(applicability.items()):
                min_km = entry.get("min_km")
                max_km = entry.get("max_km")
                src = entry.get("source", "")
                for value in r_threshold_values:
                    below = min_km is not None and value < float(min_km)
                    above = max_km is not None and value > float(max_km)
                    if below or above:
                        issues.append(
                            ValidatorIssue(
                                code="FDLT-103",
                                level="warning",
                                message=(
                                    f"r_threshold_km = {value} lies outside the "
                                    f"stated distance-applicability range "
                                    f"[{min_km}, {max_km}] km of {model} [{src}]"
                                ),
                            )
                        )

    # FDLT-002 applyToBranches resolves
    for level in spec.branching_levels:
        for bs in level.branch_sets:
            if bs.apply_to_branches is None:
                continue
            unknown = [bid for bid in bs.apply_to_branches.split() if bid not in declared_branch_ids]
            if unknown:
                issues.append(
                    ValidatorIssue(
                        code="FDLT-002",
                        level="error",
                        message=f"applyToBranches references unknown branchIDs {unknown} in {bs.branch_set_id}",
                    )
                )

    return ValidatorReport(issues=tuple(issues))


def check_r_threshold_conflict(
    spec: LogicTreeSpec,
    base_config: dict,
    ini_path: str,
    logic_tree_files: list[str],
) -> None:
    """Enforce the r_threshold_km conflict rule (no silent precedence).

    A job must choose ONE mechanism: either the scalar
    ``[calculation].r_threshold_km`` in the job INI (MODE A) or a
    ``fdhaCalcRThreshold`` branch set in the FDHA logic tree (MODE B).
    Defining both raises ConfigurationError naming both locations.

    Only the user's job INI is checked (materialised ``branch_configs/``
    INIs legitimately carry the scalar and never pass through here).
    """
    from openquake.fdha.calc.config_loader import ConfigurationError

    bs_ids = [
        bs.branch_set_id
        for level in spec.branching_levels
        for bs in level.branch_sets
        if bs.uncertainty_type == "fdhaCalcRThreshold"
    ]
    if not bs_ids:
        return

    scalar_sections = []
    for section in ("calculation", "parameters"):
        body = base_config.get(section)
        if isinstance(body, dict) and "r_threshold_km" in body:
            scalar_sections.append(f"[{section}]")
    if scalar_sections:
        raise ConfigurationError(
            "r_threshold_km is defined twice: as a scalar in "
            f"{' and '.join(scalar_sections)} of {ini_path} AND as "
            f"fdhaCalcRThreshold branch set(s) {bs_ids} in the FDHA logic "
            f"tree file(s) {list(logic_tree_files)}. These mechanisms are "
            "mutually exclusive with no precedence rule: keep the scalar "
            "(single value) or the branch set (epistemic alternatives), "
            "not both."
        )


def _class_is_registered(class_name: str) -> bool:
    try:
        from openquake.fdha import primary_surf_rup, primary_surf_displ, secondary_surf_rup, secondary_surf_displ

        return (
            hasattr(primary_surf_rup, class_name)
            or hasattr(primary_surf_displ, class_name)
            or hasattr(secondary_surf_rup, class_name)
            or hasattr(secondary_surf_displ, class_name)
        )
    except Exception:
        return False
