from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Optional

from openquake.fdha.logic_tree.param_parser import parse_uncertainty_model
from openquake.fdha.logic_tree.types import (
    ALLOWED_STYLES,
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


def validate_spec(spec: LogicTreeSpec, source_ids: Optional[set[str]] = None) -> ValidatorReport:
    issues: list[ValidatorIssue] = []

    declared_branch_ids: set[str] = set()
    pfd_meta = _primary_fd_metadata()

    for level in spec.branching_levels:
        for bs in level.branch_sets:
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
