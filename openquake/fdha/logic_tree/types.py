from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional


FDHA_UNCERTAINTY_TYPES = {
    "fdhaPrimarySRModel",
    "fdhaPrimaryFDModel",
    "fdhaSecondarySRModel",
    "fdhaSecondaryFDModel",
}

FDHA_SLOTS_BY_UTYPE = {
    "fdhaPrimarySRModel": "primary_surf_rup",
    "fdhaPrimaryFDModel": "primary_surf_displ",
    "fdhaSecondarySRModel": "secondary_surf_rup",
    "fdhaSecondaryFDModel": "secondary_surf_displ",
}

ALLOWED_STYLES = {"strike-slip", "reverse", "normal"}


@dataclass(frozen=True)
class Branch:
    branch_id: str
    uncertainty_model: str  # raw text inside <uncertaintyModel>
    uncertainty_weight: str  # raw text inside <uncertaintyWeight>


@dataclass(frozen=True)
class BranchSet:
    branch_set_id: str
    uncertainty_type: str
    apply_to_sources: Optional[str] = None
    apply_to_branches: Optional[str] = None
    apply_to_style: Optional[str] = None
    branches: tuple[Branch, ...] = ()


@dataclass(frozen=True)
class BranchingLevel:
    branching_level_id: str
    branch_sets: tuple[BranchSet, ...]


@dataclass(frozen=True)
class LogicTreeSpec:
    logic_tree_id: str
    branching_levels: tuple[BranchingLevel, ...]
    basepath: str  # directory of the XML file


@dataclass(frozen=True)
class ModelChoice:
    class_name: str
    params: dict[str, Any]
    branch_id: str
    weight: float


@dataclass(frozen=True)
class EndBranch:
    source_id: str
    style: str
    weight: float
    selections: dict[str, ModelChoice]  # slot -> choice


@dataclass(frozen=True)
class ValidatorIssue:
    code: str
    message: str
    level: str  # 'error' or 'warning'


@dataclass(frozen=True)
class ValidatorReport:
    issues: tuple[ValidatorIssue, ...]

    @property
    def errors(self) -> tuple[ValidatorIssue, ...]:
        return tuple(i for i in self.issues if i.level == "error")

    @property
    def warnings(self) -> tuple[ValidatorIssue, ...]:
        return tuple(i for i in self.issues if i.level == "warning")

    def raise_if_errors(self) -> None:
        if self.errors:
            msgs = "\n".join(f"{i.code}: {i.message}" for i in self.errors)
            raise LogicTreeValidationError(msgs)


class LogicTreeValidationError(Exception):
    pass

