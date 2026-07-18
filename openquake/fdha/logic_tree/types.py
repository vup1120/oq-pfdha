"""Shared logic-tree data types: branches, branch sets and end branches."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional


FDHA_UNCERTAINTY_TYPES = {
    "fdhaPrimarySRModel",
    "fdhaPrimaryFDModel",
    "fdhaSecondarySRModel",
    "fdhaSecondaryFDModel",
    "fdhaCalcRSigma",
}

FDHA_SLOTS_BY_UTYPE = {
    "fdhaPrimarySRModel": "primary_surf_rup",
    "fdhaPrimaryFDModel": "primary_surf_displ",
    "fdhaSecondarySRModel": "secondary_surf_rup",
    "fdhaSecondaryFDModel": "secondary_surf_displ",
}

# Calculation-parameter uncertainty types: their <uncertaintyModel> carries a
# scalar calculation parameter (epistemic alternative), not an FDHA model
# class. They are validated by a type-specific grammar and must never fall
# through the model-class checks (FDLT-006) or the [models.*] materialisation.
#
# fdhaCalcRSigma: alternative values of [calculation].r_sigma_km, the
# two-sided mapping-accuracy sigma of the rupture-location term fr(r) of
# Petersen et al. (2011, BSSA 101, 805-825, doi:10.1785/0120100035, Tables
# 2-3). 0 selects the boxcar W_p path (half-width r_threshold_km); > 0
# selects the pure-Gaussian path. Treating the mapping-accuracy class as
# weighted logic-tree branches follows Petersen et al. (2011, p. 811): "this
# epistemic uncertainty should be considered as alternative branches in a
# logic tree". Branch sets may be scoped per correlation group of sources via
# applyToSources; each source may be covered by at most one set (FDLT-012).
FDHA_CALC_PARAM_UTYPES = {
    "fdhaCalcRSigma",
}

# Pseudo-slot names used for calc-param uncertainty types inside
# ``EndBranch.selections``. Keeping them in ``selections`` (rather than a
# parallel structure) means fingerprinting, dedup, weight multiplication and
# manifest branch paths treat sigma branches as ordinary realizations.
# ``build_config`` materialises these into the branch INI's ``[calculation]``
# section instead of ``[models.*]``.
CALC_SLOTS_BY_UTYPE = {
    "fdhaCalcRSigma": "calc_r_sigma",
}

CALC_R_SIGMA_SLOT = CALC_SLOTS_BY_UTYPE["fdhaCalcRSigma"]

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

