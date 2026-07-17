from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

from openquake.fdha.calc.contexts import classify_style
from openquake.fdha.logic_tree.param_parser import (
    R_SIGMA_KM_KEY,
    parse_r_sigma_model,
    parse_uncertainty_model,
)
from openquake.fdha.logic_tree.types import (
    EndBranch,
    LogicTreeSpec,
    ModelChoice,
    CALC_SLOTS_BY_UTYPE,
    FDHA_SLOTS_BY_UTYPE,
)


@dataclass(frozen=True)
class SourceInfo:
    source_id: str
    rake: float

    @property
    def style(self) -> str:
        return classify_style(self.rake)


def enumerate_end_branches(spec: LogicTreeSpec, sources: Iterable[SourceInfo]) -> list[EndBranch]:
    end_branches: list[EndBranch] = []
    for src in sources:
        partials: list[tuple[dict[str, ModelChoice], set[str], float]] = [({}, set(), 1.0)]

        for level in spec.branching_levels:
            for bs in level.branch_sets:
                next_partials: list[tuple[dict[str, ModelChoice], set[str], float]] = []
                for sel, sel_ids, sel_w in partials:
                    if not _branchset_applies(bs, src, sel_ids):
                        next_partials.append((sel, sel_ids, sel_w))
                        continue

                    slot = FDHA_SLOTS_BY_UTYPE.get(bs.uncertainty_type)
                    calc_slot = CALC_SLOTS_BY_UTYPE.get(bs.uncertainty_type)
                    if slot is None and calc_slot is None:
                        next_partials.append((sel, sel_ids, sel_w))
                        continue

                    for br in bs.branches:
                        if calc_slot is not None:
                            # Calc-param branch (engine-style bare scalar):
                            # Cartesian-combines like any model branch but
                            # materialises into [calculation], not [models.*].
                            slot = calc_slot
                            value = parse_r_sigma_model(br.uncertainty_model)
                            class_name = R_SIGMA_KM_KEY
                            params = {R_SIGMA_KM_KEY: value}
                        else:
                            class_name, params = parse_uncertainty_model(br.uncertainty_model)
                        try:
                            w = float(br.uncertainty_weight)
                        except Exception:
                            w = float("nan")
                        new_sel = dict(sel)
                        new_sel[slot] = ModelChoice(
                            class_name=class_name,
                            params=params,
                            branch_id=br.branch_id,
                            weight=w,
                        )
                        new_ids = set(sel_ids)
                        new_ids.add(br.branch_id)
                        next_partials.append((new_sel, new_ids, sel_w * w))
                partials = next_partials

        for sel, _, w in partials:
            end_branches.append(
                EndBranch(
                    source_id=src.source_id,
                    style=src.style,
                    weight=w,
                    selections=sel,
                )
            )
    return end_branches


def _branchset_applies(bs, src: SourceInfo, chosen_branch_ids: set[str]) -> bool:
    if bs.apply_to_style is not None and bs.apply_to_style != src.style:
        return False
    if bs.apply_to_sources is not None:
        allowed = set(bs.apply_to_sources.split())
        if src.source_id not in allowed:
            return False
    if bs.apply_to_branches is not None:
        required = set(bs.apply_to_branches.split())
        if not (chosen_branch_ids & required):
            return False
    return True

