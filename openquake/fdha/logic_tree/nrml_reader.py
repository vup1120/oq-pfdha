from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import xml.etree.ElementTree as ET

from openquake.fdha.logic_tree.types import (
    Branch,
    BranchSet,
    BranchingLevel,
    LogicTreeSpec,
)


NRML_NS = {"nrml": "http://openquake.org/xmlns/nrml/0.4"}


def _get_attr(elem: ET.Element, name: str) -> Optional[str]:
    v = elem.get(name)
    return v if v is not None else None


def parse(xml_path: str | Path) -> LogicTreeSpec:
    xml_path = Path(xml_path)
    basepath = str(xml_path.parent)
    tree = ET.parse(xml_path)
    root = tree.getroot()
    lt = root.find("nrml:logicTree", NRML_NS)
    if lt is None:
        raise ValueError("Missing <logicTree> element")

    logic_tree_id = lt.get("logicTreeID") or ""
    branching_levels: list[BranchingLevel] = []

    for bl in lt.findall("nrml:logicTreeBranchingLevel", NRML_NS):
        bl_id = bl.get("branchingLevelID") or ""
        branch_sets: list[BranchSet] = []
        for bs in bl.findall("nrml:logicTreeBranchSet", NRML_NS):
            bs_id = bs.get("branchSetID") or ""
            utype = bs.get("uncertaintyType") or ""
            ats = _get_attr(bs, "applyToSources")
            atb = _get_attr(bs, "applyToBranches")
            atstyle = _get_attr(bs, "applyToStyle")
            branches: list[Branch] = []
            for br in bs.findall("nrml:logicTreeBranch", NRML_NS):
                bid = br.get("branchID") or ""
                um = br.find("nrml:uncertaintyModel", NRML_NS)
                uw = br.find("nrml:uncertaintyWeight", NRML_NS)
                if um is None or uw is None:
                    raise ValueError(f"Branch {bid} missing uncertaintyModel/uncertaintyWeight")
                uncertainty_model = (um.text or "").strip("\n")
                uncertainty_weight = (uw.text or "").strip()
                branches.append(
                    Branch(
                        branch_id=bid,
                        uncertainty_model=uncertainty_model,
                        uncertainty_weight=uncertainty_weight,
                    )
                )
            branch_sets.append(
                BranchSet(
                    branch_set_id=bs_id,
                    uncertainty_type=utype,
                    apply_to_sources=ats,
                    apply_to_branches=atb,
                    apply_to_style=atstyle,
                    branches=tuple(branches),
                )
            )
        branching_levels.append(BranchingLevel(branching_level_id=bl_id, branch_sets=tuple(branch_sets)))

    return LogicTreeSpec(
        logic_tree_id=logic_tree_id,
        branching_levels=tuple(branching_levels),
        basepath=basepath,
    )

