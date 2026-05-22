from __future__ import annotations

from openquake.fdha.logic_tree.enumerator import SourceInfo, enumerate_end_branches
from openquake.fdha.logic_tree.types import Branch, BranchSet, BranchingLevel, LogicTreeSpec


def test_style_filter_strike_slip_only():
    spec = LogicTreeSpec(
        logic_tree_id="pfdha",
        basepath=".",
        branching_levels=(
            BranchingLevel(
                branching_level_id="bl1",
                branch_sets=(
                    BranchSet(
                        branch_set_id="bs_ss",
                        uncertainty_type="fdhaPrimarySRModel",
                        apply_to_style="strike-slip",
                        branches=(Branch("SS", "Pizza2023PrimarySR", "1.0"),),
                    ),
                    BranchSet(
                        branch_set_id="bs_nm",
                        uncertainty_type="fdhaPrimarySRModel",
                        apply_to_style="normal",
                        branches=(Branch("NM", "Pizza2023PrimarySR", "1.0"),),
                    ),
                ),
            ),
        ),
    )
    src = SourceInfo(source_id="6", rake=0.0)  # strike-slip
    ebs = enumerate_end_branches(spec, [src])
    assert len(ebs) == 1
    assert list(ebs[0].selections.values())[0].branch_id == "SS"

