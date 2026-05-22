from openquake.fdha.logic_tree.enumerator import SourceInfo, enumerate_end_branches
from openquake.fdha.logic_tree.types import Branch, BranchSet, BranchingLevel, LogicTreeSpec


def test_enumerator_cartesian_product_single_source():
    spec = LogicTreeSpec(
        logic_tree_id="pfdha",
        basepath=".",
        branching_levels=(
            BranchingLevel(
                branching_level_id="bl1",
                branch_sets=(
                    BranchSet(
                        branch_set_id="bs1",
                        uncertainty_type="fdhaPrimarySRModel",
                        apply_to_style="strike-slip",
                        branches=(
                            Branch("b1", "Pizza2023PrimarySR", "0.5"),
                            Branch("b2", "Pizza2023PrimarySR", "0.5"),
                        ),
                    ),
                ),
            ),
        ),
    )
    src = SourceInfo(source_id="6", rake=0.0)  # strike-slip
    ebs = enumerate_end_branches(spec, [src])
    assert len(ebs) == 2
    assert all("primary_surf_rup" in eb.selections for eb in ebs)

