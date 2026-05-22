from pathlib import Path

from openquake.fdha.logic_tree.nrml_reader import parse


def _repo_root() -> Path:
    p = Path(__file__).resolve()
    for parent in p.parents:
        if (parent / "examples").exists():
            return parent
    raise RuntimeError("Cannot locate repo root (missing examples/)")


def test_parse_fdha_logic_tree_ss_example():
    xml = (
        _repo_root()
        / "openquake/fdha/test/fixtures/examples_archive/logic_tree"
        / "CharacteristicFaultSourceCase2ClassicalPSHA/fdha_logic_tree_SS.xml"
    )
    spec = parse(xml)
    assert spec.logic_tree_id == "pfdha"
    assert len(spec.branching_levels) >= 1
    bs0 = spec.branching_levels[0].branch_sets[0]
    assert bs0.uncertainty_type == "fdhaPrimarySRModel"
    assert bs0.apply_to_style == "strike-slip"
