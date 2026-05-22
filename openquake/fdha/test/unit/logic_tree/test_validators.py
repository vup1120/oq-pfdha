from pathlib import Path

from openquake.fdha.logic_tree.nrml_reader import parse
from openquake.fdha.logic_tree.validators import validate_spec


def _repo_root() -> Path:
    p = Path(__file__).resolve()
    for parent in p.parents:
        if (parent / "examples").exists():
            return parent
    raise RuntimeError("Cannot locate repo root (missing examples/)")


def test_validator_rejects_placeholder_weight():
    xml = (
        _repo_root()
        / "openquake/fdha/test/fixtures/examples_archive/logic_tree"
        / "CharacteristicFaultSourceCase2ClassicalPSHA/fdha_logic_tree_NM.xml"
    )
    spec = parse(xml)
    report = validate_spec(spec, source_ids={"6"})
    assert any(i.code == "FDLT-007" for i in report.issues)
