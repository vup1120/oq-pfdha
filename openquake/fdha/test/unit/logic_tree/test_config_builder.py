from openquake.fdha.calc.config_loader import load_config
from openquake.fdha.logic_tree.config_builder import build_config, dump_config_to_ini
from openquake.fdha.logic_tree.types import EndBranch, ModelChoice


def test_dump_config_to_ini_roundtrip(tmp_path):
    base = {"calculation": {"source_model_file": "x.xml"}}
    eb = EndBranch(
        source_id="6",
        style="strike-slip",
        weight=1.0,
        selections={
            "primary_surf_rup": ModelChoice("Pizza2023PrimarySR", {"style": "all"}, "b1", 1.0),
        },
    )
    cfg = build_config(base, eb)
    text = dump_config_to_ini(cfg)
    branch_dir = tmp_path / "branch_configs"
    branch_dir.mkdir()
    p = branch_dir / "t.ini"
    p.write_text(text)
    parsed = load_config(str(p))
    assert parsed["models"]["primary_surf_rup"]["type"] == "Pizza2023PrimarySR"
    assert parsed["models"]["primary_surf_rup"]["parameters"]["style"] == "all"

