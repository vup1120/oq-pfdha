"""Smoke tests for the engine-backed public displacement examples."""
from pathlib import Path

from openquake.hazardlib.pfd_lt import PFDLogicTree
from openquake.pfd.registry import get_available


ROOT = Path(__file__).parents[4]
EXAMPLES = ROOT / "examples"


def test_public_pfd_logic_trees_are_engine_compatible():
    """The public PFD trees parse with the engine's TOML-backed reader."""
    for name in ("hazard_curve_minimal", "hazard_map_minimal"):
        tree = PFDLogicTree(
            str(EXAMPLES / f"{name}_fdha_logic_tree.xml"))
        assert tree.branchsets


def test_public_jobs_select_engine_displacement():
    """The public jobs use the engine calculation surface."""
    for name in ("hazard_curve_minimal", "hazard_map_minimal"):
        ini = (EXAMPLES / f"{name}.ini").read_text()
        assert "calculation_mode = displacement" in ini
        assert "pfd_logic_tree_file" in ini
        assert '"Disp"' in ini


def test_engine_registry_contains_public_models():
    """All models used by the examples resolve from the engine registry."""
    assert "Youngs2003PrimarySR" in get_available("primary_surf_rup")
    assert "Youngs2003PrimaryFD" in get_available("primary_surf_displ")
    assert "Youngs2003SecondarySR" in get_available("secondary_surf_rup")
    assert "Youngs2003SecondaryFD" in get_available("secondary_surf_displ")
