"""Pedagogical guard: running the example `job.ini` must halt at FDLT-007.

The template `job.ini` loads multi-style `fdha_logic_tree.xml` with normal-style
branches whose weights are ``PLACEHOLDER_WEIGHT``. The validator must reject
that (FDLT-007) and no calculator must execute before the error.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from openquake.fdha.logic_tree.driver import FdhaLogicTree
from openquake.fdha.logic_tree.types import LogicTreeValidationError


def _repo_root() -> Path:
    for p in [Path(__file__).resolve(), *Path(__file__).resolve().parents]:
        if (p / "examples").is_dir() and (p / "openquake").is_dir():
            return p
    raise RuntimeError("Could not locate repository root (missing examples/ + openquake/).")


def test_fdha_ini_template_halts_at_fdlt007(monkeypatch):
    example_dir = (
        _repo_root()
        / "openquake/fdha/test/fixtures/examples_archive/logic_tree"
        / "CharacteristicFaultSourceCase2ClassicalPSHA"
    )
    ini = example_dir / "job.ini"
    assert ini.is_file(), f"Template ini not found: {ini}"

    # Hard guarantee that no calculator runs: if any end-branch calculation is
    # attempted, fail the test immediately.
    import openquake.fdha.logic_tree.driver as driver_mod

    def _must_not_run(*_args, **_kwargs):  # pragma: no cover - should never fire
        raise AssertionError(
            "Calculator was invoked before FDLT-007 halted the logic tree."
        )

    monkeypatch.setattr(driver_mod, "_run_single", _must_not_run)

    with pytest.raises(LogicTreeValidationError) as excinfo:
        FdhaLogicTree.from_ini(str(ini)).run(outdir=example_dir / "out_template_halt")

    msg = str(excinfo.value)
    assert "FDLT-007" in msg, f"Expected FDLT-007 in error, got: {msg!r}"
    # Confirm the offending document is the NM file (by the PLACEHOLDER wording).
    assert "PLACEHOLDER_WEIGHT" in msg, msg
