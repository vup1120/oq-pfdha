#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Regenerate paper / reference comparison figures for bundled benchmark suites.

Executes (same code paths as standalone plot scripts):
  1. Valentini et al. - Kumamoto case 2 (Chiou 2025) vs REF_EXCEED
  2. Norcia sensitivity - Youngs2003 AD 85 vs paper rates
  3. Visini et al. - Fig.13 cases 1–3 vs CSV reference data

Outputs go to ``openquake/fdha/test/figures/``.

Usage (from repo root)::

    PYTHONPATH=. python openquake/fdha/test/benchmark/run_all_benchmark_paper_comparison_plots.py
"""

from __future__ import annotations

import importlib.util
from pathlib import Path


def _load_and_run(py_file: Path) -> None:
    name = py_file.stem
    spec = importlib.util.spec_from_file_location(name, py_file)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {py_file}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.main()


def main() -> None:
    here = Path(__file__).resolve().parent
    steps = (
        here / "valentini_et_al_2025" / "plot_kumamoto_logic_tree_vs_paper.py",
        here / "norcia_sensitivity_youngs2003" / "plot_norcia_logic_tree_vs_paper.py",
        here / "visini_et_al_2025" / "fig13" / "plot_fig13_cases_vs_reference.py",
    )
    for pyf in steps:
        print("=" * 72)
        print(pyf.relative_to(here.parent.parent.parent.parent))
        print("=" * 72)
        if not pyf.is_file():
            raise FileNotFoundError(str(pyf))
        _load_and_run(pyf)
    print("All benchmark comparison plots finished.")


if __name__ == "__main__":
    main()
