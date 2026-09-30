#!/usr/bin/env python
"""Diagnostic: what the P11 chain would give for the Suizenji exercise.

TECDOC-2092 Fig. 15(d) shows a Petersen et al. (2011) curve for the 3rd
Kumamoto distributed exercise, whose source is the Suizenji fault at
Mw 5.8. ``Petersen2011SecondaryFD`` is calibrated for M 6-8 and raises
rather than extrapolating below M6, so ``job_distributed_P11_sens1.ini``
cannot complete and the exercise carries no P11 entry in ``manifest.py``.

That refusal is deliberate and stays. It does, however, leave an open
question that matters for validation: is the framework's P11 *unable* to
reproduce the published curve, or merely *unwilling* to try? This script
answers it, by evaluating the same chain with the magnitude check
bypassed at run time.

The bypass is confined to this process: the module source is recompiled
with the guard removed and only the resulting ``get_prob`` is bound onto
the class. Nothing on disk is modified, and importing the library
normally is unaffected. The result is written to

    computed_supplementary/kumamoto_distributed_P11_sens1_extrapolated/

and plotted by ``plot_report_cases.py`` as a dotted, explicitly labelled
diagnostic curve - never as an asserted benchmark entry.

Usage (repo root)::

    PYTHONPATH=. python openquake/fdha/test/benchmark/IAEA/kumamoto/diagnose_P11_suizenji.py
"""

from __future__ import annotations

import csv
import inspect
import os
import sys
import types
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import numpy as np

HERE = Path(__file__).resolve().parent
IAEA = HERE.parent
OUTDIR = IAEA / "computed_supplementary" / \
    "kumamoto_distributed_P11_sens1_extrapolated"

GUARD = ('raise ValueError("Magnitude must be between 6 and 8 '
         'for strike-slip faults")')


def bind_extrapolating_get_prob():
    """Rebind Petersen2011SecondaryFD.get_prob without the M6-8 guard.

    Returns the original method so the caller can restore it.
    """
    from openquake.fdha.secondary_surf_displ import petersen2011 as mod
    from openquake.fdha.secondary_surf_displ.petersen2011 import (
        Petersen2011SecondaryFD)

    src = inspect.getsource(mod)
    if GUARD not in src:
        raise SystemExit(
            "the M6-8 guard was not found in petersen2011.py; the model has "
            "changed and this diagnostic needs revisiting")
    shadow = types.ModuleType("petersen2011_extrapolated")
    shadow.__dict__["__name__"] = "petersen2011_extrapolated"
    exec(compile(src.replace(GUARD, "pass  # guard bypassed for diagnostic"),
                 "<petersen2011-extrapolated>", "exec"), shadow.__dict__)

    original = Petersen2011SecondaryFD.get_prob
    Petersen2011SecondaryFD.get_prob = \
        shadow.Petersen2011SecondaryFD.get_prob
    return Petersen2011SecondaryFD, original


def main():
    from openquake.fdha.logic_tree.driver import FdhaLogicTree

    cls, original = bind_extrapolating_get_prob()
    try:
        result = FdhaLogicTree.from_ini(
            HERE / "job_distributed_P11_sens1.ini").run(outdir=OUTDIR)
    finally:
        cls.get_prob = original          # always restore

    d0 = np.asarray(result.d0, dtype=float)
    rates = np.asarray(result.mean_rates, dtype=float)
    if rates.ndim == 2:
        rates = rates[0]

    csv_path = OUTDIR / "curve.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["disp_m", "rate"])
        for d, v in zip(d0, rates):
            w.writerow([f"{d:g}", f"{v:.6e}"])
    print(f"wrote {csv_path}")

    ref_rows = [r for r in csv.DictReader(
        (IAEA / "reference" /
         "fig15d_kumamoto_distributed_suizenji.csv").open())
        if r["P11"].strip()]
    rd = np.array([float(r["disp_cm"]) for r in ref_rows]) / 100.0
    rv = np.array([float(r["P11"]) for r in ref_rows])
    ci = 10.0 ** np.interp(np.log10(rd), np.log10(d0),
                           np.log10(np.maximum(rates, 1e-300)))
    rel = np.abs(ci / rv - 1.0)

    print("\n  disp (m)   extrapolated   published      ratio")
    for a, b, c in zip(rd, ci, rv):
        print(f"  {a:8.3f}   {b:.3e}      {c:.3e}   {b / c:6.2f}")
    print(f"\n  max |rel err| = {rel.max():.1%}   "
          f"median ratio = {np.median(ci / rv):.2f}")
    print("\nDiagnostic only: the shipped model still refuses M < 6, and this "
          "curve is not a benchmark assertion.")


if __name__ == "__main__":
    sys.exit(main())
