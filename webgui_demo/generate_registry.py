#!/usr/bin/env python
"""Generate registry.json from the actual oq-pfdha model registry.

The effective model registry is the public namespace of the four model
packages (see openquake/fdha/logic_tree/validators.py:184-195,
_class_is_registered). This script imports those packages and snapshots
every exported model class, so the GUI never hard-codes a model list.

Curation (per project owner's decision at the Phase 0 gate):
- utility constant-value classes are excluded: FixedPrimarySR,
  FixedSecondarySR (openquake/fdha/*/fixed.py)
- alias subclasses are excluded: MammarellaEtAl2024PrimarySR
  (alias of Mammarella2024PrimarySR, primary_surf_rup/mammarella2024.py:395),
  Petersen2011SecondarySR_default (alias variant,
  secondary_surf_rup/petersen2011.py:152)

For each class we also record whether its __init__ accepts the per-model
``n_sigma`` truncation parameter and its default value (e.g.
Youngs2003PrimaryFD n_sigma=6, primary_surf_displ/youngs2003.py:53-58;
Takao2013PrimaryFD n_sigma=3, primary_surf_displ/takao2013.py:43-44).

Run from the repository root with the package installed:

    python webgui_demo/generate_registry.py
"""
from __future__ import annotations

import importlib
import inspect
import json
from pathlib import Path

# Slot names follow FDHA_SLOTS_BY_UTYPE in openquake/fdha/logic_tree/types.py:14-19
SLOTS = {
    "fdhaPrimarySRModel": ("primary_surf_rup", "BasePrimarySurfRup"),
    "fdhaPrimaryFDModel": ("primary_surf_displ", "BasePrimarySurfDispl"),
    "fdhaSecondarySRModel": ("secondary_surf_rup", "BaseSecondarySurfRup"),
    "fdhaSecondaryFDModel": ("secondary_surf_displ", "BaseSecondarySurfDispl"),
}

EXCLUDED = {
    "FixedPrimarySR",            # utility constant, primary_surf_rup/fixed.py:27
    "FixedSecondarySR",          # utility constant, secondary_surf_rup/fixed.py:27
    "MammarellaEtAl2024PrimarySR",       # alias, primary_surf_rup/mammarella2024.py:395
    "Petersen2011SecondarySR_default",   # alias, secondary_surf_rup/petersen2011.py:152
}


def snapshot() -> dict:
    out: dict = {"slots": {}, "excluded_by_curation": sorted(EXCLUDED)}
    for utype, (pkg_name, _base) in SLOTS.items():
        pkg = importlib.import_module(f"openquake.fdha.{pkg_name}")
        models = []
        for name in sorted(dir(pkg)):
            obj = getattr(pkg, name)
            if not (inspect.isclass(obj) and name[0].isupper()):
                continue
            if name.startswith("Base") or name in EXCLUDED:
                continue
            entry = {
                "class_name": name,
                "import_path": f"openquake.fdha.{pkg_name}.{obj.__module__.split('.')[-1]}",
                "module": obj.__module__,
            }
            try:
                params = inspect.signature(obj.__init__).parameters
                if "n_sigma" in params:
                    default = params["n_sigma"].default
                    entry["n_sigma_default"] = (
                        None if default is inspect.Parameter.empty else default
                    )
            except (TypeError, ValueError):
                pass
            models.append(entry)
        out["slots"][utype] = {"package": f"openquake.fdha.{pkg_name}", "models": models}
    return out


if __name__ == "__main__":
    data = snapshot()
    dest = Path(__file__).parent / "registry.json"
    dest.write_text(json.dumps(data, indent=2) + "\n")
    n = sum(len(v["models"]) for v in data["slots"].values())
    print(f"wrote {dest} with {n} model classes")
