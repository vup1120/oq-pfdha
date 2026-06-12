#!/usr/bin/env python
"""Generate registry.json from the actual oq-pfdha model registry.

The effective model registry is the public namespace of the four model
packages (see openquake/fdha/logic_tree/validators.py:184-195,
_class_is_registered). This script imports those packages and snapshots
every exported model class, so the GUI never hard-codes a model list.

It additionally harvests, for each class:

- the per-model parameter table from the User Manual
  (docs/UserManual_Enhanced/models/{primary,secondary}/<Model>.md,
  rows under the "| Name | Type | Units | Default | Allowed | ..." header),
  so the GUI can render authoritative parameter documentation;
- any constructor parameters exposed with explicit defaults
  (e.g. n_sigma: Youngs2003PrimaryFD=6, primary_surf_displ/youngs2003.py:53-58;
  Takao2013PrimaryFD=3, primary_surf_displ/takao2013.py:43-44);
- a parameter prefill string for the GUI editor, taken verbatim from the
  shipped logic trees (Norcia Case 3 benchmark and Taiwan example), i.e.
  values with in-repo provenance.

Curation (per project owner's decision): utility constant classes
(FixedPrimarySR, FixedSecondarySR) and alias subclasses
(MammarellaEtAl2024PrimarySR, Petersen2011SecondarySR_default) are excluded.

Run from the repository root with the package installed:

    python webgui_demo/generate_registry.py
"""
from __future__ import annotations

import importlib
import inspect
import json
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DOCS = REPO / "docs" / "UserManual_Enhanced" / "models"

# Slot names follow FDHA_SLOTS_BY_UTYPE in openquake/fdha/logic_tree/types.py:14-19
SLOTS = {
    "fdhaPrimarySRModel": "primary_surf_rup",
    "fdhaPrimaryFDModel": "primary_surf_displ",
    "fdhaSecondarySRModel": "secondary_surf_rup",
    "fdhaSecondaryFDModel": "secondary_surf_displ",
}

EXCLUDED = {
    "FixedPrimarySR",                  # utility constant, primary_surf_rup/fixed.py:27
    "FixedSecondarySR",                # utility constant, secondary_surf_rup/fixed.py:27
    "MammarellaEtAl2024PrimarySR",     # alias, primary_surf_rup/mammarella2024.py:395
    "Petersen2011SecondarySR_default", # alias, secondary_surf_rup/petersen2011.py:152
}

# Class -> User-Manual page holding its parameter table.
DOC_MAP = {
    "Youngs2003PrimaryFD": "primary/Youngs2003.md",
    "Petersen2011PrimaryFD": "primary/Petersen2011.md",
    "Petersen2011PrimaryFD_bilinear": "primary/Petersen2011.md",
    "Petersen2011PrimaryFD_elliptical": "primary/Petersen2011.md",
    "Petersen2011PrimaryFD_quadratic": "primary/Petersen2011.md",
    "Chiou2025PrimaryFD": "primary/Chiou2025.md",
    "Kuehn2024PrimaryFD": "primary/Kuehn2024.md",
    "Lavrentiadis2023PrimaryFD": "primary/Lavrentiadis2023.md",
    "Moss2022PrimaryFD": "primary/Moss2022.md",
    "Moss2024PrimaryFD": "primary/Moss2024.md",
    "Takao2013PrimaryFD": "primary/Takao2013.md",
    "Mammarella2024PrimarySR": "primary/MammarellaEtAl2024.md",
    "Yang2021PrimarySR": "primary/Yang2021.md",
    "Petersen2011SecondarySR": "secondary/Petersen2011.md",
    "Petersen2011SecondaryFD": "secondary/Petersen2011.md",
    "Visini2025SecondarySR": "secondary/VisiniEtAl2025.md",
    "Visini2025SecondaryFD": "secondary/VisiniEtAl2025.md",
    "Youngs2003SecondarySR": "secondary/Youngs2003.md",
    "Youngs2003SecondaryFD": "secondary/Youngs2003.md",
}

# Editor prefill for <uncertaintyModel> parameter lines. Values are taken
# verbatim from logic trees shipped in this repository:
#   [N] openquake/fdha/test/benchmark/norcia_case3_iaea/
#       config_norcia_case3_iaea_fdha_logic_tree.xml
#   [T] openquake/fdha/test/fixtures/examples_archive/
#       logic_tree_validation_taiwan/fdha_logic_tree_*.xml
PREFILL = {
    # 'style = all' on the Youngs2003 models: the adapter rejects styles
    # other than all/normal unless set explicitly (calc/model_adapter.py:
    # 126-137 and analogues), and the User Manual recommends the
    # Wells & Coppersmith "all styles" coefficients ("recommended,
    # consistent with paper and fdhpy", primary/Youngs2003.md).
    "Youngs2003PrimarySR": "style = all",
    "Youngs2003PrimaryFD": "norm_disp_type = AD\nstyle = all",  # [N]
    "Takao2013PrimaryFD": "norm_disp_type = AD",        # [T]
    "Moss2024PrimaryFD": "version = AD",                # [T]
    "Youngs2003SecondarySR": "version = 3\nstyle = all",   # [N]
    "Youngs2003SecondaryFD": "percentile = 85\nstyle = all",  # [N]
    "Visini2025SecondarySR": "pixel_size = 100",        # [N] (pixel_size is Required)
    "Visini2025SecondaryFD": "scaling_model = WC1994",  # [N]
}

PARAM_HEADER = re.compile(
    r"^\|\s*Name\s*\|\s*Type\s*\|\s*Units\s*\|\s*Default\s*\|\s*Allowed\s*\|"
)

# Some classes share a doc page with their partner model and the page's
# parameter table does not apply to them. Override with code-derived rows.
CLASS_PARAM_OVERRIDES = {
    # secondary/Youngs2003.md documents only the FD 'percentile' parameter;
    # the SR class takes 'version' (get_prob signature and docstring,
    # openquake/fdha/secondary_surf_rup/youngs2003.py:36-41).
    "Youngs2003SecondarySR": [{
        "name": "version", "type": "string", "units": "–",
        "default": '"3"', "allowed": '"1", "2", "3"', "required": False,
        "description": ('Model equation: "1" original formulation, "2" '
                        'average-site formulation, "3" 50/50 weighted '
                        'average of both (default).'),
    }],
}


def parse_doc_params(md_rel: str) -> list[dict]:
    """Extract the parameter table rows from a User-Manual model page."""
    path = DOCS / md_rel
    if not path.exists():
        return []
    rows, in_table = [], False
    for line in path.read_text().splitlines():
        if PARAM_HEADER.match(line):
            in_table = True
            continue
        if in_table:
            if not line.startswith("|"):
                break
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if not cells or set(cells[0]) <= {"-", " ", ":"}:
                continue
            # Name | Type | Units | Default | Allowed | Required? | Description
            rows.append({
                "name": cells[0].strip("`"),
                "type": cells[1] if len(cells) > 1 else "",
                "units": cells[2] if len(cells) > 2 else "",
                "default": cells[3] if len(cells) > 3 else "",
                "allowed": cells[4] if len(cells) > 4 else "",
                "required": (cells[5].lower().startswith("yes")
                             if len(cells) > 5 else False),
                "description": cells[6] if len(cells) > 6 else "",
            })
    return rows


def snapshot() -> dict:
    out: dict = {"slots": {}, "excluded_by_curation": sorted(EXCLUDED)}
    for utype, pkg_name in SLOTS.items():
        pkg = importlib.import_module(f"openquake.fdha.{pkg_name}")
        models = []
        for name in sorted(dir(pkg)):
            obj = getattr(pkg, name)
            if not (inspect.isclass(obj) and name[0].isupper()):
                continue
            if name.startswith("Base") or name in EXCLUDED:
                continue
            entry: dict = {
                "class_name": name,
                "module": obj.__module__,
                "doc_page": DOC_MAP.get(name),
                "doc_params": (CLASS_PARAM_OVERRIDES.get(name)
                               or (parse_doc_params(DOC_MAP[name])
                                   if name in DOC_MAP else [])),
                "prefill": PREFILL.get(name, ""),
            }
            # Constructor params with explicit defaults (e.g. n_sigma)
            try:
                sig_params = inspect.signature(obj.__init__).parameters
                ctor = {
                    p.name: p.default
                    for p in list(sig_params.values())[1:]
                    if p.kind not in (p.VAR_POSITIONAL, p.VAR_KEYWORD)
                    and p.default is not inspect.Parameter.empty
                    and isinstance(p.default, (int, float, str, bool))
                }
                if ctor:
                    entry["ctor_defaults"] = ctor
            except (TypeError, ValueError):
                pass
            models.append(entry)
        out["slots"][utype] = {"package": f"openquake.fdha.{pkg_name}",
                               "models": models}
    return out


if __name__ == "__main__":
    data = snapshot()
    dest = Path(__file__).parent / "registry.json"
    dest.write_text(json.dumps(data, indent=2) + "\n")
    n = sum(len(v["models"]) for v in data["slots"].values())
    nd = sum(1 for v in data["slots"].values() for m in v["models"] if m["doc_params"])
    print(f"wrote {dest}: {n} model classes, {nd} with documented parameter tables")
