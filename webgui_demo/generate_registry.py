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
    # Petersen2011PrimaryFD shape variants are convenience subclasses that only
    # hardcode `version` (primary_surf_displ/petersen2011.py:190-222). The GUI
    # exposes the base class with `version` as a selectable parameter instead,
    # so these are excluded to avoid redundant duplicate entries.
    "Petersen2011PrimaryFD_bilinear",
    "Petersen2011PrimaryFD_elliptical",
    "Petersen2011PrimaryFD_quadratic",
}

# Class -> User-Manual page holding its parameter table.
DOC_MAP = {
    "Youngs2003PrimaryFD": "primary/Youngs2003.md",
    "Petersen2011PrimaryFD": "primary/Petersen2011.md",
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
    # Youngs et al. (2003) is a normal-faulting model; its `style` parameter
    # selects the coefficient dataset ("all" = WC1994 all-styles, "normal" =
    # normal-only), NOT the source mechanism. It is left unset here so the
    # engine uses normal coefficients on normal sources; the GUI warns when
    # the source is reverse/strike-slip (see STYLE_CONSTRAINTS in app.py).
    "Youngs2003PrimaryFD": "norm_disp_type = AD",       # [N]
    "Takao2013PrimaryFD": "norm_disp_type = AD",        # [T]
    "MossRoss2011PrimaryFD": "norm_disp_type = AD",     # required (no default)
    "Moss2024PrimaryFD": "version = AD",                # [T]
    "Youngs2003SecondarySR": "version = 3",             # [N]
    "Youngs2003SecondaryFD": "percentile = 85",         # [N]
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
    # MossRoss2011PrimaryFD shares Takao2013's signature
    # get_prob(d, X_L_ratio, mag, norm_disp_type) but has no doc page, so the
    # required norm_disp_type knob was undocumented (moss_ross2011.py:58).
    "MossRoss2011PrimaryFD": [{
        "name": "norm_disp_type", "type": "string", "units": "–",
        "default": "–", "allowed": '`"AD"`, `"MD"`', "required": True,
        "description": ('Normalization displacement type. `"AD"` uses average '
                        'displacement normalization (Gamma distribution), '
                        '`"MD"` uses maximum displacement normalization (Beta '
                        'distribution).'),
    }],
}

# Code-derived optional parameters that are absent from the User-Manual tables
# (these classes either have no doc page, or the page documents only a subset).
# Each list is APPENDED to whatever doc_params were parsed/overridden, so it
# complements rather than replaces. All rows are optional (the model supplies a
# working default); allowed values are taken verbatim from the get_prob
# signatures and validation in the model source.
CLASS_PARAM_EXTRA = {
    # primary_surf_displ/moss2022.py:67-69
    "Moss2022PrimaryFD": [
        {"name": "version", "type": "string", "units": "–", "default": '"MD"',
         "allowed": '`"AD"`, `"MD"`', "required": False,
         "description": "Normalization: average (AD) or maximum (MD) displacement."},
        {"name": "completeness", "type": "string", "units": "–",
         "default": '"complete"', "allowed": '`"complete"`, `"incomplete"` (MD only), `"all"`',
         "required": False, "description": "Measurement completeness subset."},
        {"name": "gamma_mode", "type": "string", "units": "–",
         "default": '"regression"', "allowed": '`"regression"`, `"global"`',
         "required": False, "description": "Gamma shape source: x/L regression or global (Eqs 4.2-4.3)."},
        {"name": "sigma_type", "type": "string", "units": "–",
         "default": '"recommended"', "allowed": '`"recommended"`, `"regression"`',
         "required": False, "description": "Sigma model: revised Table 4.4 or regression."},
    ],
    # secondary_surf_displ/moss2022.py:56-59
    "Moss2022SecondaryFD": [
        {"name": "version", "type": "string", "units": "–", "default": '"MD"',
         "allowed": '`"AD"`, `"MD"`', "required": False,
         "description": "Normalization: average (AD) or maximum (MD) displacement."},
        {"name": "completeness", "type": "string", "units": "–",
         "default": '"complete"', "allowed": '`"complete"`, `"incomplete"` (MD only), `"all"`',
         "required": False, "description": "Measurement completeness subset."},
        {"name": "sigma_type", "type": "string", "units": "–",
         "default": '"recommended"', "allowed": '`"recommended"`, `"regression"`',
         "required": False, "description": "Sigma model: revised or regression."},
        {"name": "faulting", "type": "string", "units": "–", "default": '"simple"',
         "allowed": '`"simple"`, `"complex"`', "required": False,
         "description": "Faulting complexity branch."},
        {"name": "percentile", "type": "string", "units": "–", "default": '"85"',
         "allowed": '`"50"`, `"85"`', "required": False,
         "description": "Percentile of the d/MD envelope."},
        {"name": "method", "type": "string", "units": "–", "default": '"gamma"',
         "allowed": '`"gamma"`, `"envelope"`', "required": False,
         "description": "Distribution method for the displacement model."},
        {"name": "gamma_a", "type": "float", "units": "–", "default": "global",
         "allowed": ">0", "required": False,
         "description": "Override the global gamma shape parameter (gamma method only)."},
    ],
    # secondary_surf_rup/petersen2011.py:53 (appended to documented cell_size)
    "Petersen2011SecondarySR": [
        {"name": "version", "type": "string", "units": "–", "default": '"default"',
         "allowed": '`"default"`, `"near_field"`', "required": False,
         "description": "Far-field power function (Table 4) or near-field interpolation (Table 5)."},
    ],
    # secondary_surf_rup/ferrario2021.py:53
    "FerrarioLivio2021SecondarySR": [
        {"name": "version", "type": "string", "units": "–", "default": '"regular"',
         "allowed": '`"regular"`, `"conservative"`', "required": False,
         "description": "Standard or conservative (higher-probability) coefficient set."},
    ],
    # secondary_surf_rup/moss2022.py:67
    "Moss2022SecondarySR": [
        {"name": "method", "type": "string", "units": "–", "default": '"simple"',
         "allowed": '`"simple"`, `"biexp"`', "required": False,
         "description": "Decay model: simple exponential or biexponential (needs mag)."},
    ],
    # secondary_surf_rup/rodriguez2023.py:45
    "Rodriguez2023SecondarySR": [
        {"name": "pixel_size", "type": "integer", "units": "meters", "default": "1",
         "allowed": "1", "required": False,
         "description": "Analysis pixel width; the model is calibrated for 1 m only."},
    ],
    # secondary_surf_rup/takao2013.py:46
    "Takao2013SecondarySR": [
        {"name": "pixel_size", "type": "integer", "units": "meters", "default": "500",
         "allowed": "500", "required": False,
         "description": "Analysis pixel width; the model is calibrated for 500 m only."},
    ],
    # secondary_surf_rup/takao2014.py:43
    "Takao2014SecondarySR": [
        {"name": "pixel_size", "type": "integer", "units": "meters", "default": "100",
         "allowed": "500, 250, 100, 50", "required": False,
         "description": "Analysis pixel width for the P(across) coefficient set."},
    ],
    # secondary_surf_displ/visini2025.py:105 (appended to the calculator table)
    "Visini2025SecondaryFD": [
        {"name": "scaling_model", "type": "string", "units": "–", "default": '"WC1994"',
         "allowed": '`"WC1994"`, `"THINGBAIJAM2017"`, `"LEONARD2010"`', "required": False,
         "description": "Scaling relation used to derive TPFm when not supplied explicitly."},
    ],
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
            doc_params = list(CLASS_PARAM_OVERRIDES.get(name)
                              or (parse_doc_params(DOC_MAP[name])
                                  if name in DOC_MAP else []))
            doc_params += CLASS_PARAM_EXTRA.get(name, [])
            entry: dict = {
                "class_name": name,
                "module": obj.__module__,
                "doc_page": DOC_MAP.get(name),
                "doc_params": doc_params,
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
