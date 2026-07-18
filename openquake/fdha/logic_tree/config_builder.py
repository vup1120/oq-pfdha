"""
Build the per-branch calculation configuration from an end branch and
the base job configuration.
"""
from __future__ import annotations

import copy
import json
from configparser import RawConfigParser
from io import StringIO
from typing import Any

from openquake.fdha.logic_tree.types import CALC_SLOTS_BY_UTYPE, EndBranch

# Pseudo-slots holding calculation parameters rather than model choices.
CALC_SLOTS = frozenset(CALC_SLOTS_BY_UTYPE.values())


def build_config(base_config: dict[str, Any], end_branch: EndBranch) -> dict[str, Any]:
    cfg = copy.deepcopy(base_config)
    cfg.setdefault("models", {})
    for slot, choice in end_branch.selections.items():
        if slot in CALC_SLOTS:
            # Calc-param branch: write into [calculation] so the branch INI
            # feeds the exact same runtime read as a MODE A scalar job
            # (calculators.BaseFaultRuptureCalculator._initialize_calculation_params).
            # The conflict rule guarantees the job INI did not set the key.
            cfg.setdefault("calculation", {}).update(choice.params)
            continue
        cfg["models"][slot] = {"type": choice.class_name, "parameters": choice.params}
    # Logic tree NRML pointers are driver-level concepts, not branch calculators.
    if "calculation" in cfg and isinstance(cfg["calculation"], dict):
        cfg["calculation"].pop("fdha_logic_tree_files", None)
        cfg["calculation"].pop("fdha_logic_tree_file", None)
    return cfg


def dump_config_to_ini(config: dict[str, Any]) -> str:
    """
    Dump nested config dict to INI text compatible with config_loader.load_config.
    """
    cp = RawConfigParser()
    cp.optionxform = str

    for section, body in config.items():
        if not isinstance(body, dict):
            # treat scalars at top-level as [section] key 'value'
            cp[section] = {"value": _dump_value(body)}
            continue
        if section != "models":
            cp[section] = {k: _dump_value(v) for k, v in body.items()}
            continue

        # models.<slot>
        for slot, m in body.items():
            if not isinstance(m, dict):
                continue
            sec = f"models.{slot}"
            cp[sec] = {}
            if "type" in m:
                cp[sec]["type"] = str(m["type"])
            params = m.get("parameters", {})
            if isinstance(params, dict):
                psec = f"models.{slot}.parameters"
                cp[psec] = {pk: _dump_value(pv) for pk, pv in params.items()}

    buf = StringIO()
    cp.write(buf)
    return buf.getvalue()


def _dump_value(v: Any) -> str:
    if isinstance(v, (dict, list)):
        return json.dumps(v)
    if v is None:
        return "null"
    if isinstance(v, bool):
        return "true" if v else "false"
    return str(v)

