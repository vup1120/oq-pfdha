from __future__ import annotations

import ast
import configparser
import math
from typing import Any

R_THRESHOLD_KM_KEY = "r_threshold_km"


def parse_uncertainty_model(text: str) -> tuple[str, dict[str, Any]]:
    """
    Parse <uncertaintyModel> content.

    Supports:
    - Plain class name: "Chiou2025PrimaryFD"
    - INI-block: "[ClassName]\\nkey = value\\n..."
    """
    raw = (text or "").strip()
    if not raw:
        raise ValueError("Empty uncertaintyModel")

    # Plain class name
    if raw.startswith("[") and "]" in raw.splitlines()[0]:
        return _parse_ini_block(raw)

    if "\n" in raw:
        # tolerate wrapped text but treat as plain token after stripping whitespace
        raw = "".join(line.strip() for line in raw.splitlines() if line.strip())

    return raw, {}


def _parse_ini_block(raw: str) -> tuple[str, dict[str, Any]]:
    # NRML files indent the block to the XML nesting depth (oq-engine GMPE
    # logic-tree style) and values are single-line, so strip every line:
    # otherwise configparser would treat a deeper-indented "key = value" as
    # a continuation of the previous value.
    lines = [ln.strip() for ln in raw.splitlines() if ln.strip()]
    if not lines:
        raise ValueError("Empty INI block")

    header = lines[0].strip()
    if not (header.startswith("[") and header.endswith("]")):
        raise ValueError("INI block must start with [ClassName]")
    class_name = header[1:-1].strip()
    if not class_name:
        raise ValueError("Empty class name in INI header")

    # ConfigParser requires at least one section header; we already have it.
    cp = configparser.RawConfigParser()
    cp.optionxform = str
    cp.read_string("\n".join(lines))

    params: dict[str, Any] = {}
    if cp.has_section(class_name):
        for k, v in cp.items(class_name):
            params[k] = _parse_value(v)
    return class_name, params


def parse_r_threshold_model(text: str) -> float:
    """
    Parse the <uncertaintyModel> of a ``fdhaCalcRThreshold`` branch.

    Follows the OpenQuake engine convention for scalar uncertainty types
    (:func:`openquake.hazardlib.lt.parse_uncertainty` fallback): the element
    text is a single bare float, and the parameter is identified by the
    ``uncertaintyType`` itself. Here the value is the principal/distributed
    distance threshold in kilometres and must be a positive finite float.

    The branch value is one epistemic alternative for the hard-step
    simplification of the rupture-location term fr(r) of Petersen et al.
    (2011, BSSA 101, 805-825, doi:10.1785/0120100035); representing the
    choice as weighted branches follows Petersen et al. (2011, p. 810) and
    IAEA-TECDOC-2092 (2025, Section 3.3).

    Raises ValueError for anything else (empty text, several tokens, a
    ``key = value`` line, a non-numeric / non-finite / non-positive value).
    """
    raw = (text or "").strip()
    try:
        value = float(raw)
    except (TypeError, ValueError):
        raise ValueError(
            "fdhaCalcRThreshold: expected single positive float value (km) "
            f"in <uncertaintyModel>, got {raw!r}"
        )
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(
            f"fdhaCalcRThreshold value must be a positive finite float (km); got {raw!r}"
        )
    return value


def _parse_value(v: str) -> Any:
    s = v.strip()
    if s == "":
        return ""
    try:
        # literal_eval handles numbers, lists, dicts, strings with quotes, booleans, None
        return ast.literal_eval(s)
    except Exception:
        return s

