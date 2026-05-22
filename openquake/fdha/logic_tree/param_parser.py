from __future__ import annotations

import ast
import configparser
from typing import Any


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
    lines = [ln.rstrip() for ln in raw.splitlines()]
    # drop leading/trailing empty
    while lines and not lines[0].strip():
        lines.pop(0)
    while lines and not lines[-1].strip():
        lines.pop()
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


def _parse_value(v: str) -> Any:
    s = v.strip()
    if s == "":
        return ""
    try:
        # literal_eval handles numbers, lists, dicts, strings with quotes, booleans, None
        return ast.literal_eval(s)
    except Exception:
        return s

