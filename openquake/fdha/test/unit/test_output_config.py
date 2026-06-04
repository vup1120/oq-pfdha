"""Unit tests for the OpenQuake-style ``[output]`` config (quantiles + mean).

These cover the pure configuration layer (no numpy/engine required) plus the
OpenQuake-style fractile label helper.
"""
import pytest

from openquake.fdha.calc.config_loader import (
    ConfigValidationError,
    parse_quantiles,
    resolve_output_mean,
    resolve_output_quantiles,
    _normalize_output_section,
)


def test_parse_quantiles_whitespace_string():
    assert parse_quantiles("0.05 0.16 0.5 0.84 0.95") == [0.05, 0.16, 0.5, 0.84, 0.95]


def test_parse_quantiles_sorted_and_deduped():
    assert parse_quantiles("0.9 0.1 0.5 0.5") == [0.1, 0.5, 0.9]


def test_parse_quantiles_list_and_scalar():
    assert parse_quantiles([0.16, 0.84]) == [0.16, 0.84]
    assert parse_quantiles(0.5) == [0.5]


def test_parse_quantiles_empty_is_none_requested():
    assert parse_quantiles("") == []
    assert parse_quantiles(None) == []


@pytest.mark.parametrize("bad", ["1.0", "0", "-0.1", "1.5", "abc"])
def test_parse_quantiles_rejects_out_of_range(bad):
    with pytest.raises(ConfigValidationError):
        parse_quantiles(bad)


def test_normalize_absent_section_uses_defaults():
    cfg = {}
    _normalize_output_section(cfg)
    assert cfg["output"] == {"quantiles": None, "mean": True}
    assert resolve_output_quantiles(cfg) == [0.05, 0.16, 0.5, 0.84, 0.95]
    assert resolve_output_mean(cfg) is True


def test_normalize_explicit_quantiles_and_mean():
    cfg = {"output": {"quantiles": "0.1 0.9", "mean": False}}
    _normalize_output_section(cfg)
    assert resolve_output_quantiles(cfg) == [0.1, 0.9]
    assert resolve_output_mean(cfg) is False


def test_normalize_empty_quantiles_means_no_fractiles():
    cfg = {"output": {"quantiles": ""}}
    _normalize_output_section(cfg)
    assert resolve_output_quantiles(cfg) == []
    # mean still defaults to True
    assert resolve_output_mean(cfg) is True


def test_normalize_mean_string_coercion():
    cfg = {"output": {"mean": "no"}}
    _normalize_output_section(cfg)
    assert resolve_output_mean(cfg) is False


def test_quantile_label_openquake_style():
    # io imports numpy; skip cleanly where the scientific stack is absent.
    io = pytest.importorskip("openquake.fdha.logic_tree.io")
    assert io.quantile_label(0.05) == "quantile-0.05"
    assert io.quantile_label(0.5) == "quantile-0.5"
    assert io.quantile_label(0.025) == "quantile-0.025"
    # Default labels are derived from FRACTILE_QS via quantile_label.
    assert io.FRACTILE_LABELS == tuple(io.quantile_label(q) for q in io.FRACTILE_QS)
