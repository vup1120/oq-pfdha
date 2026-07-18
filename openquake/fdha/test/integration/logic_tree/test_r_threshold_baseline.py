"""MODE A regression baseline for the scalar ``[calculation].r_threshold_km``.

These tests pin the outputs of the two public example jobs against golden
fixtures (``test/fixtures/r_threshold_baseline/``). They guarantee that jobs
using the INI scalar (``curve_explicit``) or relying on the implementation
default when the key is absent (``map_default``) keep producing bit-identical
results.

Tolerance contract (documented per the baseline task spec):

- All CSV outputs and the materialised branch INI are deterministic text and
  are compared **byte-exact** (run-to-run byte-identity of both jobs was
  verified when the fixtures were frozen).
- ``manifest.json`` is deterministic except for one embedded absolute
  filesystem path (``source_model_file``); it is compared structurally after
  normalising absolute paths to basenames. Weights are part of the structure
  and therefore compared exactly.
- Rates stored in HDF5 (map mode) are compared at ``<= 1e-12`` **relative**
  tolerance (``rtol=1e-12, atol=0``): HDF5 round-trips float64 exactly, but
  byte-level file identity is not guaranteed across h5py/libhdf5 versions.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from openquake.fdha.logic_tree.driver import FdhaLogicTree

REPO_ROOT = Path(__file__).resolve().parents[5]
EXAMPLES_DIR = REPO_ROOT / "examples"
FIXTURES_DIR = (
    Path(__file__).resolve().parents[2] / "fixtures" / "r_threshold_baseline"
)


def _normalize_paths(obj):
    """Replace absolute filesystem paths in JSON-like data with basenames."""
    if isinstance(obj, dict):
        return {k: _normalize_paths(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_normalize_paths(v) for v in obj]
    if isinstance(obj, str) and obj.startswith("/"):
        return Path(obj).name
    return obj


def _assert_manifest_matches(fresh: Path, frozen: Path) -> None:
    got = _normalize_paths(json.loads(fresh.read_text()))
    want = _normalize_paths(json.loads(frozen.read_text()))
    assert got == want


@pytest.mark.integration
@pytest.mark.regression
def test_curve_job_with_explicit_r_threshold_matches_baseline(tmp_path):
    """hazard_curve_minimal.ini sets r_threshold_km = 0.1 explicitly."""
    lt = FdhaLogicTree.from_ini(str(EXAMPLES_DIR / "hazard_curve_minimal.ini"))
    res = lt.run(outdir=tmp_path / "out")
    assert res.mode == "hazard_curve"

    fix = FIXTURES_DIR / "curve_explicit"
    out = tmp_path / "out"

    # Deterministic text outputs: byte-exact.
    assert (out / "aggregate_hazard.csv").read_bytes() == (
        fix / "aggregate_hazard.csv"
    ).read_bytes()
    assert (out / "hazard_curves" / "branch_0000.csv").read_bytes() == (
        fix / "branch_0000.csv"
    ).read_bytes()

    _assert_manifest_matches(out / "manifest.json", fix / "manifest.json")

    # The materialised branch INI must carry the scalar through unchanged.
    branch_ini = (
        out
        / "source_model_branches"
        / "00_b_sm_hazard_curve_minimal"
        / "branch_configs"
        / "branch_0000.ini"
    )
    assert "r_threshold_km = 0.1" in branch_ini.read_text()


@pytest.mark.integration
@pytest.mark.regression
def test_map_job_with_default_r_threshold_matches_baseline(tmp_path):
    """hazard_map_minimal.ini omits r_threshold_km (implementation default)."""
    lt = FdhaLogicTree.from_ini(str(EXAMPLES_DIR / "hazard_map_minimal.ini"))
    res = lt.run(outdir=tmp_path / "out")
    assert res.mode == "hazard_map"

    fix = FIXTURES_DIR / "map_default"
    out = tmp_path / "out"

    # Deterministic text outputs: byte-exact.
    for name in (
        "displacement_map_mean.csv",
        "displacement_map_quantile-0.05.csv",
        "displacement_map_quantile-0.16.csv",
        "displacement_map_quantile-0.5.csv",
        "displacement_map_quantile-0.84.csv",
        "displacement_map_quantile-0.95.csv",
    ):
        assert (out / "aggregate" / name).read_bytes() == (
            fix / name
        ).read_bytes(), f"aggregate/{name} deviates from frozen baseline"

    _assert_manifest_matches(out / "manifest.json", fix / "manifest.json")

    # The job omits the key: the materialised branch INI must NOT introduce it.
    branch_dir = out / "source_model_branches" / "00_b_sm_hazard_map_minimal"
    branch_ini = branch_dir / "branch_configs" / "branch_0000.ini"
    assert "r_threshold_km" not in branch_ini.read_text()

    # HDF5-borne rates: <= 1e-12 relative.
    h5py = pytest.importorskip("h5py")
    frozen = np.load(fix / "rates_baseline.npz")

    with h5py.File(branch_dir / "branches" / "branch_0000.h5") as f:
        np.testing.assert_allclose(
            f["rates"][:], frozen["branch_rates"], rtol=1e-12, atol=0
        )
        np.testing.assert_array_equal(f["d0"][:], frozen["d0"])
        np.testing.assert_array_equal(f["site_lons"][:], frozen["site_lons"])
        np.testing.assert_array_equal(f["site_lats"][:], frozen["site_lats"])
        assert float(f.attrs["weight"]) == float(frozen["branch_weight"][0])
        assert str(f.attrs["fingerprint"]) == str(frozen["branch_fingerprint"][0])

    with h5py.File(out / "aggregate" / "rates_mean.h5") as f:
        np.testing.assert_allclose(
            f["rates_mean"][:], frozen["rates_mean"], rtol=1e-12, atol=0
        )
    with h5py.File(out / "aggregate" / "rates_fractiles.h5") as f:
        np.testing.assert_allclose(
            f["rates_fractiles"][:], frozen["rates_fractiles"], rtol=1e-12, atol=0
        )
        np.testing.assert_array_equal(f["quantiles"][:], frozen["quantiles"])
