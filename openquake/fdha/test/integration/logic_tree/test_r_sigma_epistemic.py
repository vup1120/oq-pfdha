"""C3 tests: fdhaCalcRSigma enumeration + execution (MODE B).

Battery for r_sigma_km (mapping accuracy) as epistemic logic-tree branches:

a. EQUIVALENCE — MODE B with a single sigma branch (weight 1.0, value 0) is
   bit-identical to MODE A with no r_sigma_km at all (sigma = 0 is the
   default: the boxcar W_p path), compared against the frozen fixtures
   (``test/fixtures/r_threshold_baseline/``). Covered for the curve job and
   the map job.
b. ANALYTIC AGGREGATION (mean linearity, verification V7 of
   docs/design/rupture_location_uncertainty.md) — for two branches
   (S1,w1),(S2,w2) the aggregated mean rate at EVERY displacement level
   equals ``w1*rate_single(S1) + w2*rate_single(S2)`` where
   ``rate_single(Si)`` are independent MODE A runs. Tolerance <= 1e-12
   relative.
c. QUANTILES — fractile outputs equal weighted empirical quantiles computed
   independently in this file from the per-branch rates and manifest weights.
d. MANIFEST — all end branches listed; combined weights sum to 1.0;
   sigma branch IDs and values are visible under fdha_calc_params.
Plus the end-to-end conflict rule (scalar + branch set -> ConfigurationError).

Sigma values used here are TEST DATA: the example site sits ~0.39 km off the
trace, outside the 0.1 km boxcar (sigma = 0 branch: no principal term) but
inside the +-2 sigma support of the sigma = 0.3 km branch
(W_p ~ exp(-0.39^2/(2*0.3^2)) ~ 0.43), so the two branches produce genuinely
different rate curves (asserted in-test). They are not recommendations of any
kind — Petersen's own two-sided classes are 0.027-0.116 km.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import pytest

from openquake.fdha.calc.config_loader import ConfigurationError
from openquake.fdha.logic_tree.driver import FdhaLogicTree

REPO_ROOT = Path(__file__).resolve().parents[5]
EXAMPLES_DIR = REPO_ROOT / "examples"
FIXTURES_DIR = (
    Path(__file__).resolve().parents[2] / "fixtures" / "r_threshold_baseline"
)

S1, W1 = 0.0, 0.4
S2, W2 = 0.3, 0.6


def _sigma_branching_level(branches: list[tuple[str, str, str]]) -> str:
    lines = [
        '    <logicTreeBranchingLevel branchingLevelID="bl_5_r_sigma">',
        '      <logicTreeBranchSet branchSetID="bs_5_r_sigma" '
        'uncertaintyType="fdhaCalcRSigma">',
    ]
    for bid, value, weight in branches:
        lines += [
            f'        <logicTreeBranch branchID="{bid}">',
            f'          <uncertaintyModel>{value}</uncertaintyModel>',
            f'          <uncertaintyWeight>{weight}</uncertaintyWeight>',
            '        </logicTreeBranch>',
        ]
    lines += ['      </logicTreeBranchSet>', '    </logicTreeBranchingLevel>']
    return "\n".join(lines)


def _make_curve_job(
    workdir: Path,
    *,
    scalar: str | None,
    branches: list[tuple[str, str, str]] | None,
) -> Path:
    """Materialise a copy of examples/hazard_curve_minimal.ini in workdir.

    ``scalar`` is an ``r_sigma_km`` value to insert into [calculation]
    (MODE A); the fixed ``r_threshold_km = 0.1`` scalar of the example INI is
    never touched (it is the sigma=0 boxcar half-width, not the epistemic
    quantity, and may legitimately coexist with a sigma branch set).
    """
    workdir.mkdir(parents=True, exist_ok=True)
    shutil.copy(EXAMPLES_DIR / "source_model.xml", workdir)
    shutil.copy(
        EXAMPLES_DIR / "hazard_curve_minimal_source_model_logic_tree.xml", workdir
    )
    lt_text = (EXAMPLES_DIR / "hazard_curve_minimal_fdha_logic_tree.xml").read_text()
    if branches:
        lt_text = lt_text.replace(
            "  </logicTree>", _sigma_branching_level(branches) + "\n  </logicTree>"
        )
    (workdir / "hazard_curve_minimal_fdha_logic_tree.xml").write_text(lt_text)

    ini_text = (EXAMPLES_DIR / "hazard_curve_minimal.ini").read_text()
    if scalar is not None:
        ini_text = ini_text.replace(
            "r_threshold_km = 0.1",
            f"r_threshold_km = 0.1\nr_sigma_km = {scalar}",
        )
    ini = workdir / "job.ini"
    ini.write_text(ini_text)
    return ini


def _make_map_job(
    workdir: Path, *, branches: list[tuple[str, str, str]] | None
) -> Path:
    """Materialise a copy of examples/hazard_map_minimal.ini in workdir."""
    workdir.mkdir(parents=True, exist_ok=True)
    shutil.copy(EXAMPLES_DIR / "source_model_char.xml", workdir)
    shutil.copy(
        EXAMPLES_DIR / "hazard_map_minimal_source_model_logic_tree.xml", workdir
    )
    lt_text = (EXAMPLES_DIR / "hazard_map_minimal_fdha_logic_tree.xml").read_text()
    if branches:
        lt_text = lt_text.replace(
            "  </logicTree>", _sigma_branching_level(branches) + "\n  </logicTree>"
        )
    (workdir / "hazard_map_minimal_fdha_logic_tree.xml").write_text(lt_text)
    shutil.copy(EXAMPLES_DIR / "hazard_map_minimal.ini", workdir / "job.ini")
    return workdir / "job.ini"


def _branch_rates(outdir: Path, idx: int) -> np.ndarray:
    return np.loadtxt(
        outdir / "hazard_curves" / f"branch_{idx:04d}.csv",
        delimiter=",", skiprows=1, usecols=1,
    )


# ------------------------------------------------------------ a. EQUIVALENCE

@pytest.mark.integration
@pytest.mark.regression
def test_mode_b_single_zero_branch_equals_mode_a_curve_fixture(tmp_path):
    """MODE B (one branch, weight 1.0, sigma 0) == the frozen curve fixture
    (which never defined r_sigma_km, i.e. the sigma = 0 boxcar default)."""
    ini = _make_curve_job(
        tmp_path / "job", scalar=None, branches=[("RS_V", "0.0", "1.0")]
    )
    lt = FdhaLogicTree.from_ini(str(ini))
    lt.run(outdir=tmp_path / "out")

    fix = FIXTURES_DIR / "curve_explicit"
    assert (tmp_path / "out" / "aggregate_hazard.csv").read_bytes() == (
        fix / "aggregate_hazard.csv"
    ).read_bytes()
    assert (
        tmp_path / "out" / "hazard_curves" / "branch_0000.csv"
    ).read_bytes() == (fix / "branch_0000.csv").read_bytes()

    # The branch INI must carry the materialised value on the same key the
    # MODE A scalar uses (single shared consumption point).
    branch_ini = (
        tmp_path / "out" / "source_model_branches"
        / "00_b_sm_hazard_curve_minimal" / "branch_configs" / "branch_0000.ini"
    )
    assert "r_sigma_km = 0.0" in branch_ini.read_text()


@pytest.mark.integration
@pytest.mark.regression
def test_mode_b_single_zero_branch_equals_mode_a_map_fixture(tmp_path):
    """MODE B sigma-0 branch == the map fixture (no r_sigma_km, boxcar path):
    bit-identical CSVs / <=1e-12 HDF5 rates."""
    ini = _make_map_job(tmp_path / "job", branches=[("RS_V", "0.0", "1.0")])
    lt = FdhaLogicTree.from_ini(str(ini))
    lt.run(outdir=tmp_path / "out")

    fix = FIXTURES_DIR / "map_default"
    for name in (
        "displacement_map_mean.csv",
        "displacement_map_quantile-0.05.csv",
        "displacement_map_quantile-0.16.csv",
        "displacement_map_quantile-0.5.csv",
        "displacement_map_quantile-0.84.csv",
        "displacement_map_quantile-0.95.csv",
    ):
        assert (tmp_path / "out" / "aggregate" / name).read_bytes() == (
            fix / name
        ).read_bytes(), f"aggregate/{name} deviates from MODE A baseline"

    h5py = pytest.importorskip("h5py")
    frozen = np.load(fix / "rates_baseline.npz")
    branch_h5 = (
        tmp_path / "out" / "source_model_branches"
        / "00_b_sm_hazard_map_minimal" / "branches" / "branch_0000.h5"
    )
    with h5py.File(branch_h5) as f:
        np.testing.assert_allclose(
            f["rates"][:], frozen["branch_rates"], rtol=1e-12, atol=0
        )


# ------------------- b. ANALYTIC AGGREGATION (mean linearity, V7)

@pytest.fixture(scope="module")
def two_branch_run(tmp_path_factory):
    """One MODE B two-branch run + two independent MODE A runs (curve)."""
    root = tmp_path_factory.mktemp("rs_two_branch")

    ini_b = _make_curve_job(
        root / "mode_b", scalar=None,
        branches=[("RS_A", str(S1), str(W1)), ("RS_B", str(S2), str(W2))],
    )
    res_b = FdhaLogicTree.from_ini(str(ini_b)).run(outdir=root / "out_b")

    singles = {}
    for tag, s in (("s1", S1), ("s2", S2)):
        ini_a = _make_curve_job(
            root / f"mode_a_{tag}", scalar=str(s), branches=None,
        )
        FdhaLogicTree.from_ini(str(ini_a)).run(outdir=root / f"out_a_{tag}")
        singles[tag] = _branch_rates(root / f"out_a_{tag}", 0)

    return {
        "result": res_b,
        "outdir": root / "out_b",
        "rate_s1": singles["s1"],
        "rate_s2": singles["s2"],
    }


@pytest.mark.integration
def test_analytic_aggregation_mean(two_branch_run):
    rate_s1 = two_branch_run["rate_s1"]
    rate_s2 = two_branch_run["rate_s2"]
    # Guard against a vacuous test: the two sigmas must matter (the site is
    # outside the boxcar but inside the Gaussian's +-2 sigma support).
    assert not np.allclose(rate_s1, rate_s2), (
        "sigma branches produce identical rates; probe site/values "
        "no longer discriminate"
    )

    mean_lt = np.asarray(two_branch_run["result"].mean_rates, dtype=float)[0]
    expected = W1 * rate_s1 + W2 * rate_s2
    np.testing.assert_allclose(mean_lt, expected, rtol=1e-12, atol=0)


@pytest.mark.integration
def test_mode_b_branches_match_independent_mode_a_runs(two_branch_run):
    """Each MODE B per-branch curve equals its independent MODE A run."""
    outdir = two_branch_run["outdir"]
    manifest = json.loads((outdir / "manifest.json").read_text())
    by_value = {
        b["fdha_calc_params"]["r_sigma_km"]: b["global_index"]
        for b in manifest["branches"]
    }
    np.testing.assert_allclose(
        _branch_rates(outdir, by_value[S1]), two_branch_run["rate_s1"],
        rtol=0, atol=0,
    )
    np.testing.assert_allclose(
        _branch_rates(outdir, by_value[S2]), two_branch_run["rate_s2"],
        rtol=0, atol=0,
    )


# ----------------------------------------------------------- c. QUANTILES

def _weighted_empirical_quantile(values, weights, q):
    """Independent reference: weighted empirical quantile, linear in CDF."""
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    order = np.argsort(values)
    xs, ws = values[order], weights[order] / weights.sum()
    cdf = np.cumsum(ws)
    if q <= cdf[0]:
        return xs[0]
    if q >= cdf[-1]:
        return xs[-1]
    j = int(np.searchsorted(cdf, q, side="left"))
    t = (q - cdf[j - 1]) / (cdf[j] - cdf[j - 1])
    return xs[j - 1] + t * (xs[j] - xs[j - 1])


@pytest.mark.integration
def test_quantiles_equal_weighted_empirical_quantiles(two_branch_run):
    outdir = two_branch_run["outdir"]
    manifest = json.loads((outdir / "manifest.json").read_text())
    branch_rates = np.stack(
        [_branch_rates(outdir, b["global_index"]) for b in manifest["branches"]]
    )
    weights = np.array(
        [b["combined_branch_weight"] for b in manifest["branches"]], dtype=float
    )

    fractiles = two_branch_run["result"].fractiles
    for q, cube in fractiles.items():
        got = np.asarray(cube, dtype=float)[0]
        expected = np.array(
            [
                _weighted_empirical_quantile(branch_rates[:, j], weights, q)
                for j in range(branch_rates.shape[1])
            ]
        )
        np.testing.assert_allclose(got, expected, rtol=1e-12, atol=0)


# ------------------------------------------------------------ d. MANIFEST

@pytest.mark.integration
def test_manifest_lists_all_end_branches_and_weights(two_branch_run):
    outdir = two_branch_run["outdir"]
    manifest = json.loads((outdir / "manifest.json").read_text())

    assert manifest["n_realizations"] == 2
    weights = [b["combined_branch_weight"] for b in manifest["branches"]]
    assert abs(sum(weights) - 1.0) < 1e-12
    assert manifest["total_combined_weight"] == pytest.approx(1.0, abs=1e-12)

    # Sigma branch IDs appear in the composed branch path; the chosen
    # value is recorded under fdha_calc_params.
    ids = {b["fdha_branch_id"].split("|")[0] for b in manifest["branches"]}
    assert ids == {"RS_A", "RS_B"}
    values = {b["fdha_calc_params"]["r_sigma_km"] for b in manifest["branches"]}
    assert values == {S1, S2}
    # Model slots keep listing only real model classes.
    for b in manifest["branches"]:
        assert "calc_r_sigma" not in b["fdha_models"]
        assert set(b["fdha_models"]) == {
            "primary_surf_rup", "primary_surf_displ",
            "secondary_surf_rup", "secondary_surf_displ",
        }


# --------------------------------------------------- conflict rule (driver)

@pytest.mark.integration
def test_conflict_scalar_plus_branch_set_fails_at_startup(tmp_path):
    ini = _make_curve_job(
        tmp_path / "job", scalar="0.1", branches=[("RS_V", "0.2", "1.0")]
    )
    lt = FdhaLogicTree.from_ini(str(ini))
    with pytest.raises(ConfigurationError, match="mutually exclusive"):
        lt.run(outdir=tmp_path / "out")
    # No branch may have been executed.
    assert not (tmp_path / "out" / "hazard_curves").exists()
