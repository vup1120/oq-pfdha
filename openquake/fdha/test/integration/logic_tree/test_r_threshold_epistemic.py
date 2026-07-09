"""Phase 3 tests: fdhaCalcRThreshold enumeration + execution (MODE B).

Battery required by the r_threshold-as-epistemic-uncertainty task:

a. EQUIVALENCE — MODE B with a single threshold branch (weight 1.0, value V)
   is bit-identical to MODE A with the scalar ``r_threshold_km = V``,
   compared against the frozen Phase 1 fixtures
   (``test/fixtures/r_threshold_baseline/``). Covered for the curve job
   (fixture sets the scalar explicitly) and the map job (fixture omits the
   scalar, exercising the implementation default 0.1 == V).
b. ANALYTIC AGGREGATION — for two branches (T1,w1),(T2,w2) the aggregated
   mean rate at EVERY displacement level equals
   ``w1*rate_single(T1) + w2*rate_single(T2)`` where ``rate_single(Ti)`` are
   independent MODE A runs. Tolerance <= 1e-12 relative.
c. QUANTILES — fractile outputs equal weighted empirical quantiles computed
   independently in this file from the per-branch rates and manifest weights.
d. MANIFEST — all end branches listed; combined weights sum to 1.0;
   threshold branch IDs and values are visible.
Plus the end-to-end conflict rule (scalar + branch set -> ConfigurationError).

Threshold values used here are TEST DATA chosen so the example site
(offset to fall between them in fault distance) is routed to the
distributed zone for T1 and the principal zone for T2 — they are not
recommendations of any kind.
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

# Off-trace probe site: fault distance lies between T1 and T2 below, so the
# two branches produce genuinely different rate curves (asserted in-test).
OFF_TRACE_SITE = "16.18 39.64704451"
T1, W1 = 0.5, 0.4
T2, W2 = 2.0, 0.6


def _threshold_branching_level(branches: list[tuple[str, str, str]]) -> str:
    lines = [
        '    <logicTreeBranchingLevel branchingLevelID="bl_5_r_threshold">',
        '      <logicTreeBranchSet branchSetID="bs_5_r_threshold" '
        'uncertaintyType="fdhaCalcRThreshold">',
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
    site: str | None = None,
) -> Path:
    """Materialise a copy of examples/hazard_curve_minimal.ini in workdir."""
    workdir.mkdir(parents=True, exist_ok=True)
    shutil.copy(EXAMPLES_DIR / "source_model.xml", workdir)
    shutil.copy(
        EXAMPLES_DIR / "hazard_curve_minimal_source_model_logic_tree.xml", workdir
    )
    lt_text = (EXAMPLES_DIR / "hazard_curve_minimal_fdha_logic_tree.xml").read_text()
    if branches:
        lt_text = lt_text.replace(
            "  </logicTree>", _threshold_branching_level(branches) + "\n  </logicTree>"
        )
    (workdir / "hazard_curve_minimal_fdha_logic_tree.xml").write_text(lt_text)

    ini_text = (EXAMPLES_DIR / "hazard_curve_minimal.ini").read_text()
    if scalar is None:
        ini_text = ini_text.replace("r_threshold_km = 0.1\n", "")
    else:
        ini_text = ini_text.replace(
            "r_threshold_km = 0.1", f"r_threshold_km = {scalar}"
        )
    if site is not None:
        ini_text = ini_text.replace(
            "sites = 16.16573727 39.64704451", f"sites = {site}"
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
            "  </logicTree>", _threshold_branching_level(branches) + "\n  </logicTree>"
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
def test_mode_b_single_branch_equals_mode_a_curve_fixture(tmp_path):
    """MODE B (one branch, weight 1.0, value 0.1) == Phase 1 curve fixture."""
    ini = _make_curve_job(
        tmp_path / "job", scalar=None, branches=[("RT_V", "0.1", "1.0")]
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
    assert "r_threshold_km = 0.1" in branch_ini.read_text()


@pytest.mark.integration
@pytest.mark.regression
def test_mode_b_single_branch_equals_mode_a_map_fixture(tmp_path):
    """MODE B branch value == the implementation default of the map fixture.

    The Phase 1 map fixture omits the scalar, so the run used the
    implementation default (0.1 km). A single MODE B branch with the same
    value must reproduce it bit-identically (CSVs) / <=1e-12 (HDF5 rates).
    """
    ini = _make_map_job(tmp_path / "job", branches=[("RT_V", "0.1", "1.0")])
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


# ------------------------------------------- b. ANALYTIC AGGREGATION (mean)

@pytest.fixture(scope="module")
def two_branch_run(tmp_path_factory):
    """One MODE B two-branch run + two independent MODE A runs (curve)."""
    root = tmp_path_factory.mktemp("rt_two_branch")

    ini_b = _make_curve_job(
        root / "mode_b", scalar=None,
        branches=[("RT_A", str(T1), str(W1)), ("RT_B", str(T2), str(W2))],
        site=OFF_TRACE_SITE,
    )
    res_b = FdhaLogicTree.from_ini(str(ini_b)).run(outdir=root / "out_b")

    singles = {}
    for tag, t in (("t1", T1), ("t2", T2)):
        ini_a = _make_curve_job(
            root / f"mode_a_{tag}", scalar=str(t), branches=None,
            site=OFF_TRACE_SITE,
        )
        FdhaLogicTree.from_ini(str(ini_a)).run(outdir=root / f"out_a_{tag}")
        singles[tag] = _branch_rates(root / f"out_a_{tag}", 0)

    return {
        "result": res_b,
        "outdir": root / "out_b",
        "rate_t1": singles["t1"],
        "rate_t2": singles["t2"],
    }


@pytest.mark.integration
def test_analytic_aggregation_mean(two_branch_run):
    rate_t1 = two_branch_run["rate_t1"]
    rate_t2 = two_branch_run["rate_t2"]
    # Guard against a vacuous test: the two thresholds must matter.
    assert not np.allclose(rate_t1, rate_t2), (
        "threshold branches produce identical rates; probe site/values "
        "no longer discriminate"
    )

    mean_lt = np.asarray(two_branch_run["result"].mean_rates, dtype=float)[0]
    expected = W1 * rate_t1 + W2 * rate_t2
    np.testing.assert_allclose(mean_lt, expected, rtol=1e-12, atol=0)


@pytest.mark.integration
def test_mode_b_branches_match_independent_mode_a_runs(two_branch_run):
    """Each MODE B per-branch curve equals its independent MODE A run."""
    outdir = two_branch_run["outdir"]
    manifest = json.loads((outdir / "manifest.json").read_text())
    by_value = {
        b["fdha_calc_params"]["r_threshold_km"]: b["global_index"]
        for b in manifest["branches"]
    }
    np.testing.assert_allclose(
        _branch_rates(outdir, by_value[T1]), two_branch_run["rate_t1"],
        rtol=0, atol=0,
    )
    np.testing.assert_allclose(
        _branch_rates(outdir, by_value[T2]), two_branch_run["rate_t2"],
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

    # Threshold branch IDs appear in the composed branch path; the chosen
    # value is recorded under fdha_calc_params.
    ids = {b["fdha_branch_id"].split("|")[0] for b in manifest["branches"]}
    assert ids == {"RT_A", "RT_B"}
    values = {b["fdha_calc_params"]["r_threshold_km"] for b in manifest["branches"]}
    assert values == {T1, T2}
    # Model slots keep listing only real model classes.
    for b in manifest["branches"]:
        assert "calc_r_threshold" not in b["fdha_models"]
        assert set(b["fdha_models"]) == {
            "primary_surf_rup", "primary_surf_displ",
            "secondary_surf_rup", "secondary_surf_displ",
        }


# --------------------------------------------------- conflict rule (driver)

@pytest.mark.integration
def test_conflict_scalar_plus_branch_set_fails_at_startup(tmp_path):
    ini = _make_curve_job(
        tmp_path / "job", scalar="0.1", branches=[("RT_V", "0.2", "1.0")]
    )
    lt = FdhaLogicTree.from_ini(str(ini))
    with pytest.raises(ConfigurationError, match="mutually exclusive"):
        lt.run(outdir=tmp_path / "out")
    # No branch may have been executed.
    assert not (tmp_path / "out" / "hazard_curves").exists()
