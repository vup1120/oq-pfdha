"""End-to-end logic-tree runs verified against closed-form expected values.

These are the FDHA analogs of the OpenQuake engine's full-enumeration QA
tests (``CollapseTestCase`` in ``openquake/hazardlib/tests/lt_test.py`` and
the ``qa_tests_data/logictree`` cases): real hazard calculations whose
aggregate statistics are checked against values derived *outside* the
logic-tree framework.

Two independent anchors are used:

1. ``FixedPrimarySR`` makes the hazard integral exactly linear in the fixed
   surface-rupture probability p: every branch curve must equal
   ``p * base_curve``, so the weighted mean and every weighted fractile of
   the ensemble have closed forms computable by hand.

2. An ``abGRAbsolute`` SMLT branch must produce the same curve as a plain
   run whose source-model XML has the modified (a, b) baked in - a
   non-circular check that realisation modifications reach the calculator
   with the same numbers the XML parser would deliver.
"""
from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import numpy as np
import pytest

from openquake.fdha.logic_tree.driver import FdhaLogicTree

DEMO_DIR = Path(__file__).resolve().parents[3] / "demo" / "hazard_curve"

FD_CHAIN = {
    "primary_fd": (
        "[Youngs2003PrimaryFD]\nstyle = 'normal'\n"
        "scaling_model = 'WC1994'\nnorm_disp_type = 'AD'"
    ),
    "secondary_sr": "[FixedSecondarySR]\nvalue = 0.0",
    "secondary_fd": "[Youngs2003SecondaryFD]\nstyle = 'all'",
}


def _fdha_lt_xml(primary_sr_branches: list[tuple[str, str, float]]) -> str:
    """FDHA tree: variable primary-SR branch set + fixed single-branch rest."""
    sr_branches = "\n".join(
        f'        <logicTreeBranch branchID="{bid}">\n'
        f"          <uncertaintyModel>{model}</uncertaintyModel>\n"
        f"          <uncertaintyWeight>{w}</uncertaintyWeight>\n"
        f"        </logicTreeBranch>"
        for bid, model, w in primary_sr_branches
    )
    sr_ids = " ".join(bid for bid, _, _ in primary_sr_branches)
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="pfdha">
    <logicTreeBranchSet branchSetID="bs1" uncertaintyType="fdhaPrimarySRModel">
{sr_branches}
    </logicTreeBranchSet>
    <logicTreeBranchSet branchSetID="bs2" uncertaintyType="fdhaPrimaryFDModel" applyToBranches="{sr_ids}">
      <logicTreeBranch branchID="FD">
        <uncertaintyModel>{FD_CHAIN["primary_fd"]}</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
    <logicTreeBranchSet branchSetID="bs3" uncertaintyType="fdhaSecondarySRModel" applyToBranches="FD">
      <logicTreeBranch branchID="SSR">
        <uncertaintyModel>{FD_CHAIN["secondary_sr"]}</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
    <logicTreeBranchSet branchSetID="bs4" uncertaintyType="fdhaSecondaryFDModel" applyToBranches="SSR">
      <logicTreeBranch branchID="SFD">
        <uncertaintyModel>{FD_CHAIN["secondary_fd"]}</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
  </logicTree>
</nrml>
"""


def _smlt_xml(source_model_name: str, extra_branchsets: str = "") -> str:
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="smlt">
    <logicTreeBranchSet branchSetID="bs_sm" uncertaintyType="sourceModel">
      <logicTreeBranch branchID="SM">
        <uncertaintyModel>{source_model_name}</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
{extra_branchsets}
  </logicTree>
</nrml>
"""


def _write_job(job_dir: Path, *, smlt_body: str, fdha_lt_body: str,
               source_xml: Path) -> Path:
    job_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source_xml, job_dir / "source_model.xml")
    (job_dir / "smlt.xml").write_text(smlt_body)
    (job_dir / "fdha_lt.xml").write_text(fdha_lt_body)
    ini = job_dir / "job.ini"
    ini.write_text(
        "[general]\ndescription = closed_form_numerical\n\n"
        # On-trace site (a fault-trace vertex, same as the demo job):
        # principal hazard is zero at off-trace sites, which would make
        # every closed-form check vacuous.
        "[geometry]\nsites = 16.1455213236 39.6196231258\n\n"
        "[site_params]\nreference_vs30_value = 760.0\n\n"
        "[erf]\nrupture_mesh_spacing = 1.0\nwidth_of_mfd_bin = 0.1\n\n"
        "[calculation]\n"
        "source_model_logic_tree_file = smlt.xml\n"
        'displacement_measure_levels = {"FD": [0.01, 0.1]}\n'
        "fdha_logic_tree_file = fdha_lt.xml\n"
    )
    return ini


def _run(ini: Path):
    res = FdhaLogicTree.from_ini(str(ini)).run(outdir=ini.parent / "out")
    manifest = json.loads((ini.parent / "out" / "manifest.json").read_text())
    return res, manifest


def _branch_curve(outdir: Path, curve_file: str) -> np.ndarray:
    """Total annual rate column of a single-site branch CSV."""
    return np.loadtxt(outdir / curve_file, delimiter=",", skiprows=1,
                      usecols=1, ndmin=1)


def _read_aggregate_csv(outdir: Path) -> dict[str, np.ndarray]:
    with (outdir / "aggregate_hazard.csv").open() as f:
        rows = list(csv.DictReader(f))
    return {
        name: np.array([float(r[name]) for r in rows])
        for name in rows[0]
    }


# ---------------------------------------------------------------------------
# 1. FixedPrimarySR: linearity + closed-form mean and fractiles
# ---------------------------------------------------------------------------

def test_fixed_sr_ensemble_matches_closed_form(tmp_path):
    """Three FixedPrimarySR branches p = 0.2 / 0.5 / 1.0, w = 0.2 / 0.3 / 0.5.

    Physics: the principal hazard integral is linear in P(SR), so
    ``curve(p) = p * curve(1)``. With the secondary chain fixed at zero the
    distributed contribution vanishes. On such a proportional ensemble the
    weighted statistics have closed forms:

    * mean  = (0.2*0.2 + 0.3*0.5 + 0.5*1.0) * base = 0.69 * base
    * sorted branch values [.2b, .5b, b] have cumulative weights
      [.2, .5, 1.0], so the weighted fractile at q is
      ``max(q, 0.2) * base`` for every q <= 1 (linear CDF interpolation).
    """
    branches = [
        ("SR02", "[FixedPrimarySR]\nvalue = 0.2", 0.2),
        ("SR05", "[FixedPrimarySR]\nvalue = 0.5", 0.3),
        ("SR10", "[FixedPrimarySR]\nvalue = 1.0", 0.5),
    ]
    ini = _write_job(
        tmp_path,
        smlt_body=_smlt_xml("source_model.xml"),
        fdha_lt_body=_fdha_lt_xml(branches),
        source_xml=DEMO_DIR / "source_model.xml",
    )
    res, manifest = _run(ini)
    outdir = Path(res.outdir)

    assert manifest["n_realizations"] == 3
    recs = manifest["branches"]
    assert sum(r["combined_branch_weight"] for r in recs) == pytest.approx(1.0)

    # Map realisations to their SR value via the (distinct) branch weights.
    by_weight = {round(r["combined_branch_weight"], 6): r for r in recs}
    curve = {
        0.2: _branch_curve(outdir, by_weight[0.2]["curve_file"]),
        0.5: _branch_curve(outdir, by_weight[0.3]["curve_file"]),
        1.0: _branch_curve(outdir, by_weight[0.5]["curve_file"]),
    }
    base = curve[1.0]
    assert np.all(base > 0), "base hazard curve is zero - test is vacuous"

    # --- physics: hazard is linear in the fixed P(SR)
    np.testing.assert_allclose(curve[0.2], 0.2 * base, rtol=1e-9)
    np.testing.assert_allclose(curve[0.5], 0.5 * base, rtol=1e-9)

    # --- closed-form weighted mean
    mean = np.asarray(res.mean_rates, dtype=float)[0]
    np.testing.assert_allclose(mean, 0.69 * base, rtol=1e-9)

    # --- closed-form weighted fractiles: max(q, 0.2) * base
    for q, expected_factor in [(0.05, 0.2), (0.16, 0.2), (0.5, 0.5),
                               (0.84, 0.84), (0.95, 0.95)]:
        got = np.asarray(res.fractiles[q], dtype=float)[0]
        np.testing.assert_allclose(
            got, expected_factor * base, rtol=1e-9,
            err_msg=f"fractile q={q}",
        )

    # --- with the secondary chain at zero the split is all-principal
    mean_principal = np.asarray(res.mean_rates_principal, dtype=float)[0]
    mean_distributed = np.asarray(res.mean_rates_distributed, dtype=float)[0]
    np.testing.assert_allclose(mean_principal, mean, rtol=1e-12)
    np.testing.assert_allclose(mean_distributed, 0.0, atol=1e-300)

    # --- the aggregate CSV must reproduce the in-memory result exactly
    cols = _read_aggregate_csv(outdir)
    np.testing.assert_allclose(cols["mean"], mean, rtol=1e-12)
    np.testing.assert_allclose(cols["mean_principal"], mean_principal,
                               rtol=1e-12)
    np.testing.assert_allclose(cols["mean_distributed"], 0.0, atol=1e-300)
    for q in (0.05, 0.16, 0.5, 0.84, 0.95):
        label = [c for c in cols if c.endswith(str(q))]
        assert label, f"missing fractile column for q={q}"
        np.testing.assert_allclose(
            cols[label[0]], np.asarray(res.fractiles[q], dtype=float)[0],
            rtol=1e-12,
        )


# ---------------------------------------------------------------------------
# 2. abGRAbsolute SMLT branch vs baked-in source model (non-circular)
# ---------------------------------------------------------------------------

AB_BASE = (4.2, 0.9)   # the demo source model's own (a, b)
AB_MOD = (4.5, 1.0)


def _bake_ab(source_xml: Path, out_xml: Path, a: float, b: float) -> None:
    text = source_xml.read_text()
    assert 'aValue="4.2000000E+00" bValue="9.0000000E-01"' in text
    out_xml.write_text(
        text.replace(
            'aValue="4.2000000E+00" bValue="9.0000000E-01"',
            f'aValue="{a}" bValue="{b}"',
        )
    )


def test_smlt_abgr_absolute_equals_baked_in_model(tmp_path):
    """SMLT: sourceModel x abGRAbsolute{(4.2,0.9) w=.35, (4.5,1.0) w=.65}.

    The identity branch must reproduce a plain run of the untouched XML,
    the modified branch must reproduce a plain run of an XML with (a', b')
    baked in, and the logic-tree mean must equal the closed-form
    combination ``0.35 * c_base + 0.65 * c_mod`` of those two *independent*
    reference runs.
    """
    single_sr = [("SR10", "[FixedPrimarySR]\nvalue = 1.0", 1.0)]

    # Reference run C: untouched demo source model, trivial tree.
    ini_c = _write_job(
        tmp_path / "ref_base",
        smlt_body=_smlt_xml("source_model.xml"),
        fdha_lt_body=_fdha_lt_xml(single_sr),
        source_xml=DEMO_DIR / "source_model.xml",
    )
    res_c, _ = _run(ini_c)
    c_base = np.asarray(res_c.mean_rates, dtype=float)[0]

    # Reference run B: (a', b') baked into the XML, trivial tree.
    baked = tmp_path / "baked_source.xml"
    _bake_ab(DEMO_DIR / "source_model.xml", baked, *AB_MOD)
    ini_b = _write_job(
        tmp_path / "ref_mod",
        smlt_body=_smlt_xml("source_model.xml"),
        fdha_lt_body=_fdha_lt_xml(single_sr),
        source_xml=baked,
    )
    res_b, _ = _run(ini_b)
    c_mod = np.asarray(res_b.mean_rates, dtype=float)[0]

    assert np.all(c_base > 0)
    assert not np.allclose(c_mod, c_base), \
        "modified (a, b) did not change the hazard - vacuous comparison"

    # Run A: the same two alternatives expressed as an abGRAbsolute
    # branch set inside the source-model logic tree.
    ab_branchset = f"""    <logicTreeBranchSet branchSetID="bs_ab"
                        uncertaintyType="abGRAbsolute"
                        applyToSources="3">
      <logicTreeBranch branchID="AB_BASE">
        <uncertaintyModel>{AB_BASE[0]} {AB_BASE[1]}</uncertaintyModel>
        <uncertaintyWeight>0.35</uncertaintyWeight>
      </logicTreeBranch>
      <logicTreeBranch branchID="AB_MOD">
        <uncertaintyModel>{AB_MOD[0]} {AB_MOD[1]}</uncertaintyModel>
        <uncertaintyWeight>0.65</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>"""
    ini_a = _write_job(
        tmp_path / "lt_ab",
        smlt_body=_smlt_xml("source_model.xml", ab_branchset),
        fdha_lt_body=_fdha_lt_xml(single_sr),
        source_xml=DEMO_DIR / "source_model.xml",
    )
    res_a, manifest = _run(ini_a)
    outdir_a = Path(res_a.outdir)

    assert manifest["n_realizations"] == 2
    recs = {r["source_model_branch_id"]: r for r in manifest["branches"]}
    assert set(recs) == {"SM|AB_BASE", "SM|AB_MOD"}
    assert recs["SM|AB_BASE"]["combined_branch_weight"] == pytest.approx(0.35)
    assert recs["SM|AB_MOD"]["combined_branch_weight"] == pytest.approx(0.65)
    utypes = [
        u["uncertainty_type"]
        for u in recs["SM|AB_MOD"]["source_model_uncertainties"]
    ]
    assert utypes == ["abGRAbsolute"]

    lt_base = _branch_curve(outdir_a, recs["SM|AB_BASE"]["curve_file"])
    lt_mod = _branch_curve(outdir_a, recs["SM|AB_MOD"]["curve_file"])

    # Identity branch == plain parse; modified branch == baked-in parse.
    np.testing.assert_allclose(lt_base, c_base, rtol=1e-9)
    np.testing.assert_allclose(lt_mod, c_mod, rtol=1e-9)

    # Logic-tree mean == closed-form combination of the two reference runs.
    mean_a = np.asarray(res_a.mean_rates, dtype=float)[0]
    np.testing.assert_allclose(mean_a, 0.35 * c_base + 0.65 * c_mod,
                               rtol=1e-9)
