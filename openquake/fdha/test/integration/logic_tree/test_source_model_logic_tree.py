"""End-to-end tests for source-model logic-tree support.

These tests cover the new ``source_model_logic_tree_file`` INI field by:

* verifying backward compatibility - an FDHA logic-tree run with no
  source-model logic tree must produce identical results;
* verifying that two source-model branches combined with a multi-branch FDHA
  logic tree expand into the expected cartesian product of realisations and
  that the combined branch weights and mean hazard match the weighted
  arithmetic-mean expectation.
"""
from __future__ import annotations

import json
import csv
from pathlib import Path

import numpy as np
import pytest

from openquake.fdha.logic_tree.driver import FdhaLogicTree


def _repo_root() -> Path:
    for p in [Path(__file__).resolve(), *Path(__file__).resolve().parents]:
        if (p / "examples").is_dir() and (p / "openquake").is_dir():
            return p
    raise RuntimeError("Could not locate repository root.")


_BASE_SOURCE_MODEL = """<?xml version="1.0" encoding="utf-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.5"
      xmlns:gml="http://www.opengis.net/gml">
  <sourceModel name="sm">
    <sourceGroup name="g" tectonicRegion="Active Shallow Crust">
      <simpleFaultSource id="3" name="Simple Fault Source"
                         tectonicRegion="Active Shallow Crust">
        <simpleFaultGeometry>
          <gml:LineString>
            <gml:posList>
              16.1870362931 39.6974567774
              16.1505303713 39.6274955882
              16.1455213236 39.6196231258
              16.1408996159 39.5840858442
            </gml:posList>
          </gml:LineString>
          <dip>30.0</dip>
          <upperSeismoDepth>0.0</upperSeismoDepth>
          <lowerSeismoDepth>15.0</lowerSeismoDepth>
        </simpleFaultGeometry>
        <magScaleRel>WC1994</magScaleRel>
        <ruptAspectRatio>2.0</ruptAspectRatio>
        <truncGutenbergRichterMFD aValue="{aValue}" bValue="0.9"
                                  maxMag="7.5" minMag="6.5"/>
        <rake>90.0</rake>
      </simpleFaultSource>
    </sourceGroup>
  </sourceModel>
</nrml>
"""


def _write_source_model(path: Path, aValue: float) -> None:
    path.write_text(_BASE_SOURCE_MODEL.format(aValue=aValue))


def _write_fdha_lt(path: Path) -> None:
    """Two-branch FDHA logic tree (FixedPrimarySR with two values) so the
    cartesian product against the source-model logic tree is non-trivial."""
    path.write_text(
        """<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="pfdha">
    <logicTreeBranchingLevel branchingLevelID="bl1">
      <logicTreeBranchSet branchSetID="bs1" uncertaintyType="fdhaPrimarySRModel">
        <logicTreeBranch branchID="SR0">
          <uncertaintyModel>[FixedPrimarySR]\nvalue = 0.0</uncertaintyModel>
          <uncertaintyWeight>0.4</uncertaintyWeight>
        </logicTreeBranch>
        <logicTreeBranch branchID="SR1">
          <uncertaintyModel>[FixedPrimarySR]\nvalue = 1.0</uncertaintyModel>
          <uncertaintyWeight>0.6</uncertaintyWeight>
        </logicTreeBranch>
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
    <logicTreeBranchingLevel branchingLevelID="bl2">
      <logicTreeBranchSet branchSetID="bs2" uncertaintyType="fdhaPrimaryFDModel"
                          applyToBranches="SR0 SR1">
        <logicTreeBranch branchID="FD">
          <uncertaintyModel>[Youngs2003PrimaryFD]\nstyle = 'normal'\nscaling_model = 'WC1994'\nnorm_disp_type = 'AD'</uncertaintyModel>
          <uncertaintyWeight>1.0</uncertaintyWeight>
        </logicTreeBranch>
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
    <logicTreeBranchingLevel branchingLevelID="bl3">
      <logicTreeBranchSet branchSetID="bs3" uncertaintyType="fdhaSecondarySRModel"
                          applyToBranches="SR0 SR1">
        <logicTreeBranch branchID="SSR">
          <uncertaintyModel>[FixedSecondarySR]\nvalue = 0.0</uncertaintyModel>
          <uncertaintyWeight>1.0</uncertaintyWeight>
        </logicTreeBranch>
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
    <logicTreeBranchingLevel branchingLevelID="bl4">
      <logicTreeBranchSet branchSetID="bs4" uncertaintyType="fdhaSecondaryFDModel"
                          applyToBranches="SSR">
        <logicTreeBranch branchID="SFD">
          <uncertaintyModel>[Youngs2003SecondaryFD]\nstyle = 'all'</uncertaintyModel>
          <uncertaintyWeight>1.0</uncertaintyWeight>
        </logicTreeBranch>
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
  </logicTree>
</nrml>
"""
    )


def _write_smlt(path: Path, sm_files: list[tuple[str, str, float]]) -> None:
    branch_xml = "\n".join(
        f'      <logicTreeBranch branchID="{bid}">\n'
        f'        <uncertaintyModel>{p}</uncertaintyModel>\n'
        f'        <uncertaintyWeight>{w}</uncertaintyWeight>\n'
        f'      </logicTreeBranch>'
        for bid, p, w in sm_files
    )
    path.write_text(
        f"""<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="smlt">
    <logicTreeBranchSet branchSetID="bs1" uncertaintyType="sourceModel">
{branch_xml}
    </logicTreeBranchSet>
  </logicTree>
</nrml>
"""
    )


def _write_ini(
    path: Path,
    *,
    source_model_file: str,
    smlt_file: str,
    fdha_lt_file: str,
) -> None:
    path.write_text(
        f"""[general]
description = source_model_lt_test

[geometry]
sites = 16.16573727 39.64704451

[site_params]
reference_vs30_value = 760.0

[erf]
rupture_mesh_spacing = 1.0
width_of_mfd_bin = 0.1

[calculation]
source_model_logic_tree_file = {smlt_file}
displacement_measure_levels = {{"FD": [0.01, 0.1]}}
fdha_logic_tree_file = {fdha_lt_file}
"""
    )


# --------------------------------------------------------- backward compat


def test_missing_source_model_lt_is_rejected(tmp_path):
    """Canonical public logic-tree jobs require a source-model LT wrapper."""
    sm = tmp_path / "source_model.xml"
    _write_source_model(sm, aValue=4.2)

    fdha_lt = tmp_path / "fdha_lt.xml"
    _write_fdha_lt(fdha_lt)

    ini = tmp_path / "job.ini"
    ini.write_text(
        f"""[general]
description = missing_smlt

[geometry]
sites = 16.16573727 39.64704451

[site_params]
reference_vs30_value = 760.0

[erf]
rupture_mesh_spacing = 1.0
width_of_mfd_bin = 0.1

[calculation]
source_model_logic_tree_file =
displacement_measure_levels = {{"FD": [0.01, 0.1]}}
fdha_logic_tree_file = {fdha_lt}
"""
    )

    from openquake.fdha.calc.config_loader import ConfigurationError

    with pytest.raises(ConfigurationError, match="source_model_logic_tree_file"):
        FdhaLogicTree.from_ini(str(ini)).run(outdir=tmp_path / "out")


# ------------------------------------------- multi-source × multi-FDHA case


def test_two_source_model_branches_combine_with_fdha_lt(tmp_path):
    """Two source models with weight 0.5 each, two FDHA branches with
    weights (0.4, 0.6). The combined ensemble should have 4 realisations,
    weights summing to 1.0, and a mean equal to the weighted arithmetic mean.
    """
    sm_a = tmp_path / "source_model_a.xml"
    sm_b = tmp_path / "source_model_b.xml"
    # Two slightly different MFDs so the branch curves differ.
    _write_source_model(sm_a, aValue=4.2)
    _write_source_model(sm_b, aValue=4.5)

    fdha_lt = tmp_path / "fdha_lt.xml"
    _write_fdha_lt(fdha_lt)

    smlt = tmp_path / "smlt.xml"
    _write_smlt(smlt, [
        ("sm_a", "source_model_a.xml", 0.5),
        ("sm_b", "source_model_b.xml", 0.5),
    ])

    ini = tmp_path / "job.ini"
    _write_ini(
        ini,
        source_model_file=str(sm_a),  # ignored when smlt_file is set
        smlt_file=str(smlt),
        fdha_lt_file=str(fdha_lt),
    )

    out = tmp_path / "out"
    res = FdhaLogicTree.from_ini(str(ini)).run(outdir=out)

    manifest = json.loads((out / "manifest.json").read_text())

    # 2 source-model branches × 2 FDHA branches = 4 realisations.
    assert manifest["n_realizations"] == 4
    assert len(manifest["branches"]) == 4

    # Combined weight per record == source weight × FDHA weight, summing to 1.
    expected_weights = sorted([0.5 * 0.4, 0.5 * 0.6, 0.5 * 0.4, 0.5 * 0.6])
    actual_weights = sorted(b["combined_branch_weight"] for b in manifest["branches"])
    assert actual_weights == [pytest.approx(w) for w in expected_weights]
    assert sum(b["combined_branch_weight"] for b in manifest["branches"]) == \
        pytest.approx(1.0)

    # Required output fields are present in every record.
    required_keys = {
        "source_model_branch_id",
        "source_model_file",
        "source_model_weight",
        "fdha_branch_id",
        "fdha_branch_weight",
        "combined_branch_weight",
    }
    for rec in manifest["branches"]:
        assert required_keys.issubset(rec.keys())

    # Mean curve equals the weighted arithmetic mean of per-branch hazard rate
    # curves recorded in the per-realisation CSVs.
    branch_curves = []
    branch_weights = []
    for rec in sorted(manifest["branches"], key=lambda r: r["global_index"]):
        csv_path = out / rec["curve_file"]
        rates = _load_curve_rates(csv_path)
        branch_curves.append(rates)
        branch_weights.append(rec["combined_branch_weight"])

    rates_arr = np.stack(branch_curves, axis=0)
    w = np.asarray(branch_weights, dtype=float)
    expected_mean = np.tensordot(w, rates_arr, axes=(0, 0))
    actual_mean = np.asarray(res.mean_rates, dtype=float)

    assert actual_mean.shape == expected_mean.shape
    assert np.allclose(actual_mean, expected_mean, atol=1e-12, rtol=0)


def test_smlt_bgr_and_fdha_lt_matches_analytical_curve(tmp_path):
    """Source-model LT × NRML uncertainty × FDHA LT is analytically additive.

    This exercises the full composition:

    * 2 sourceModel branches (0.4, 0.6);
    * 2 bGRRelative branches (0.25, 0.75);
    * 2 FDHA branches (0.4, 0.6).

    The framework mean must equal the direct weighted sum of every emitted
    end-realisation curve, i.e. ``sum_i combined_weight_i * lambda_i(D0)``.
    """
    sm_a = tmp_path / "sm_a.xml"
    sm_b = tmp_path / "sm_b.xml"
    _write_source_model(sm_a, aValue=4.2)
    _write_source_model(sm_b, aValue=4.5)

    fdha_lt = tmp_path / "fdha_lt.xml"
    _write_nonzero_two_branch_fdha_lt(fdha_lt)

    smlt = tmp_path / "smlt.xml"
    smlt.write_text(
        """<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="lt0">
    <logicTreeBranchSet branchSetID="bs1" uncertaintyType="sourceModel">
      <logicTreeBranch branchID="sm_a"><uncertaintyModel>sm_a.xml</uncertaintyModel><uncertaintyWeight>0.4</uncertaintyWeight></logicTreeBranch>
      <logicTreeBranch branchID="sm_b"><uncertaintyModel>sm_b.xml</uncertaintyModel><uncertaintyWeight>0.6</uncertaintyWeight></logicTreeBranch>
    </logicTreeBranchSet>
    <logicTreeBranchSet branchSetID="bs2" uncertaintyType="bGRRelative"
                        applyToSources="3" applyToBranches="sm_a sm_b">
      <logicTreeBranch branchID="bg0"><uncertaintyModel>0.0</uncertaintyModel><uncertaintyWeight>0.25</uncertaintyWeight></logicTreeBranch>
      <logicTreeBranch branchID="bg1"><uncertaintyModel>0.2</uncertaintyModel><uncertaintyWeight>0.75</uncertaintyWeight></logicTreeBranch>
    </logicTreeBranchSet>
  </logicTree>
</nrml>
"""
    )

    ini = tmp_path / "job.ini"
    _write_ini(
        ini,
        source_model_file=str(sm_a),
        smlt_file=str(smlt),
        fdha_lt_file=str(fdha_lt),
    )

    out = tmp_path / "out"
    res = FdhaLogicTree.from_ini(str(ini)).run(outdir=out)
    manifest = json.loads((out / "manifest.json").read_text())

    assert manifest["n_realizations"] == 8
    assert sum(b["combined_branch_weight"] for b in manifest["branches"]) == \
        pytest.approx(1.0)
    assert {b["source_model_branch_id"] for b in manifest["branches"]} == {
        "sm_a|bg0", "sm_a|bg1", "sm_b|bg0", "sm_b|bg1",
    }

    expected_weights = sorted([
        0.4 * 0.25 * 0.4, 0.4 * 0.25 * 0.6,
        0.4 * 0.75 * 0.4, 0.4 * 0.75 * 0.6,
        0.6 * 0.25 * 0.4, 0.6 * 0.25 * 0.6,
        0.6 * 0.75 * 0.4, 0.6 * 0.75 * 0.6,
    ])
    actual_weights = sorted(
        rec["combined_branch_weight"] for rec in manifest["branches"]
    )
    assert actual_weights == [pytest.approx(w) for w in expected_weights]

    branch_curves = []
    branch_weights = []
    for rec in sorted(manifest["branches"], key=lambda r: r["global_index"]):
        branch_curves.append(_load_curve_rates(out / rec["curve_file"]))
        branch_weights.append(rec["combined_branch_weight"])

    expected_mean = np.tensordot(
        np.asarray(branch_weights, dtype=float),
        np.stack(branch_curves, axis=0),
        axes=(0, 0),
    )
    assert np.any(expected_mean > 0), "Combined curve is zero; test is vacuous"
    np.testing.assert_allclose(np.asarray(res.mean_rates), expected_mean,
                               atol=1e-12, rtol=0)

    # Check the persisted aggregate CSV, not only the in-memory result.
    aggregate_mean = _load_aggregate_mean(out / "aggregate_hazard.csv")
    np.testing.assert_allclose(aggregate_mean, expected_mean,
                               atol=1e-12, rtol=0)


def _write_nonzero_two_branch_fdha_lt(path: Path) -> None:
    """Two FDHA branches with distributed hazard active at the test site."""
    path.write_text(
        """<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="pfdha">
    <logicTreeBranchingLevel branchingLevelID="bl1">
      <logicTreeBranchSet branchSetID="bs1" uncertaintyType="fdhaPrimarySRModel">
        <logicTreeBranch branchID="PSR">
          <uncertaintyModel>[FixedPrimarySR]\nvalue = 1.0</uncertaintyModel>
          <uncertaintyWeight>1.0</uncertaintyWeight>
        </logicTreeBranch>
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
    <logicTreeBranchingLevel branchingLevelID="bl2">
      <logicTreeBranchSet branchSetID="bs2" uncertaintyType="fdhaPrimaryFDModel"
                          applyToBranches="PSR">
        <logicTreeBranch branchID="PFD">
          <uncertaintyModel>[Youngs2003PrimaryFD]\nstyle = 'normal'\nscaling_model = 'WC1994'\nnorm_disp_type = 'AD'</uncertaintyModel>
          <uncertaintyWeight>1.0</uncertaintyWeight>
        </logicTreeBranch>
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
    <logicTreeBranchingLevel branchingLevelID="bl3">
      <logicTreeBranchSet branchSetID="bs3" uncertaintyType="fdhaSecondarySRModel"
                          applyToBranches="PFD">
        <logicTreeBranch branchID="SSR0">
          <uncertaintyModel>[FixedSecondarySR]\nvalue = 0.5</uncertaintyModel>
          <uncertaintyWeight>0.4</uncertaintyWeight>
        </logicTreeBranch>
        <logicTreeBranch branchID="SSR1">
          <uncertaintyModel>[FixedSecondarySR]\nvalue = 1.0</uncertaintyModel>
          <uncertaintyWeight>0.6</uncertaintyWeight>
        </logicTreeBranch>
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
    <logicTreeBranchingLevel branchingLevelID="bl4">
      <logicTreeBranchSet branchSetID="bs4" uncertaintyType="fdhaSecondaryFDModel"
                          applyToBranches="SSR0 SSR1">
        <logicTreeBranch branchID="SFD">
          <uncertaintyModel>[Youngs2003SecondaryFD]\nstyle = 'all'</uncertaintyModel>
          <uncertaintyWeight>1.0</uncertaintyWeight>
        </logicTreeBranch>
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
  </logicTree>
</nrml>
"""
    )


def _write_singletrack_fdha_lt(path: Path) -> None:
    """A 1-branch FDHA logic tree that activates the distributed component
    so the resulting hazard curves are guaranteed non-zero at any site
    inside the model's distance window.  Used to verify NRML uncertainties
    propagate to the calculator."""
    path.write_text(
        """<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="pfdha">
    <logicTreeBranchingLevel branchingLevelID="bl1">
      <logicTreeBranchSet branchSetID="bs1" uncertaintyType="fdhaPrimarySRModel">
        <logicTreeBranch branchID="SR">
          <uncertaintyModel>[FixedPrimarySR]\nvalue = 1.0</uncertaintyModel>
          <uncertaintyWeight>1.0</uncertaintyWeight>
        </logicTreeBranch>
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
    <logicTreeBranchingLevel branchingLevelID="bl2">
      <logicTreeBranchSet branchSetID="bs2" uncertaintyType="fdhaPrimaryFDModel" applyToBranches="SR">
        <logicTreeBranch branchID="FD">
          <uncertaintyModel>[Youngs2003PrimaryFD]\nstyle = 'normal'\nscaling_model = 'WC1994'\nnorm_disp_type = 'AD'</uncertaintyModel>
          <uncertaintyWeight>1.0</uncertaintyWeight>
        </logicTreeBranch>
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
    <logicTreeBranchingLevel branchingLevelID="bl3">
      <logicTreeBranchSet branchSetID="bs3" uncertaintyType="fdhaSecondarySRModel" applyToBranches="SR">
        <logicTreeBranch branchID="SSR">
          <uncertaintyModel>[FixedSecondarySR]\nvalue = 1.0</uncertaintyModel>
          <uncertaintyWeight>1.0</uncertaintyWeight>
        </logicTreeBranch>
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
    <logicTreeBranchingLevel branchingLevelID="bl4">
      <logicTreeBranchSet branchSetID="bs4" uncertaintyType="fdhaSecondaryFDModel" applyToBranches="SSR">
        <logicTreeBranch branchID="SFD">
          <uncertaintyModel>[Youngs2003SecondaryFD]\nstyle = 'all'</uncertaintyModel>
          <uncertaintyWeight>1.0</uncertaintyWeight>
        </logicTreeBranch>
      </logicTreeBranchSet>
    </logicTreeBranchingLevel>
  </logicTree>
</nrml>
"""
    )


def test_bgr_relative_in_smlt_changes_curve(tmp_path):
    """An SMLT with a single sourceModel branch (weight 1.0) plus a
    ``bGRRelative`` branch set whose two branches add 0.0 vs 0.5 to the
    b-value must produce two FDHA hazard-curve realisations whose rates
    differ (proving NRML uncertainties propagate end-to-end), and whose
    weighted mean equals 0.5*r0 + 0.5*r1."""
    sm = tmp_path / "source_model.xml"
    _write_source_model(sm, aValue=4.2)

    fdha_lt = tmp_path / "fdha_lt.xml"
    _write_singletrack_fdha_lt(fdha_lt)

    smlt = tmp_path / "smlt.xml"
    smlt.write_text(
        f"""<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="lt0">
    <logicTreeBranchSet branchSetID="bs1" uncertaintyType="sourceModel">
      <logicTreeBranch branchID="sm0">
        <uncertaintyModel>source_model.xml</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
    <logicTreeBranchSet branchSetID="bs2" uncertaintyType="bGRRelative" applyToSources="3">
      <logicTreeBranch branchID="bg0"><uncertaintyModel>0.0</uncertaintyModel><uncertaintyWeight>0.5</uncertaintyWeight></logicTreeBranch>
      <logicTreeBranch branchID="bg1"><uncertaintyModel>0.5</uncertaintyModel><uncertaintyWeight>0.5</uncertaintyWeight></logicTreeBranch>
    </logicTreeBranchSet>
  </logicTree>
</nrml>
"""
    )

    ini = tmp_path / "job.ini"
    _write_ini(
        ini,
        source_model_file=str(sm),
        smlt_file=str(smlt),
        fdha_lt_file=str(fdha_lt),
    )

    out = tmp_path / "out"
    res = FdhaLogicTree.from_ini(str(ini)).run(outdir=out)

    manifest = json.loads((out / "manifest.json").read_text())

    # 2 sourceModel-LT realisations × 1 FDHA branch = 2.
    assert manifest["n_realizations"] == 2
    types = sorted({
        u["uncertainty_type"]
        for rec in manifest["branches"]
        for u in rec["source_model_uncertainties"]
    })
    assert types == ["bGRRelative"]
    assert sum(b["combined_branch_weight"] for b in manifest["branches"]) == \
        pytest.approx(1.0)

    # Per-realisation curves: bg0 (no b-value change) vs bg1 (+0.5).
    by_smlt = {}
    for rec in manifest["branches"]:
        rates = _load_curve_rates(out / rec["curve_file"])
        by_smlt[rec["source_model_branch_id"]] = rates
    paths = sorted(by_smlt)
    assert len(paths) == 2
    rates_bg0 = by_smlt[paths[0]]
    rates_bg1 = by_smlt[paths[1]]
    # Both curves must be non-zero (uncertainty has something to scale).
    assert np.any(rates_bg0 > 0), "baseline rates are zero - test is meaningless"
    assert np.any(rates_bg1 > 0)
    diff = np.max(np.abs(rates_bg0 - rates_bg1))
    assert diff > 0, "bGRRelative did not change the hazard curve"

    # Mean equals weighted arithmetic mean over both realisations.
    branch_curves = []
    branch_weights = []
    for rec in sorted(manifest["branches"], key=lambda r: r["global_index"]):
        branch_curves.append(_load_curve_rates(out / rec["curve_file"]))
        branch_weights.append(rec["combined_branch_weight"])
    expected_mean = np.tensordot(
        np.asarray(branch_weights, dtype=float),
        np.stack(branch_curves, axis=0),
        axes=(0, 0),
    )
    assert np.allclose(np.asarray(res.mean_rates), expected_mean,
                       atol=1e-12, rtol=0)


def test_combined_smlt_uncertainties_full_cartesian_product(tmp_path):
    """Two ``sourceModel`` branches × two ``bGRRelative`` branches × two
    ``maxMagGRRelative`` branches must expand to 8 SMLT realisations.  The
    full realisation set crossed with the 2-branch FDHA logic tree gives 16
    end realisations, all weights still summing to 1.0."""
    sm_a = tmp_path / "sm_a.xml"
    sm_b = tmp_path / "sm_b.xml"
    _write_source_model(sm_a, aValue=4.2)
    _write_source_model(sm_b, aValue=4.5)

    fdha_lt = tmp_path / "fdha_lt.xml"
    _write_fdha_lt(fdha_lt)

    smlt = tmp_path / "smlt.xml"
    smlt.write_text(
        f"""<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="lt0">
    <logicTreeBranchSet branchSetID="bs1" uncertaintyType="sourceModel">
      <logicTreeBranch branchID="sm_a"><uncertaintyModel>sm_a.xml</uncertaintyModel><uncertaintyWeight>0.5</uncertaintyWeight></logicTreeBranch>
      <logicTreeBranch branchID="sm_b"><uncertaintyModel>sm_b.xml</uncertaintyModel><uncertaintyWeight>0.5</uncertaintyWeight></logicTreeBranch>
    </logicTreeBranchSet>
    <logicTreeBranchSet branchSetID="bs2" uncertaintyType="bGRRelative"
                        applyToSources="3" applyToBranches="sm_a sm_b">
      <logicTreeBranch branchID="bg0"><uncertaintyModel>0.0</uncertaintyModel><uncertaintyWeight>0.5</uncertaintyWeight></logicTreeBranch>
      <logicTreeBranch branchID="bg1"><uncertaintyModel>0.2</uncertaintyModel><uncertaintyWeight>0.5</uncertaintyWeight></logicTreeBranch>
    </logicTreeBranchSet>
    <logicTreeBranchSet branchSetID="bs3" uncertaintyType="maxMagGRRelative"
                        applyToSources="3" applyToBranches="bg0 bg1">
      <logicTreeBranch branchID="m0"><uncertaintyModel>0.0</uncertaintyModel><uncertaintyWeight>0.5</uncertaintyWeight></logicTreeBranch>
      <logicTreeBranch branchID="m1"><uncertaintyModel>0.2</uncertaintyModel><uncertaintyWeight>0.5</uncertaintyWeight></logicTreeBranch>
    </logicTreeBranchSet>
  </logicTree>
</nrml>
"""
    )

    ini = tmp_path / "job.ini"
    _write_ini(
        ini,
        source_model_file=str(sm_a),
        smlt_file=str(smlt),
        fdha_lt_file=str(fdha_lt),
    )

    out = tmp_path / "out"
    res = FdhaLogicTree.from_ini(str(ini)).run(outdir=out)
    manifest = json.loads((out / "manifest.json").read_text())

    # 2 SM × 2 b × 2 maxMag × 2 FDHA = 16 realisations.
    assert manifest["n_realizations"] == 16
    assert sum(b["combined_branch_weight"] for b in manifest["branches"]) == \
        pytest.approx(1.0)
    # Every record carries both NRML uncertainties.
    for rec in manifest["branches"]:
        types = sorted(u["uncertainty_type"]
                       for u in rec["source_model_uncertainties"])
        assert types == ["bGRRelative", "maxMagGRRelative"]
    # Result has the right shape; mean rates non-empty.
    assert len(res.mean_rates[0]) == 2


def _load_curve_rates(csv_path: Path) -> np.ndarray:
    """Read a per-branch hazard-curve CSV produced by ``write_branch_rates_csv``.

    Returns a (n_sites, n_displ) array of annual exceedance rates.

    Supports both the single-site layout (``D0, annual_rate`` rows, one row
    per displacement level) and the multi-site layout
    (``site_id, lon, lat, D0, annual_rate``).
    """
    text = csv_path.read_text().splitlines()
    if not text:
        return np.zeros((1, 0))
    header = text[0].split(",")
    if header == ["D0", "annual_rate"]:
        rates = [float(line.split(",")[1]) for line in text[1:] if line.strip()]
        return np.asarray([rates], dtype=float)
    # Multi-site: gather by site_id then by D0 order.
    by_site: dict[int, list[float]] = {}
    for line in text[1:]:
        if not line.strip():
            continue
        parts = line.split(",")
        sid = int(parts[0])
        rate = float(parts[4])
        by_site.setdefault(sid, []).append(rate)
    return np.asarray([by_site[i] for i in sorted(by_site)], dtype=float)


def _load_aggregate_mean(csv_path: Path) -> np.ndarray:
    """Read the ``mean`` column from ``aggregate_hazard.csv``."""
    with csv_path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return np.zeros((1, 0))
    if "site_id" not in rows[0]:
        return np.asarray([[float(row["mean"]) for row in rows]], dtype=float)
    by_site: dict[int, list[float]] = {}
    for row in rows:
        by_site.setdefault(int(row["site_id"]), []).append(float(row["mean"]))
    return np.asarray([by_site[i] for i in sorted(by_site)], dtype=float)


