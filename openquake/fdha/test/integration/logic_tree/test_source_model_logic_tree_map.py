"""End-to-end tests for ``source_model_logic_tree_file`` in hazard-map mode.

These tests guard the multi-source map path: with an SMLT configured the
driver must run map mode without raising ``NotImplementedError`` and produce
SMLT × FDHA realisations whose mean rate cube is the SMLT-weighted
arithmetic mean of the per-realisation maps.  Without an SMLT the existing
single-source map path must still produce identical outputs.
"""
from __future__ import annotations

import json
import csv
from pathlib import Path

import h5py
import numpy as np
import pytest

from openquake.fdha.calc.utils.interpolation import get_map_from_curves
from openquake.fdha.logic_tree.driver import FdhaLogicTree


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


def _write_singletrack_fdha_lt(path: Path) -> None:
    """One-branch FDHA logic tree that activates the distributed component
    so map outputs are non-zero on every site within the model's window."""
    path.write_text(
        """<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="pfdha">
    <logicTreeBranchSet branchSetID="bs1" uncertaintyType="fdhaPrimarySRModel">
      <logicTreeBranch branchID="SR">
        <uncertaintyModel>[FixedPrimarySR]\nvalue = 1.0</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
    <logicTreeBranchSet branchSetID="bs2" uncertaintyType="fdhaPrimaryFDModel" applyToBranches="SR">
      <logicTreeBranch branchID="FD">
        <uncertaintyModel>[Youngs2003PrimaryFD]\nstyle = 'normal'\nscaling_model = 'WC1994'\nnorm_disp_type = 'AD'</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
    <logicTreeBranchSet branchSetID="bs3" uncertaintyType="fdhaSecondarySRModel" applyToBranches="SR">
      <logicTreeBranch branchID="SSR">
        <uncertaintyModel>[FixedSecondarySR]\nvalue = 1.0</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
    <logicTreeBranchSet branchSetID="bs4" uncertaintyType="fdhaSecondaryFDModel" applyToBranches="SSR">
      <logicTreeBranch branchID="SFD">
        <uncertaintyModel>[Youngs2003SecondaryFD]\nstyle = 'all'</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
  </logicTree>
</nrml>
"""
    )


def _write_two_branch_fdha_lt(path: Path) -> None:
    """Two FDHA branches so map-mode SMLT aggregation uses combined weights."""
    path.write_text(
        """<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="pfdha">
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
    <logicTreeBranchSet branchSetID="bs2" uncertaintyType="fdhaPrimaryFDModel"
                        applyToBranches="SR0 SR1">
      <logicTreeBranch branchID="FD">
        <uncertaintyModel>[Youngs2003PrimaryFD]\nstyle = 'normal'\nscaling_model = 'WC1994'\nnorm_disp_type = 'AD'</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
    <logicTreeBranchSet branchSetID="bs3" uncertaintyType="fdhaSecondarySRModel"
                        applyToBranches="SR0 SR1">
      <logicTreeBranch branchID="SSR">
        <uncertaintyModel>[FixedSecondarySR]\nvalue = 0.0</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
    <logicTreeBranchSet branchSetID="bs4" uncertaintyType="fdhaSecondaryFDModel"
                        applyToBranches="SSR">
      <logicTreeBranch branchID="SFD">
        <uncertaintyModel>[Youngs2003SecondaryFD]\nstyle = 'all'</uncertaintyModel>
        <uncertaintyWeight>1.0</uncertaintyWeight>
      </logicTreeBranch>
    </logicTreeBranchSet>
  </logicTree>
</nrml>
"""
    )


def _write_smlt(path: Path, branches: list[tuple[str, str, float]]) -> None:
    branch_xml = "\n".join(
        f'      <logicTreeBranch branchID="{bid}">\n'
        f'        <uncertaintyModel>{p}</uncertaintyModel>\n'
        f'        <uncertaintyWeight>{w}</uncertaintyWeight>\n'
        f'      </logicTreeBranch>'
        for bid, p, w in branches
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


def _write_map_ini(
    path: Path,
    *,
    source_model_file: str,
    smlt_file: str,
    fdha_lt_file: str,
    return_period: float = 2475.0,
) -> None:
    path.write_text(
        f"""[general]
description = source_model_lt_map_test
calculation_mode = fdha_classical

[geometry]
region = 16.10 39.55, 16.22 39.55, 16.22 39.72, 16.10 39.72
region_grid_spacing = 0.05
max_distance_km = 15.0

[site_params]
reference_vs30_value = 760.0

[logic_tree]
number_of_logic_tree_samples = 0

[erf]
rupture_mesh_spacing = 2.0
width_of_mfd_bin = 0.1

[calculation]
source_model_logic_tree_file = {smlt_file}
investigation_time = 1.0
displacement_measure_levels = {{"FD": [0.01, 0.1, 1.0]}}
return_period = {return_period}
r_threshold_km = 0.5
principal_distance_km = 2.0
fdha_logic_tree_file = {fdha_lt_file}
"""
    )


def _read_rates_mean(out_dir: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    path = out_dir / "aggregate" / "rates_mean.h5"
    with h5py.File(path, "r") as f:
        rates = np.asarray(f["rates_mean"][()], dtype=float)
        d0 = np.asarray(f["d0"][()], dtype=float)
        lons = np.asarray(f["site_lons"][()], dtype=float)
        lats = np.asarray(f["site_lats"][()], dtype=float)
    return rates, d0, lons, lats


def _read_branch_rates_h5(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as f:
        rates = np.asarray(f["rates"][()], dtype=float)
        d0 = np.asarray(f["d0"][()], dtype=float)
        lons = np.asarray(f["site_lons"][()], dtype=float)
        lats = np.asarray(f["site_lats"][()], dtype=float)
    return rates, d0, lons, lats


def _read_displacement_mean_csv(
    path: Path,
) -> tuple[float | None, np.ndarray, np.ndarray]:
    rp = None
    with path.open() as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith("#") and "return_period" in line:
            rp = float(line.split("=")[1].strip())
    rows = [ln for ln in lines if not ln.startswith("#")]
    reader = csv.DictReader(rows)
    parsed = list(reader)
    return (
        rp,
        np.asarray([float(row["displ_mean"]) for row in parsed], dtype=float),
        np.asarray([bool(int(row["is_trace"])) for row in parsed], dtype=bool),
    )


# ---------------------------------------------------------- backward compat


def test_map_mode_without_smlt_is_rejected(tmp_path):
    """Canonical public map jobs require ``source_model_logic_tree_file``."""
    sm = tmp_path / "source_model.xml"
    _write_source_model(sm, aValue=4.2)
    fdha_lt = tmp_path / "fdha_lt.xml"
    _write_singletrack_fdha_lt(fdha_lt)

    ini = tmp_path / "job.ini"
    ini.write_text(
        f"""[general]
description = missing_smlt_map
calculation_mode = fdha_classical

[geometry]
region = 16.10 39.55, 16.22 39.55, 16.22 39.72, 16.10 39.72
region_grid_spacing = 0.05

[site_params]
reference_vs30_value = 760.0

[erf]
rupture_mesh_spacing = 2.0
width_of_mfd_bin = 0.1

[calculation]
investigation_time = 1.0
displacement_measure_levels = {{"FD": [0.01, 0.1, 1.0]}}
return_period = 2475.0
fdha_logic_tree_file = {fdha_lt}
"""
    )

    from openquake.fdha.calc.config_loader import ConfigurationError

    with pytest.raises(ConfigurationError, match="source_model_logic_tree_file"):
        FdhaLogicTree.from_ini(str(ini)).run(outdir=tmp_path / "out")


# -------------------------------------------------------- SMLT × map mode


def test_map_mode_with_smlt_no_longer_raises(tmp_path):
    """A map-mode INI with ``source_model_logic_tree_file`` configured runs
    end-to-end and yields a top-level rates_mean.h5 plus per-SMLT subdirs."""
    sm_a = tmp_path / "sm_a.xml"
    sm_b = tmp_path / "sm_b.xml"
    _write_source_model(sm_a, aValue=4.2)
    _write_source_model(sm_b, aValue=4.5)

    fdha_lt = tmp_path / "fdha_lt.xml"
    _write_singletrack_fdha_lt(fdha_lt)

    smlt = tmp_path / "smlt.xml"
    _write_smlt(smlt, [
        ("sm_a", "sm_a.xml", 0.4),
        ("sm_b", "sm_b.xml", 0.6),
    ])

    ini = tmp_path / "job.ini"
    _write_map_ini(
        ini, source_model_file=str(sm_a), smlt_file=str(smlt),
        fdha_lt_file=str(fdha_lt),
    )

    out = tmp_path / "out"
    res = FdhaLogicTree.from_ini(str(ini)).run(outdir=out)

    assert res.mode == "hazard_map"
    # Top-level aggregate exists.
    assert (out / "aggregate" / "rates_mean.h5").is_file()
    # Per-SMLT subdirectories exist.
    sub_a = out / "source_model_branches" / "00_sm_a"
    sub_b = out / "source_model_branches" / "01_sm_b"
    assert (sub_a / "aggregate" / "rates_mean.h5").is_file()
    assert (sub_b / "aggregate" / "rates_mean.h5").is_file()

    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["mode"] == "hazard_map"
    assert manifest["n_realizations"] == 2  # 2 SMLT × 1 FDHA
    assert manifest["total_combined_weight"] == pytest.approx(1.0)
    sm_branch_ids = {b["branch_id"] for b in manifest["source_model_branches"]}
    assert sm_branch_ids == {"sm_a", "sm_b"}


def test_map_smlt_mean_equals_weighted_per_branch_mean(tmp_path):
    """The combined map mean must equal sm_w[A] * mean_A + sm_w[B] * mean_B
    where mean_A and mean_B are the per-SMLT-branch mean cubes."""
    sm_a = tmp_path / "sm_a.xml"
    sm_b = tmp_path / "sm_b.xml"
    _write_source_model(sm_a, aValue=4.2)
    _write_source_model(sm_b, aValue=4.5)

    fdha_lt = tmp_path / "fdha_lt.xml"
    _write_singletrack_fdha_lt(fdha_lt)

    smlt = tmp_path / "smlt.xml"
    _write_smlt(smlt, [
        ("sm_a", "sm_a.xml", 0.4),
        ("sm_b", "sm_b.xml", 0.6),
    ])

    ini = tmp_path / "job.ini"
    _write_map_ini(
        ini, source_model_file=str(sm_a), smlt_file=str(smlt),
        fdha_lt_file=str(fdha_lt),
    )

    out = tmp_path / "out"
    FdhaLogicTree.from_ini(str(ini)).run(outdir=out)

    rates_total, _, lons_total, _ = _read_rates_mean(out)
    rates_a, _, lons_a, _ = _read_rates_mean(
        out / "source_model_branches" / "00_sm_a"
    )
    rates_b, _, lons_b, _ = _read_rates_mean(
        out / "source_model_branches" / "01_sm_b"
    )

    # All three runs must use the same site grid for the SMLT-weighted
    # arithmetic-mean equality to hold site-by-site.
    np.testing.assert_allclose(lons_total, lons_a)
    np.testing.assert_allclose(lons_total, lons_b)

    expected = 0.4 * rates_a + 0.6 * rates_b
    assert rates_total.max() > 0.0, "Combined rates are zero - test is meaningless"
    np.testing.assert_allclose(rates_total, expected, atol=1e-12, rtol=0)


def test_map_smlt_and_fdha_lt_match_analytical_rates_at_rp_10000(tmp_path):
    """Hazard-map SMLT × FDHA-LT aggregation matches analytic rate blending.

    The return period is fixed to 10000 yr to verify that the persisted
    displacement map is the inversion of the analytically blended rate cube,
    not a direct average of branch displacement maps.
    """
    sm_a = tmp_path / "sm_a.xml"
    sm_b = tmp_path / "sm_b.xml"
    _write_source_model(sm_a, aValue=4.2)
    _write_source_model(sm_b, aValue=4.5)

    fdha_lt = tmp_path / "fdha_lt.xml"
    _write_two_branch_fdha_lt(fdha_lt)

    smlt = tmp_path / "smlt.xml"
    _write_smlt(smlt, [
        ("sm_a", "sm_a.xml", 0.4),
        ("sm_b", "sm_b.xml", 0.6),
    ])

    ini = tmp_path / "job.ini"
    _write_map_ini(
        ini,
        source_model_file=str(sm_a),
        smlt_file=str(smlt),
        fdha_lt_file=str(fdha_lt),
        return_period=10000.0,
    )

    out = tmp_path / "out"
    res = FdhaLogicTree.from_ini(str(ini)).run(outdir=out)
    manifest = json.loads((out / "manifest.json").read_text())

    assert res.mode == "hazard_map"
    assert manifest["return_period"] == pytest.approx(10000.0)
    assert manifest["n_realizations"] == 4  # 2 SMLT x 2 FDHA
    assert sum(b["combined_branch_weight"] for b in manifest["branches"]) == \
        pytest.approx(1.0)

    expected_weights = sorted([
        0.4 * 0.4, 0.4 * 0.6,
        0.6 * 0.4, 0.6 * 0.6,
    ])
    actual_weights = sorted(
        rec["combined_branch_weight"] for rec in manifest["branches"]
    )
    assert actual_weights == [pytest.approx(w) for w in expected_weights]

    branch_rates = []
    branch_weights = []
    d0_ref = lons_ref = lats_ref = None
    for rec in sorted(manifest["branches"], key=lambda r: r["global_index"]):
        rates, d0, lons, lats = _read_branch_rates_h5(out / rec["branch_h5"])
        if d0_ref is None:
            d0_ref, lons_ref, lats_ref = d0, lons, lats
        else:
            np.testing.assert_allclose(d0, d0_ref)
            np.testing.assert_allclose(lons, lons_ref)
            np.testing.assert_allclose(lats, lats_ref)
        branch_rates.append(rates)
        branch_weights.append(rec["combined_branch_weight"])

    expected_rates = np.tensordot(
        np.asarray(branch_weights, dtype=float),
        np.stack(branch_rates, axis=0),
        axes=(0, 0),
    )
    assert expected_rates.max() > 0.0, "Combined rates are zero; test is vacuous"

    rates_total, d0_total, lons_total, lats_total = _read_rates_mean(out)
    np.testing.assert_allclose(d0_total, d0_ref)
    np.testing.assert_allclose(lons_total, lons_ref)
    np.testing.assert_allclose(lats_total, lats_ref)
    np.testing.assert_allclose(rates_total, expected_rates,
                               atol=1e-12, rtol=0)
    np.testing.assert_allclose(np.asarray(res.mean_rates), expected_rates,
                               atol=1e-12, rtol=0)

    manual_displ = get_map_from_curves(d0_total, expected_rates, 1.0 / 10000.0)
    rp_csv, csv_displ, csv_is_trace = _read_displacement_mean_csv(
        out / "aggregate" / "displacement_map_mean.csv"
    )
    assert rp_csv == pytest.approx(10000.0)
    assert csv_is_trace.any(), "SMLT map aggregate lost trace-site markers"
    np.testing.assert_allclose(np.asarray(res.displ_mean), manual_displ,
                               atol=1e-10, rtol=0)
    np.testing.assert_allclose(csv_displ, manual_displ, atol=1e-10, rtol=0)
    assert float(csv_displ.max()) > 0.0, (
        "All sites report zero displacement; inversion check is degenerate."
    )


def test_map_smlt_with_bgr_relative_changes_rates(tmp_path):
    """A map-mode SMLT with one ``sourceModel`` branch and a 2-branch
    ``bGRRelative`` set must produce two SMLT realisations whose per-site
    rate cubes differ (NRML uncertainties propagate to map mode)."""
    sm = tmp_path / "sm.xml"
    _write_source_model(sm, aValue=4.2)

    fdha_lt = tmp_path / "fdha_lt.xml"
    _write_singletrack_fdha_lt(fdha_lt)

    smlt = tmp_path / "smlt.xml"
    smlt.write_text(
        """<?xml version="1.0" encoding="UTF-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">
  <logicTree logicTreeID="lt0">
    <logicTreeBranchSet branchSetID="bs1" uncertaintyType="sourceModel">
      <logicTreeBranch branchID="sm0">
        <uncertaintyModel>sm.xml</uncertaintyModel>
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
    _write_map_ini(
        ini, source_model_file=str(sm), smlt_file=str(smlt),
        fdha_lt_file=str(fdha_lt),
    )

    out = tmp_path / "out"
    FdhaLogicTree.from_ini(str(ini)).run(outdir=out)

    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["n_realizations"] == 2
    types = sorted({
        u["uncertainty_type"]
        for rec in manifest["branches"]
        for u in rec["source_model_uncertainties"]
    })
    assert types == ["bGRRelative"]

    sub_dirs = sorted((out / "source_model_branches").iterdir())
    assert len(sub_dirs) == 2
    rates_first, _, _, _ = _read_rates_mean(sub_dirs[0])
    rates_second, _, _, _ = _read_rates_mean(sub_dirs[1])
    assert rates_first.max() > 0.0 and rates_second.max() > 0.0
    assert np.max(np.abs(rates_first - rates_second)) > 0.0, (
        "bGRRelative did not change map-mode hazard"
    )
