# -*- coding: utf-8 -*-
"""
Curve-mode logic tree with per-source model assignment (applyToSources).

Independent sources ADD hazard. When branch selections differ per source,
each end-branch must run against its own source group only, and the final
mean must be the SUM over source groups of each group's weighted LT mean -
NOT a pooled weighted mean of every-branch-times-every-source runs (which
both averages away half the hazard and applies each group's models to the
other group's sources).

Regression for the bug where curve mode passed the full ``fault_sources``
dict to every end-branch (map mode already filtered correctly).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from openquake.fdha.logic_tree.driver import FdhaLogicTree
from openquake.fdha.logic_tree.aggregation import weighted_fractiles

pytestmark = pytest.mark.integration


# Two parallel ~16.5 km normal faults; the site sits between them,
# ~5.5 km from each trace, so both contribute distributed hazard.
_SOURCE_TEMPLATE = """\
            <characteristicFaultSource id="{sid}" name="fault_{sid}"
                                       tectonicRegion="Active Shallow Crust">
                <incrementalMFD binWidth="0.1" minMag="6.5">
                    <occurRates>1.0e-4 5.0e-5 2.0e-5</occurRates>
                </incrementalMFD>
                <rake>-90.0</rake>
                <surface>
                    <simpleFaultGeometry>
                        <gml:LineString>
                            <gml:posList>13.00 {lat} 13.20 {lat}</gml:posList>
                        </gml:LineString>
                        <dip>60.0</dip>
                        <upperSeismoDepth>0.0</upperSeismoDepth>
                        <lowerSeismoDepth>15.0</lowerSeismoDepth>
                    </simpleFaultGeometry>
                </surface>
            </characteristicFaultSource>
"""

_MODEL_TEMPLATE = """<?xml version="1.0" encoding="utf-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.5"
      xmlns:gml="http://www.opengis.net/gml">
    <sourceModel name="two normal faults">
        <sourceGroup name="g1" tectonicRegion="Active Shallow Crust">
{sources}
        </sourceGroup>
    </sourceModel>
</nrml>
"""

_SFD_YOUNGS_85 = "[Youngs2003SecondaryFD]\npercentile = 85"
_SFD_YOUNGS_95 = "[Youngs2003SecondaryFD]\npercentile = 95"
_SFD_PETERSEN = "[Petersen2011SecondaryFD]\npixel_size = 25"

_COMMON_LEVELS = """\
  <logicTreeBranchSet branchSetID="bs1" uncertaintyType="fdhaPrimarySRModel">
    <logicTreeBranch branchID="PSR">
      <uncertaintyModel><![CDATA[[Youngs2003PrimarySR]
style = all]]></uncertaintyModel>
      <uncertaintyWeight>1.0</uncertaintyWeight>
    </logicTreeBranch>
  </logicTreeBranchSet>
  <logicTreeBranchSet branchSetID="bs2" uncertaintyType="fdhaPrimaryFDModel">
    <logicTreeBranch branchID="PFD">
      <uncertaintyModel><![CDATA[[Youngs2003PrimaryFD]
style = normal
norm_disp_type = AD]]></uncertaintyModel>
      <uncertaintyWeight>1.0</uncertaintyWeight>
    </logicTreeBranch>
  </logicTreeBranchSet>
  <logicTreeBranchSet branchSetID="bs3" uncertaintyType="fdhaSecondarySRModel">
    <logicTreeBranch branchID="SSR">
      <uncertaintyModel><![CDATA[[Youngs2003SecondarySR]
version = 3]]></uncertaintyModel>
      <uncertaintyWeight>1.0</uncertaintyWeight>
    </logicTreeBranch>
  </logicTreeBranchSet>
"""


def _sfd_level(bl_id, bs_id, branches, apply_to=None):
    ats = f' applyToSources="{apply_to}"' if apply_to else ""
    body = "".join(
        f"""      <logicTreeBranch branchID="{bid}">
        <uncertaintyModel><![CDATA[{model}]]></uncertaintyModel>
        <uncertaintyWeight>{weight}</uncertaintyWeight>
      </logicTreeBranch>\n"""
        for bid, model, weight in branches
    )
    return (
        f'    <logicTreeBranchSet branchSetID="{bs_id}" '
        f'uncertaintyType="fdhaSecondaryFDModel"{ats}>\n'
        f"{body}"
        f"    </logicTreeBranchSet>\n"
    )


def _write_lt(path: Path, sfd_levels: str):
    path.write_text(
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">\n'
        '  <logicTree logicTreeID="lt">\n'
        + _COMMON_LEVELS + sfd_levels
        + "  </logicTree>\n</nrml>\n"
    )


def _write_job(workdir: Path, name: str, source_bodies: str, sfd_levels: str):
    workdir.mkdir(parents=True, exist_ok=True)
    sm = workdir / "source_model.xml"
    sm.write_text(_MODEL_TEMPLATE.format(sources=source_bodies))
    smlt = workdir / "smlt.xml"
    smlt.write_text(
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">\n'
        '  <logicTree logicTreeID="smlt">\n'
        '    <logicTreeBranchSet branchSetID="bs_sm" uncertaintyType="sourceModel">\n'
        '      <logicTreeBranch branchID="sm0">\n'
        f'        <uncertaintyModel>{sm.name}</uncertaintyModel>\n'
        '        <uncertaintyWeight>1.0</uncertaintyWeight>\n'
        '      </logicTreeBranch>\n'
        '    </logicTreeBranchSet>\n'
        '  </logicTree>\n</nrml>\n'
    )
    lt = workdir / "fdha_lt.xml"
    _write_lt(lt, sfd_levels)
    ini = workdir / "job.ini"
    ini.write_text(
        "[general]\n"
        f"description = {name}\n"
        "calculation_mode = fdha_classical\n\n"
        "[geometry]\n"
        "sites = 13.10 42.05\n\n"
        "[erf]\n"
        "rupture_mesh_spacing = 2.0\n"
        "width_of_mfd_bin = 0.1\n\n"
        "[calculation]\n"
        "investigation_time = 1.0\n"
        'displacement_measure_levels = {"FD": [0.01, 0.1]}\n'
        "r_threshold_km = 0.1\n"
        f"fdha_logic_tree_file = {lt.name}\n"
        f"source_model_logic_tree_file = {smlt.name}\n\n"
        "[output]\n"
        "mean = true\n"
        "quantiles = 0.5\n"
    )
    return ini


def _run(ini: Path):
    res = FdhaLogicTree(str(ini)).run(outdir=ini.parent / "out")
    return res


SRC1 = _SOURCE_TEMPLATE.format(sid="1", lat="42.00")
SRC2 = _SOURCE_TEMPLATE.format(sid="2", lat="42.10")


def test_per_source_assignment_sums_sources(tmp_path):
    # --- job under test: two sources, per-source secondary-FD selections;
    # source 1 additionally carries a weighted two-branch choice.
    lt_levels = (
        _sfd_level("bl4a", "bs4a",
                   [("SFD1_85", _SFD_YOUNGS_85, 0.6),
                    ("SFD1_95", _SFD_YOUNGS_95, 0.4)],
                   apply_to="1")
        + _sfd_level("bl4b", "bs4b",
                     [("SFD2_PET", _SFD_PETERSEN, 1.0)],
                     apply_to="2")
    )
    res = _run(_write_job(tmp_path / "two", "two_sources",
                          SRC1 + SRC2, lt_levels))
    mean_two = np.asarray(res.mean_rates, dtype=float)

    # --- references: single-source runs with a single, unrestricted model
    def ref(name, src_body, model):
        levels = _sfd_level("bl4", "bs4", [("SFD", model, 1.0)])
        r = _run(_write_job(tmp_path / name, name, src_body, levels))
        return np.asarray(r.mean_rates, dtype=float)

    r1_85 = ref("only1_85", SRC1, _SFD_YOUNGS_85)
    r1_95 = ref("only1_95", SRC1, _SFD_YOUNGS_95)
    r2_pet = ref("only2_pet", SRC2, _SFD_PETERSEN)

    expected_mean = 0.6 * r1_85 + 0.4 * r1_95 + r2_pet
    assert expected_mean.max() > 0.0, "reference hazard is zero - vacuous"
    np.testing.assert_allclose(mean_two, expected_mean, rtol=1e-9)

    # --- fractiles: per-group weighted fractile, summed across groups
    # (group 2 has a single branch, so its p50 equals its mean curve)
    stacked = np.stack([r1_85, r1_95], axis=0)
    p50_group1 = weighted_fractiles(stacked, np.array([0.6, 0.4]),
                                    qs=[0.5])[0.5]
    expected_p50 = p50_group1 + r2_pet
    p50_two = np.asarray(res.fractiles[0.5], dtype=float)
    np.testing.assert_allclose(p50_two, expected_p50, rtol=1e-9)


def test_same_selection_sources_unchanged(tmp_path):
    """Two same-style sources sharing one selection (the dedup path) must
    keep the historical pooled behaviour: mean == sum of both sources."""
    levels = _sfd_level("bl4", "bs4", [("SFD", _SFD_YOUNGS_85, 1.0)])
    res = _run(_write_job(tmp_path / "both", "both_same",
                          SRC1 + SRC2, levels))
    mean_both = np.asarray(res.mean_rates, dtype=float)

    r1 = np.asarray(_run(_write_job(tmp_path / "s1", "s1", SRC1, levels)
                         ).mean_rates, dtype=float)
    r2 = np.asarray(_run(_write_job(tmp_path / "s2", "s2", SRC2, levels)
                         ).mean_rates, dtype=float)
    assert (r1 + r2).max() > 0.0
    np.testing.assert_allclose(mean_both, r1 + r2, rtol=1e-9)
