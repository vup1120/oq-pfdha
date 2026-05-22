from __future__ import annotations

import pytest

from openquake.fdha.logic_tree.source_model_lt import (
    IMPLICIT_BRANCH_ID,
    SourceModelBranch,
    SourceModelLogicTreeError,
    apply_realization_to_sources,
    load_source_model_branches,
    parse_source_model_logic_tree,
)


def _write_minimal_source_xml(path):
    path.write_text(
        '<?xml version="1.0" encoding="utf-8"?>\n'
        '<nrml xmlns="http://openquake.org/xmlns/nrml/0.5"\n'
        '      xmlns:gml="http://www.opengis.net/gml">\n'
        '  <sourceModel name="dummy">\n'
        '    <sourceGroup tectonicRegion="Active Shallow Crust">\n'
        '      <pointSource id="1" name="p" tectonicRegion="Active Shallow Crust"/>\n'
        '    </sourceGroup>\n'
        '  </sourceModel>\n'
        '</nrml>\n'
    )


def _write_simple_fault_source(path, *, sid="3", aValue=4.2, bValue=0.9,
                               minMag=6.5, maxMag=7.5, dip=30.0):
    path.write_text(
        f'<?xml version="1.0" encoding="utf-8"?>\n'
        f'<nrml xmlns="http://openquake.org/xmlns/nrml/0.5"\n'
        f'      xmlns:gml="http://www.opengis.net/gml">\n'
        f'  <sourceModel name="sm">\n'
        f'    <sourceGroup tectonicRegion="Active Shallow Crust">\n'
        f'      <simpleFaultSource id="{sid}" name="F"\n'
        f'                         tectonicRegion="Active Shallow Crust">\n'
        f'        <simpleFaultGeometry>\n'
        f'          <gml:LineString>\n'
        f'            <gml:posList>16.18 39.69 16.14 39.58</gml:posList>\n'
        f'          </gml:LineString>\n'
        f'          <dip>{dip}</dip>\n'
        f'          <upperSeismoDepth>0.0</upperSeismoDepth>\n'
        f'          <lowerSeismoDepth>15.0</lowerSeismoDepth>\n'
        f'        </simpleFaultGeometry>\n'
        f'        <magScaleRel>WC1994</magScaleRel>\n'
        f'        <ruptAspectRatio>2.0</ruptAspectRatio>\n'
        f'        <truncGutenbergRichterMFD aValue="{aValue}" bValue="{bValue}"\n'
        f'                                  maxMag="{maxMag}" minMag="{minMag}"/>\n'
        f'        <rake>90.0</rake>\n'
        f'      </simpleFaultSource>\n'
        f'    </sourceGroup>\n'
        f'  </sourceModel>\n'
        f'</nrml>\n'
    )


def _write_smlt(path, branches, namespace="http://openquake.org/xmlns/nrml/0.4",
                utype="sourceModel"):
    branch_xml = "\n".join(
        f'      <logicTreeBranch branchID="{bid}">'
        f'\n        <uncertaintyModel>{p}</uncertaintyModel>'
        f'\n        <uncertaintyWeight>{w}</uncertaintyWeight>'
        f'\n      </logicTreeBranch>'
        for bid, p, w in branches
    )
    path.write_text(
        f'<?xml version="1.0" encoding="UTF-8"?>\n'
        f'<nrml xmlns="{namespace}">\n'
        f'  <logicTree logicTreeID="lt0">\n'
        f'    <logicTreeBranchSet branchSetID="bs1" uncertaintyType="{utype}">\n'
        f'{branch_xml}\n'
        f'    </logicTreeBranchSet>\n'
        f'  </logicTree>\n'
        f'</nrml>\n'
    )


# --------------------------------------------------- sourceModel-only tests


def test_parses_two_branches_with_relative_paths(tmp_path):
    sm_a = tmp_path / "source_model_a.xml"
    sm_b = tmp_path / "source_model_b.xml"
    _write_minimal_source_xml(sm_a)
    _write_minimal_source_xml(sm_b)

    smlt = tmp_path / "smlt.xml"
    _write_smlt(smlt, [
        ("sm_a", "source_model_a.xml", 0.4),
        ("sm_b", "source_model_b.xml", 0.6),
    ])

    branches = parse_source_model_logic_tree(smlt)
    assert len(branches) == 2
    assert [b.weight for b in branches] == [0.4, 0.6]
    for b in branches:
        # Multi-uncertainty design: branch_id is a path, source-model branch
        # ID is recorded in metadata.
        assert b.metadata["source_model_branch_id"] in {"sm_a", "sm_b"}
        assert b.uncertainties == ()
        for fp in b.metadata["files"]:
            assert fp.startswith(str(tmp_path.resolve()))


def test_absolute_branch_path_rejected(tmp_path):
    sm = tmp_path / "abs.xml"
    _write_minimal_source_xml(sm)
    smlt = tmp_path / "smlt.xml"
    _write_smlt(smlt, [("sm0", str(sm.resolve()), 1.0)])

    with pytest.raises(SourceModelLogicTreeError, match="must be a relative path"):
        parse_source_model_logic_tree(smlt)


def test_parses_multifile_branch(tmp_path):
    sm_a = tmp_path / "a.xml"
    sm_b = tmp_path / "b.xml"
    sm_c = tmp_path / "c.xml"
    for p in (sm_a, sm_b, sm_c):
        _write_minimal_source_xml(p)
    smlt = tmp_path / "smlt.xml"
    _write_smlt(smlt, [
        ("combo", "a.xml b.xml", 0.5),
        ("solo", "c.xml", 0.5),
    ])
    branches = parse_source_model_logic_tree(smlt)
    assert [b.metadata["source_model_branch_id"] for b in branches] == ["combo", "solo"]
    assert branches[0].metadata["files"] == [str(sm_a.resolve()), str(sm_b.resolve())]
    assert branches[1].metadata["files"] == [str(sm_c.resolve())]


def test_weights_must_sum_to_one(tmp_path):
    sm_a = tmp_path / "a.xml"
    sm_b = tmp_path / "b.xml"
    _write_minimal_source_xml(sm_a)
    _write_minimal_source_xml(sm_b)
    smlt = tmp_path / "smlt.xml"
    _write_smlt(smlt, [
        ("sm_a", "a.xml", 0.5),
        ("sm_b", "b.xml", 0.4),
    ])
    with pytest.raises(SourceModelLogicTreeError, match="weights sum"):
        parse_source_model_logic_tree(smlt)


def test_missing_source_model_file_raises(tmp_path):
    smlt = tmp_path / "smlt.xml"
    _write_smlt(smlt, [("sm_a", "does_not_exist.xml", 1.0)])
    with pytest.raises(SourceModelLogicTreeError, match="No such file|not found"):
        parse_source_model_logic_tree(smlt)


def test_missing_smlt_file_raises(tmp_path):
    with pytest.raises(SourceModelLogicTreeError, match="not found"):
        parse_source_model_logic_tree(tmp_path / "nope.xml")


# ------------------------------------------------- multi-uncertainty tests


def test_bgr_relative_uncertainty_modifies_b_value(tmp_path):
    sm = tmp_path / "sm.xml"
    _write_simple_fault_source(sm, sid="3", bValue=0.9)

    smlt = tmp_path / "smlt.xml"
    smlt.write_text(
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">\n'
        '  <logicTree logicTreeID="lt0">\n'
        '    <logicTreeBranchSet branchSetID="bs1" uncertaintyType="sourceModel">\n'
        '      <logicTreeBranch branchID="sm0">\n'
        '        <uncertaintyModel>sm.xml</uncertaintyModel>\n'
        '        <uncertaintyWeight>1.0</uncertaintyWeight>\n'
        '      </logicTreeBranch>\n'
        '    </logicTreeBranchSet>\n'
        '    <logicTreeBranchSet branchSetID="bs2"\n'
        '                        uncertaintyType="bGRRelative"\n'
        '                        applyToSources="3">\n'
        '      <logicTreeBranch branchID="bg0">\n'
        '        <uncertaintyModel>0.0</uncertaintyModel>\n'
        '        <uncertaintyWeight>0.5</uncertaintyWeight>\n'
        '      </logicTreeBranch>\n'
        '      <logicTreeBranch branchID="bg1">\n'
        '        <uncertaintyModel>0.1</uncertaintyModel>\n'
        '        <uncertaintyWeight>0.5</uncertaintyWeight>\n'
        '      </logicTreeBranch>\n'
        '    </logicTreeBranchSet>\n'
        '  </logicTree>\n'
        '</nrml>\n'
    )

    branches = parse_source_model_logic_tree(smlt)
    assert len(branches) == 2
    assert [b.weight for b in branches] == pytest.approx([0.5, 0.5])
    types = [u.uncertainty_type for b in branches for u in b.uncertainties]
    assert types == ["bGRRelative", "bGRRelative"]
    values = sorted(u.value for b in branches for u in b.uncertainties)
    assert values == [0.0, 0.1]

    from openquake.fdha.calc.utils.parsing import parse_source_model_faults
    base = parse_source_model_faults(
        [str(sm)], hdf5path="", rupture_mesh_spacing=2.0, width_of_mfd_bin=0.1,
    )
    b_vals = []
    for br in branches:
        modified = apply_realization_to_sources(br, base)
        b_vals.append(modified["3"].mfd.b_val)
    # base 0.9 + relative {0.0, 0.1} -> {0.9, 1.0}
    assert sorted(b_vals) == pytest.approx([0.9, 1.0])
    # Ensure original sources are untouched.
    assert base["3"].mfd.b_val == pytest.approx(0.9)


def test_max_mag_gr_relative_modifies_max_mag(tmp_path):
    sm = tmp_path / "sm.xml"
    _write_simple_fault_source(sm, sid="3", maxMag=7.5)
    smlt = tmp_path / "smlt.xml"
    smlt.write_text(
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">\n'
        '  <logicTree logicTreeID="lt0">\n'
        '    <logicTreeBranchSet branchSetID="bs1" uncertaintyType="sourceModel">\n'
        '      <logicTreeBranch branchID="sm0">\n'
        '        <uncertaintyModel>sm.xml</uncertaintyModel>\n'
        '        <uncertaintyWeight>1.0</uncertaintyWeight>\n'
        '      </logicTreeBranch>\n'
        '    </logicTreeBranchSet>\n'
        '    <logicTreeBranchSet branchSetID="bs2"\n'
        '                        uncertaintyType="maxMagGRRelative"\n'
        '                        applyToSources="3">\n'
        '      <logicTreeBranch branchID="m0">\n'
        '        <uncertaintyModel>0.0</uncertaintyModel>\n'
        '        <uncertaintyWeight>0.5</uncertaintyWeight>\n'
        '      </logicTreeBranch>\n'
        '      <logicTreeBranch branchID="m1">\n'
        '        <uncertaintyModel>0.2</uncertaintyModel>\n'
        '        <uncertaintyWeight>0.5</uncertaintyWeight>\n'
        '      </logicTreeBranch>\n'
        '    </logicTreeBranchSet>\n'
        '  </logicTree>\n'
        '</nrml>\n'
    )

    branches = parse_source_model_logic_tree(smlt)
    from openquake.fdha.calc.utils.parsing import parse_source_model_faults
    base = parse_source_model_faults(
        [str(sm)], hdf5path="", rupture_mesh_spacing=2.0, width_of_mfd_bin=0.1,
    )
    modified_max = sorted(
        apply_realization_to_sources(br, base)["3"].mfd.max_mag for br in branches
    )
    assert modified_max == pytest.approx([7.5, 7.7])


def test_simple_fault_dip_relative_modifies_dip(tmp_path):
    sm = tmp_path / "sm.xml"
    _write_simple_fault_source(sm, sid="3", dip=30.0)
    smlt = tmp_path / "smlt.xml"
    smlt.write_text(
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">\n'
        '  <logicTree logicTreeID="lt0">\n'
        '    <logicTreeBranchSet branchSetID="bs1" uncertaintyType="sourceModel">\n'
        '      <logicTreeBranch branchID="sm0">\n'
        '        <uncertaintyModel>sm.xml</uncertaintyModel>\n'
        '        <uncertaintyWeight>1.0</uncertaintyWeight>\n'
        '      </logicTreeBranch>\n'
        '    </logicTreeBranchSet>\n'
        '    <logicTreeBranchSet branchSetID="bs2"\n'
        '                        uncertaintyType="simpleFaultDipRelative"\n'
        '                        applyToSources="3">\n'
        '      <logicTreeBranch branchID="d0">\n'
        '        <uncertaintyModel>0.0</uncertaintyModel>\n'
        '        <uncertaintyWeight>0.5</uncertaintyWeight>\n'
        '      </logicTreeBranch>\n'
        '      <logicTreeBranch branchID="d1">\n'
        '        <uncertaintyModel>10.0</uncertaintyModel>\n'
        '        <uncertaintyWeight>0.5</uncertaintyWeight>\n'
        '      </logicTreeBranch>\n'
        '    </logicTreeBranchSet>\n'
        '  </logicTree>\n'
        '</nrml>\n'
    )
    branches = parse_source_model_logic_tree(smlt)
    from openquake.fdha.calc.utils.parsing import parse_source_model_faults
    base = parse_source_model_faults(
        [str(sm)], hdf5path="", rupture_mesh_spacing=2.0, width_of_mfd_bin=0.1,
    )
    dips = sorted(
        apply_realization_to_sources(br, base)["3"].dip for br in branches
    )
    assert dips == pytest.approx([30.0, 40.0])


def test_combined_uncertainties_cartesian_product(tmp_path):
    """sourceModel × bGRRelative × maxMagGRRelative -> 1 × 2 × 2 = 4 paths."""
    sm = tmp_path / "sm.xml"
    _write_simple_fault_source(sm, sid="3", bValue=0.9, maxMag=7.5)
    smlt = tmp_path / "smlt.xml"
    smlt.write_text(
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">\n'
        '  <logicTree logicTreeID="lt0">\n'
        '    <logicTreeBranchSet branchSetID="bs1" uncertaintyType="sourceModel">\n'
        '      <logicTreeBranch branchID="sm0">\n'
        '        <uncertaintyModel>sm.xml</uncertaintyModel>\n'
        '        <uncertaintyWeight>1.0</uncertaintyWeight>\n'
        '      </logicTreeBranch>\n'
        '    </logicTreeBranchSet>\n'
        '    <logicTreeBranchSet branchSetID="bs2"\n'
        '                        uncertaintyType="bGRRelative"\n'
        '                        applyToSources="3">\n'
        '      <logicTreeBranch branchID="bg0">\n'
        '        <uncertaintyModel>0.0</uncertaintyModel>\n'
        '        <uncertaintyWeight>0.5</uncertaintyWeight>\n'
        '      </logicTreeBranch>\n'
        '      <logicTreeBranch branchID="bg1">\n'
        '        <uncertaintyModel>0.1</uncertaintyModel>\n'
        '        <uncertaintyWeight>0.5</uncertaintyWeight>\n'
        '      </logicTreeBranch>\n'
        '    </logicTreeBranchSet>\n'
        '    <logicTreeBranchSet branchSetID="bs3"\n'
        '                        uncertaintyType="maxMagGRRelative"\n'
        '                        applyToSources="3">\n'
        '      <logicTreeBranch branchID="m0">\n'
        '        <uncertaintyModel>0.0</uncertaintyModel>\n'
        '        <uncertaintyWeight>0.7</uncertaintyWeight>\n'
        '      </logicTreeBranch>\n'
        '      <logicTreeBranch branchID="m1">\n'
        '        <uncertaintyModel>0.2</uncertaintyModel>\n'
        '        <uncertaintyWeight>0.3</uncertaintyWeight>\n'
        '      </logicTreeBranch>\n'
        '    </logicTreeBranchSet>\n'
        '  </logicTree>\n'
        '</nrml>\n'
    )
    branches = parse_source_model_logic_tree(smlt)
    assert len(branches) == 4
    assert sum(b.weight for b in branches) == pytest.approx(1.0)
    expected_weights = sorted([0.5 * 0.7, 0.5 * 0.3, 0.5 * 0.7, 0.5 * 0.3])
    assert sorted(b.weight for b in branches) == pytest.approx(expected_weights)
    # Each branch has both uncertainties recorded.
    for b in branches:
        types = sorted(u.uncertainty_type for u in b.uncertainties)
        assert types == ["bGRRelative", "maxMagGRRelative"]


def test_apply_to_sources_filters_modification(tmp_path):
    """An uncertainty applyToSources=A must not modify source B."""
    # Two-fault model so we can exercise filtering.
    path = tmp_path / "sm.xml"
    path.write_text(
        '<?xml version="1.0" encoding="utf-8"?>\n'
        '<nrml xmlns="http://openquake.org/xmlns/nrml/0.5"\n'
        '      xmlns:gml="http://www.opengis.net/gml">\n'
        '  <sourceModel name="sm">\n'
        '    <sourceGroup tectonicRegion="Active Shallow Crust">\n'
        '      <simpleFaultSource id="A" name="A" tectonicRegion="Active Shallow Crust">\n'
        '        <simpleFaultGeometry>\n'
        '          <gml:LineString><gml:posList>16.18 39.69 16.14 39.58</gml:posList></gml:LineString>\n'
        '          <dip>30.0</dip><upperSeismoDepth>0.0</upperSeismoDepth><lowerSeismoDepth>15.0</lowerSeismoDepth>\n'
        '        </simpleFaultGeometry>\n'
        '        <magScaleRel>WC1994</magScaleRel><ruptAspectRatio>2.0</ruptAspectRatio>\n'
        '        <truncGutenbergRichterMFD aValue="4.2" bValue="0.9" maxMag="7.5" minMag="6.5"/>\n'
        '        <rake>90.0</rake>\n'
        '      </simpleFaultSource>\n'
        '      <simpleFaultSource id="B" name="B" tectonicRegion="Active Shallow Crust">\n'
        '        <simpleFaultGeometry>\n'
        '          <gml:LineString><gml:posList>16.30 39.70 16.30 39.60</gml:posList></gml:LineString>\n'
        '          <dip>60.0</dip><upperSeismoDepth>0.0</upperSeismoDepth><lowerSeismoDepth>20.0</lowerSeismoDepth>\n'
        '        </simpleFaultGeometry>\n'
        '        <magScaleRel>WC1994</magScaleRel><ruptAspectRatio>2.0</ruptAspectRatio>\n'
        '        <truncGutenbergRichterMFD aValue="4.0" bValue="0.85" maxMag="7.0" minMag="6.5"/>\n'
        '        <rake>90.0</rake>\n'
        '      </simpleFaultSource>\n'
        '    </sourceGroup>\n'
        '  </sourceModel>\n'
        '</nrml>\n'
    )
    smlt = tmp_path / "smlt.xml"
    smlt.write_text(
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">\n'
        '  <logicTree logicTreeID="lt0">\n'
        '    <logicTreeBranchSet branchSetID="bs1" uncertaintyType="sourceModel">\n'
        '      <logicTreeBranch branchID="sm0"><uncertaintyModel>sm.xml</uncertaintyModel><uncertaintyWeight>1.0</uncertaintyWeight></logicTreeBranch>\n'
        '    </logicTreeBranchSet>\n'
        '    <logicTreeBranchSet branchSetID="bs2" uncertaintyType="bGRRelative" applyToSources="A">\n'
        '      <logicTreeBranch branchID="bg0"><uncertaintyModel>0.5</uncertaintyModel><uncertaintyWeight>1.0</uncertaintyWeight></logicTreeBranch>\n'
        '    </logicTreeBranchSet>\n'
        '  </logicTree>\n'
        '</nrml>\n'
    )
    branches = parse_source_model_logic_tree(smlt)
    assert len(branches) == 1
    assert branches[0].uncertainties[0].apply_to_sources == ("A",)

    from openquake.fdha.calc.utils.parsing import parse_source_model_faults
    base = parse_source_model_faults(
        [str(path)], hdf5path="", rupture_mesh_spacing=2.0, width_of_mfd_bin=0.1,
    )
    modified = apply_realization_to_sources(branches[0], base)
    # Source A had +0.5 added; source B is untouched.
    assert modified["A"].mfd.b_val == pytest.approx(1.4)
    assert modified["B"].mfd.b_val == pytest.approx(0.85)
    # Original dict not mutated.
    assert base["A"].mfd.b_val == pytest.approx(0.9)


def test_apply_to_branches_only_modifies_selected_source_model(tmp_path):
    """``applyToBranches`` should restrict modifications to specific
    sourceModel branches; the other branch must remain unmodified."""
    sm_a = tmp_path / "sm_a.xml"
    sm_b = tmp_path / "sm_b.xml"
    _write_simple_fault_source(sm_a, sid="3", bValue=0.9)
    _write_simple_fault_source(sm_b, sid="3", bValue=0.9)
    smlt = tmp_path / "smlt.xml"
    smlt.write_text(
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<nrml xmlns="http://openquake.org/xmlns/nrml/0.4">\n'
        '  <logicTree logicTreeID="lt0">\n'
        '    <logicTreeBranchSet branchSetID="bs1" uncertaintyType="sourceModel">\n'
        '      <logicTreeBranch branchID="sm_a"><uncertaintyModel>sm_a.xml</uncertaintyModel><uncertaintyWeight>0.5</uncertaintyWeight></logicTreeBranch>\n'
        '      <logicTreeBranch branchID="sm_b"><uncertaintyModel>sm_b.xml</uncertaintyModel><uncertaintyWeight>0.5</uncertaintyWeight></logicTreeBranch>\n'
        '    </logicTreeBranchSet>\n'
        '    <logicTreeBranchSet branchSetID="bs2" uncertaintyType="bGRRelative"\n'
        '                        applyToSources="3" applyToBranches="sm_a">\n'
        '      <logicTreeBranch branchID="bg0"><uncertaintyModel>0.5</uncertaintyModel><uncertaintyWeight>1.0</uncertaintyWeight></logicTreeBranch>\n'
        '    </logicTreeBranchSet>\n'
        '  </logicTree>\n'
        '</nrml>\n'
    )
    branches = parse_source_model_logic_tree(smlt)
    # Two realisations: sm_a + bg0 (modified), sm_b (untouched).
    by_sm = {b.metadata["source_model_branch_id"]: b for b in branches}
    assert set(by_sm) == {"sm_a", "sm_b"}
    assert len(by_sm["sm_a"].uncertainties) == 1
    assert by_sm["sm_b"].uncertainties == ()
    assert by_sm["sm_a"].weight == pytest.approx(0.5)
    assert by_sm["sm_b"].weight == pytest.approx(0.5)


# ------------------------------------------------------ normalisation tests


def test_load_source_model_branches_implicit_when_no_lt_file(tmp_path):
    sm = tmp_path / "single.xml"
    _write_minimal_source_xml(sm)
    config = {"calculation": {"source_model_file": "single.xml"}}
    branches = load_source_model_branches(config, tmp_path)
    assert len(branches) == 1
    assert branches[0].branch_id == IMPLICIT_BRANCH_ID
    assert branches[0].weight == 1.0
    assert branches[0].uncertainties == ()


def test_load_source_model_branches_uses_lt_file_when_present(tmp_path):
    sm_a = tmp_path / "a.xml"
    sm_b = tmp_path / "b.xml"
    _write_minimal_source_xml(sm_a)
    _write_minimal_source_xml(sm_b)
    smlt = tmp_path / "smlt.xml"
    _write_smlt(smlt, [
        ("sm_a", "a.xml", 0.5),
        ("sm_b", "b.xml", 0.5),
    ])
    config = {
        "calculation": {
            "source_model_file": "ignored_when_lt_set.xml",
            "source_model_logic_tree_file": "smlt.xml",
        }
    }
    branches = load_source_model_branches(config, tmp_path)
    assert sum(b.weight for b in branches) == pytest.approx(1.0)


def test_load_source_model_branches_errors_when_nothing_set(tmp_path):
    config = {"calculation": {}}
    with pytest.raises(SourceModelLogicTreeError):
        load_source_model_branches(config, tmp_path)
