"""Weight arithmetic of the FDHA logic-tree enumeration.

Mirrors the OpenQuake engine's enumeration tests
(``BranchSetEnumerateTestCase`` and ``CompositeLogicTreeTestCase.test5`` in
``openquake/commonlib/tests/logictree_test.py`` /
``openquake/hazardlib/tests/lt_test.py``): every end-branch weight must be
the product of the weights along its path, partially-applied branch sets
(``applyToBranches``) must produce the engine's characteristic weight
pattern, and the weights of a complete enumeration must sum to one per
source.

Also covers the driver-level helpers ``_dedup_end_branches``,
``_aggregate_multi_source_mean`` and ``_aggregate_multi_source_fractiles``
with hand-computed expected values.
"""
from __future__ import annotations

import numpy as np
import pytest

from openquake.fdha.logic_tree.driver import (
    _aggregate_multi_source_fractiles,
    _aggregate_multi_source_mean,
    _dedup_end_branches,
)
from openquake.fdha.logic_tree.enumerator import SourceInfo, enumerate_end_branches
from openquake.fdha.logic_tree.types import (
    Branch,
    BranchingLevel,
    BranchSet,
    LogicTreeSpec,
)


def _spec(*branch_sets: BranchSet) -> LogicTreeSpec:
    """One branching level per branch set, matching the canonical layout."""
    levels = tuple(
        BranchingLevel(branching_level_id=f"bl{i}", branch_sets=(bs,))
        for i, bs in enumerate(branch_sets)
    )
    return LogicTreeSpec(logic_tree_id="pfdha", branching_levels=levels, basepath=".")


def _sr_branchset(branches, *, bs_id="bs_sr", style=None, sources=None,
                  apply_to_branches=None) -> BranchSet:
    return BranchSet(
        branch_set_id=bs_id,
        uncertainty_type="fdhaPrimarySRModel",
        apply_to_style=style,
        apply_to_sources=sources,
        apply_to_branches=apply_to_branches,
        branches=tuple(
            Branch(bid, f"[FixedPrimarySR]\nvalue = {v}", str(w))
            for bid, v, w in branches
        ),
    )


def _fd_branchset(branches, *, bs_id="bs_fd", apply_to_branches=None) -> BranchSet:
    return BranchSet(
        branch_set_id=bs_id,
        uncertainty_type="fdhaPrimaryFDModel",
        apply_to_branches=apply_to_branches,
        branches=tuple(
            Branch(bid, "Youngs2003PrimaryFD", str(w)) for bid, w in branches
        ),
    )


SS_SOURCE = SourceInfo(source_id="1", rake=0.0)       # strike-slip
REV_SOURCE = SourceInfo(source_id="2", rake=90.0)     # reverse


def _weight_of(eb):
    return float(eb.weight)


# ------------------------------------------------------- cartesian product


def test_cartesian_product_weights_are_products():
    # 2 SR x 3 FD = 6 end branches; weight = product of the path weights.
    spec = _spec(
        _sr_branchset([("A1", 0.3, 0.4), ("A2", 0.7, 0.6)]),
        _fd_branchset([("B1", 0.25), ("B2", 0.35), ("B3", 0.40)]),
    )
    ebs = enumerate_end_branches(spec, [SS_SOURCE])
    assert len(ebs) == 6

    expected = {
        ("A1", "B1"): 0.4 * 0.25,
        ("A1", "B2"): 0.4 * 0.35,
        ("A1", "B3"): 0.4 * 0.40,
        ("A2", "B1"): 0.6 * 0.25,
        ("A2", "B2"): 0.6 * 0.35,
        ("A2", "B3"): 0.6 * 0.40,
    }
    got = {
        (
            eb.selections["primary_surf_rup"].branch_id,
            eb.selections["primary_surf_displ"].branch_id,
        ): _weight_of(eb)
        for eb in ebs
    }
    assert got.keys() == expected.keys()
    for key in expected:
        assert got[key] == pytest.approx(expected[key], rel=1e-12)
    assert sum(got.values()) == pytest.approx(1.0, rel=1e-12)


def test_partial_extension_matches_engine_weight_pattern():
    # The oq-engine CollapseTestCase tree:
    #    ___/ A-C (w = .4*.5 = .2)
    #  _/   \ A-D (w = .4*.5 = .2)
    #   \____ B   (w = .6, FD level not applied)
    spec = _spec(
        _sr_branchset([("A", 0.5, 0.4), ("B", 1.0, 0.6)]),
        _fd_branchset([("C", 0.5), ("D", 0.5)], apply_to_branches="A"),
    )
    ebs = enumerate_end_branches(spec, [SS_SOURCE])
    assert len(ebs) == 3

    weights = sorted(_weight_of(eb) for eb in ebs)
    assert weights == pytest.approx([0.2, 0.2, 0.6], rel=1e-12)
    assert sum(weights) == pytest.approx(1.0, rel=1e-12)

    # The B path must NOT carry an FD selection; both A paths must.
    for eb in ebs:
        sr_id = eb.selections["primary_surf_rup"].branch_id
        if sr_id == "B":
            assert "primary_surf_displ" not in eb.selections
        else:
            assert eb.selections["primary_surf_displ"].branch_id in {"C", "D"}


def test_apply_to_style_partitions_sources():
    spec = _spec(
        _sr_branchset([("SS1", 0.5, 0.3), ("SS2", 1.0, 0.7)],
                      bs_id="bs_ss", style="strike-slip"),
        _sr_branchset([("RV1", 1.0, 1.0)], bs_id="bs_rev", style="reverse"),
    )
    ebs = enumerate_end_branches(spec, [SS_SOURCE, REV_SOURCE])

    ss = [eb for eb in ebs if eb.source_id == "1"]
    rev = [eb for eb in ebs if eb.source_id == "2"]
    assert len(ss) == 2 and len(rev) == 1
    assert sorted(_weight_of(eb) for eb in ss) == pytest.approx([0.3, 0.7])
    assert _weight_of(rev[0]) == pytest.approx(1.0)
    # Style classification derives from rake: 0 deg -> strike-slip, 90 -> reverse.
    assert {eb.style for eb in ss} == {"strike-slip"}
    assert {eb.style for eb in rev} == {"reverse"}
    # Each source's enumeration is a complete probability space.
    assert sum(_weight_of(eb) for eb in ss) == pytest.approx(1.0)
    assert sum(_weight_of(eb) for eb in rev) == pytest.approx(1.0)


def test_apply_to_sources_restricts_branchset():
    src_a = SourceInfo(source_id="1", rake=0.0)
    src_b = SourceInfo(source_id="2", rake=0.0)
    spec = _spec(
        _sr_branchset([("COM1", 0.5, 0.4), ("COM2", 1.0, 0.6)], bs_id="bs_all"),
        _fd_branchset([("FD1", 0.2), ("FD2", 0.8)], bs_id="bs_only1"),
    )
    # Restrict the FD branch set to source "1" only.
    restricted = BranchSet(
        branch_set_id=spec.branching_levels[1].branch_sets[0].branch_set_id,
        uncertainty_type="fdhaPrimaryFDModel",
        apply_to_sources="1",
        branches=spec.branching_levels[1].branch_sets[0].branches,
    )
    spec = _spec(spec.branching_levels[0].branch_sets[0], restricted)

    ebs = enumerate_end_branches(spec, [src_a, src_b])
    per_1 = [eb for eb in ebs if eb.source_id == "1"]
    per_2 = [eb for eb in ebs if eb.source_id == "2"]
    # Source 1: 2 SR x 2 FD = 4; source 2: SR only = 2.
    assert len(per_1) == 4 and len(per_2) == 2
    assert sum(_weight_of(eb) for eb in per_1) == pytest.approx(1.0)
    assert sum(_weight_of(eb) for eb in per_2) == pytest.approx(1.0)
    assert all("primary_surf_displ" not in eb.selections for eb in per_2)
    w2 = sorted(_weight_of(eb) for eb in per_2)
    assert w2 == pytest.approx([0.4, 0.6])


def test_calc_r_threshold_combines_cartesian():
    # fdhaCalcRThreshold branches multiply into the realisation weights like
    # any model branch but land in the calc pseudo-slot.
    thr_bs = BranchSet(
        branch_set_id="bs_thr",
        uncertainty_type="fdhaCalcRThreshold",
        branches=(Branch("T05", "0.5", "0.3"), Branch("T20", "2.0", "0.7")),
    )
    spec = _spec(_sr_branchset([("A1", 0.5, 0.4), ("A2", 1.0, 0.6)]), thr_bs)
    ebs = enumerate_end_branches(spec, [SS_SOURCE])
    assert len(ebs) == 4

    got = {
        (
            eb.selections["primary_surf_rup"].branch_id,
            eb.selections["calc_r_threshold"].branch_id,
        ): (
            _weight_of(eb),
            eb.selections["calc_r_threshold"].params["r_threshold_km"],
        )
        for eb in ebs
    }
    expected = {
        ("A1", "T05"): (0.4 * 0.3, 0.5),
        ("A1", "T20"): (0.4 * 0.7, 2.0),
        ("A2", "T05"): (0.6 * 0.3, 0.5),
        ("A2", "T20"): (0.6 * 0.7, 2.0),
    }
    assert got.keys() == expected.keys()
    for key, (w, thr) in expected.items():
        assert got[key][0] == pytest.approx(w, rel=1e-12)
        assert got[key][1] == thr
    assert sum(w for w, _ in got.values()) == pytest.approx(1.0, rel=1e-12)


# ------------------------------------------------------------------- dedup


def test_dedup_sums_weights_across_same_style_sources():
    # Two same-style sources see identical selections; the enumerator emits
    # 2 sources x 2 branches = 4 end-branches, dedup collapses to 2 whose
    # weights are summed across sources (total = n_sources).
    spec = _spec(_sr_branchset([("A1", 0.5, 0.4), ("A2", 1.0, 0.6)]))
    src_a = SourceInfo(source_id="1", rake=0.0)
    src_b = SourceInfo(source_id="2", rake=0.0)
    ebs = enumerate_end_branches(spec, [src_a, src_b])
    assert len(ebs) == 4

    deduped = _dedup_end_branches(ebs)
    assert len(deduped) == 2
    by_branch = {
        eb.selections["primary_surf_rup"].branch_id: eb for eb in deduped
    }
    assert by_branch["A1"].weight == pytest.approx(0.8, rel=1e-12)  # 0.4 + 0.4
    assert by_branch["A2"].weight == pytest.approx(1.2, rel=1e-12)  # 0.6 + 0.6
    assert by_branch["A1"].source_id == "1,2"
    assert by_branch["A2"].source_id == "1,2"


def test_dedup_keeps_distinct_selections_apart():
    # Different styles produce different selections; nothing may collapse.
    spec = _spec(
        _sr_branchset([("SS1", 1.0, 1.0)], bs_id="bs_ss", style="strike-slip"),
        _sr_branchset([("RV1", 1.0, 1.0)], bs_id="bs_rev", style="reverse"),
    )
    ebs = enumerate_end_branches(spec, [SS_SOURCE, REV_SOURCE])
    deduped = _dedup_end_branches(ebs)
    assert len(deduped) == 2
    assert {eb.source_id for eb in deduped} == {"1", "2"}


# --------------------------------------------- grouped aggregation helpers


def test_aggregate_multi_source_mean_hand_computed():
    # Group "A": branches 0, 1 with (unnormalised) weights 3, 2.
    # Group "B": branch 2, any positive weight.
    # total = mean_A + mean_B where mean_A uses weights (0.6, 0.4).
    rates = np.array(
        [
            [[10.0, 20.0]],
            [[30.0, 40.0]],
            [[5.0, 7.0]],
        ]
    )  # (3 branches, 1 site, 2 D0)
    weights = [3.0, 2.0, 5.0]
    groups = ["A", "A", "B"]
    total, per_source = _aggregate_multi_source_mean(rates, weights, groups)

    mean_a = np.array([[0.6 * 10 + 0.4 * 30, 0.6 * 20 + 0.4 * 40]])  # [[18, 28]]
    mean_b = rates[2]
    np.testing.assert_allclose(per_source["A"], mean_a, rtol=1e-12)
    np.testing.assert_allclose(per_source["B"], mean_b, rtol=1e-12)
    np.testing.assert_allclose(total, mean_a + mean_b, rtol=1e-12)


def test_aggregate_multi_source_mean_single_group_is_pooled_mean():
    rates = np.array([[[1.0]], [[3.0]]])
    total, per_source = _aggregate_multi_source_mean(
        rates, [0.25, 0.75], ["S", "S"]
    )
    np.testing.assert_allclose(total, [[0.25 * 1 + 0.75 * 3]], rtol=1e-12)
    assert set(per_source) == {"S"}


def test_aggregate_multi_source_fractiles_hand_computed():
    # Group "A": values 1 and 3, equal weights -> sorted xs=[1,3],
    # cdf=[.5,1]. q=0.84 interpolates: 1 + (0.84-0.5)/0.5*(3-1) = 2.36;
    # q=0.5 clamps left (q <= cdf[0]) -> 1.
    # Group "B": single branch -> every fractile equals its value (7).
    rates = np.array([[[1.0]], [[3.0]], [[7.0]]])
    weights = [0.5, 0.5, 1.0]
    groups = ["A", "A", "B"]
    fr = _aggregate_multi_source_fractiles(
        rates, weights, groups, qs=[0.5, 0.84]
    )
    np.testing.assert_allclose(fr[0.5], [[1.0 + 7.0]], rtol=1e-12)
    np.testing.assert_allclose(fr[0.84], [[2.36 + 7.0]], rtol=1e-12)
