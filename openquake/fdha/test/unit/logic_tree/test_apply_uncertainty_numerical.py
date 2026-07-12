"""Numerical checks of NRML uncertainty application on FDHA-parsed sources.

The OpenQuake-engine analog is ``BranchSetApplyUncertaintyTestCase`` in
``openquake/commonlib/tests/logictree_test.py``: apply each uncertainty type
and verify the modified source parameters *numerically*. Here the sources
come through FDHA's own parsing pipeline (``parse_source_model_faults``) and
the uncertainties through :func:`apply_realization_to_sources`, i.e. the
exact plumbing the logic-tree driver uses, so these tests pin down the
contract between the SMLT realisations and the calculator inputs:

* ``abGRAbsolute`` sets (a, b) exactly and the resulting occurrence rates
  match the closed-form truncated Gutenberg-Richter increment
  ``10^(a - b*m_lo) - 10^(a - b*m_hi)`` per magnitude bin;
* ``bGRRelative`` shifts b and conserves the total moment rate (OpenQuake's
  ``increment_b`` rebalances a);
* ``maxMagGRAbsolute`` truncates the magnitude range;
* ``simpleFaultDipAbsolute`` / ``simpleFaultDipRelative`` act on the dip;
* ``applyToSources`` filters exclude non-listed sources untouched;
* the input source dict is never mutated.
"""
from __future__ import annotations

import numpy as np
import pytest

from openquake.fdha.calc.utils.parsing import parse_source_model_faults
from openquake.fdha.logic_tree.source_model_lt import (
    SourceModelBranch,
    SourceModelLogicTreeError,
    UncertaintyApplication,
    apply_realization_to_sources,
)

A0, B0 = 4.2, 0.9
MIN_MAG, MAX_MAG = 6.5, 7.5
BIN_WIDTH = 0.1
DIP0 = 30.0


@pytest.fixture(scope="module")
def base_sources(tmp_path_factory):
    xml = tmp_path_factory.mktemp("srcmodel") / "source_model.xml"
    xml.write_text(
        f'<?xml version="1.0" encoding="utf-8"?>\n'
        f'<nrml xmlns="http://openquake.org/xmlns/nrml/0.5"\n'
        f'      xmlns:gml="http://www.opengis.net/gml">\n'
        f'  <sourceModel name="sm">\n'
        f'    <sourceGroup tectonicRegion="Active Shallow Crust">\n'
        f'      <simpleFaultSource id="3" name="F"\n'
        f'                         tectonicRegion="Active Shallow Crust">\n'
        f'        <simpleFaultGeometry>\n'
        f'          <gml:LineString>\n'
        f'            <gml:posList>16.18 39.69 16.14 39.58</gml:posList>\n'
        f'          </gml:LineString>\n'
        f'          <dip>{DIP0}</dip>\n'
        f'          <upperSeismoDepth>0.0</upperSeismoDepth>\n'
        f'          <lowerSeismoDepth>15.0</lowerSeismoDepth>\n'
        f'        </simpleFaultGeometry>\n'
        f'        <magScaleRel>WC1994</magScaleRel>\n'
        f'        <ruptAspectRatio>2.0</ruptAspectRatio>\n'
        f'        <truncGutenbergRichterMFD aValue="{A0}" bValue="{B0}"\n'
        f'                                  maxMag="{MAX_MAG}" minMag="{MIN_MAG}"/>\n'
        f'        <rake>90.0</rake>\n'
        f'      </simpleFaultSource>\n'
        f'    </sourceGroup>\n'
        f'  </sourceModel>\n'
        f'</nrml>\n'
    )
    return parse_source_model_faults(
        [str(xml)], hdf5path="",
        rupture_mesh_spacing=1.0, width_of_mfd_bin=BIN_WIDTH,
    )


def _branch(*uncertainties: UncertaintyApplication) -> SourceModelBranch:
    return SourceModelBranch(
        branch_id="sm0|" + "|".join(u.branch_id for u in uncertainties),
        source_model_file="source_model.xml",
        weight=1.0,
        uncertainties=tuple(uncertainties),
    )


def _ua(utype, value, *, sources=None) -> UncertaintyApplication:
    return UncertaintyApplication(
        uncertainty_type=utype,
        value=value,
        branch_set_id="bs",
        branch_id="br",
        apply_to_sources=sources,
    )


def _truncated_gr_rates(a, b, min_mag, max_mag, bin_width):
    """Closed-form truncated GR bin rates: 10^(a-b*m1) - 10^(a-b*m2)."""
    edges = np.arange(min_mag, max_mag + bin_width / 2, bin_width)
    return [
        (
            (m1 + m2) / 2,
            10 ** (a - b * m1) - 10 ** (a - b * m2),
        )
        for m1, m2 in zip(edges[:-1], edges[1:])
    ]


def test_ab_gr_absolute_sets_values_and_closed_form_rates(base_sources):
    a1, b1 = 4.8, 1.1
    out = apply_realization_to_sources(
        _branch(_ua("abGRAbsolute", (a1, b1), sources=("3",))), base_sources
    )
    mfd = out["3"].mfd
    assert mfd.a_val == a1
    assert mfd.b_val == b1

    got = mfd.get_annual_occurrence_rates()
    expected = _truncated_gr_rates(a1, b1, MIN_MAG, MAX_MAG, BIN_WIDTH)
    assert len(got) == len(expected) == 10
    for (mag_g, rate_g), (mag_e, rate_e) in zip(got, expected):
        assert mag_g == pytest.approx(mag_e, abs=1e-9)
        assert rate_g == pytest.approx(rate_e, rel=1e-12)


def test_b_gr_relative_shifts_b_and_conserves_moment_rate(base_sources):
    delta = -0.2
    tmr0 = base_sources["3"].mfd._get_total_moment_rate()
    out = apply_realization_to_sources(
        _branch(_ua("bGRRelative", delta)), base_sources
    )
    mfd = out["3"].mfd
    assert mfd.b_val == pytest.approx(B0 + delta, rel=1e-12)
    # OpenQuake's increment_b rebalances a to conserve total moment rate.
    assert mfd.a_val != A0
    assert mfd._get_total_moment_rate() == pytest.approx(tmr0, rel=1e-9)


def test_max_mag_gr_absolute_truncates_bins(base_sources):
    out = apply_realization_to_sources(
        _branch(_ua("maxMagGRAbsolute", 7.0, sources=("3",))), base_sources
    )
    mfd = out["3"].mfd
    assert mfd.max_mag == 7.0
    rates = mfd.get_annual_occurrence_rates()
    assert len(rates) == 5  # (7.0 - 6.5) / 0.1
    assert all(mag < 7.0 for mag, _ in rates)
    # a and b untouched: surviving bins keep their closed-form rates.
    expected = _truncated_gr_rates(A0, B0, MIN_MAG, 7.0, BIN_WIDTH)
    for (mag_g, rate_g), (_, rate_e) in zip(rates, expected):
        assert rate_g == pytest.approx(rate_e, rel=1e-12)


def test_simple_fault_dip_absolute_and_relative(base_sources):
    out = apply_realization_to_sources(
        _branch(_ua("simpleFaultDipAbsolute", 45.0, sources=("3",))),
        base_sources,
    )
    assert out["3"].dip == pytest.approx(45.0)

    out = apply_realization_to_sources(
        _branch(_ua("simpleFaultDipRelative", -10.0)), base_sources
    )
    assert out["3"].dip == pytest.approx(DIP0 - 10.0)


def test_stacked_uncertainties_apply_in_order(base_sources):
    # Same pattern as the engine's test_relative_uncertainty: two stacked
    # modifications on one realisation.
    out = apply_realization_to_sources(
        _branch(
            _ua("maxMagGRRelative", +0.5),
            _ua("bGRRelative", -0.2),
        ),
        base_sources,
    )
    mfd = out["3"].mfd
    assert mfd.max_mag == pytest.approx(MAX_MAG + 0.5, rel=1e-12)
    assert mfd.b_val == pytest.approx(B0 - 0.2, rel=1e-12)


def test_apply_to_sources_filter_excludes_source(base_sources):
    out = apply_realization_to_sources(
        _branch(_ua("abGRAbsolute", (5.0, 1.2), sources=("999",))),
        base_sources,
    )
    # Filter miss: the source is passed through by reference, unmodified.
    assert out["3"] is base_sources["3"]
    assert out["3"].mfd.a_val == A0
    assert out["3"].mfd.b_val == B0


def test_original_sources_never_mutated(base_sources):
    apply_realization_to_sources(
        _branch(
            _ua("abGRAbsolute", (9.9, 2.0), sources=("3",)),
            _ua("simpleFaultDipAbsolute", 89.0, sources=("3",)),
        ),
        base_sources,
    )
    assert base_sources["3"].mfd.a_val == A0
    assert base_sources["3"].mfd.b_val == B0
    assert base_sources["3"].dip == DIP0


def test_unknown_uncertainty_type_raises(base_sources):
    with pytest.raises(SourceModelLogicTreeError, match="Unknown uncertainty"):
        apply_realization_to_sources(
            _branch(_ua("noSuchUncertainty", 1.0)), base_sources
        )


def test_no_uncertainties_returns_equal_dict(base_sources):
    out = apply_realization_to_sources(
        SourceModelBranch(
            branch_id="sm0", source_model_file="x.xml", weight=1.0
        ),
        base_sources,
    )
    assert out == dict(base_sources)
    assert out is not base_sources  # a copy, safe for the caller to edit
