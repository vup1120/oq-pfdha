# -*- coding: utf-8 -*-
"""
Tests for the ``version`` parameter of Youngs2003PrimarySR.

The parameter used to be called ``style``, which was misleading: it never
selected a faulting style, only which published data set Equation 4's
coefficients were fitted to. Two of the four data sets tabulated in the
Youngs et al. (2003) Appendix ("Coefficients for Equation 4 shown on
Figure 4") were also missing.

These tests pin:

1. all four Appendix coefficient pairs;
2. that ``version="WC93"`` is the Wells & Coppersmith (1993) regression,
   i.e. identical to WC1993PrimarySR - the two must never be treated as
   independent models on one logic-tree branch set;
3. the ``style`` -> ``version`` backward-compatible aliases, so inputs
   written before the rename reproduce bit-for-bit;
4. the selection precedence. This one matters in practice: the adapter
   injects a ``style`` derived from the rupture rake on *every* call
   (``calc/model_adapter.py:_resolve_style``), so a branch pinning
   ``version`` must not be overridden by a style it never asked for.
"""

import numpy as np
import pytest

from openquake.fdha.primary_surf_rup.youngs2003 import Youngs2003PrimarySR
from openquake.fdha.primary_surf_rup.wells_coppersmith1993 import (
    WC1993PrimarySR)


# Youngs et al. (2003), Appendix, "Coefficients for Equation 4 shown on
# Figure 4". WC93 is the Wells & Coppersmith (1993) worldwide all-slip-type
# regression; the other three are the Pezzopane & Dawson (1996)
# normal-faulting data sets.
APPENDIX_COEFFS = {
    "WC93": (-12.51, 2.053),
    "GreatBasin": (-16.02, 2.685),
    "NorthernBasinAndRange": (-18.71, 3.041),
    "ExtensionalCordillera": (-12.53, 1.921),
}

MAGS = np.linspace(4.0, 8.5, 91)


def logistic(a, b, mag):
    fx = a + b * np.asarray(mag, dtype=float)
    return np.exp(fx) / (1.0 + np.exp(fx))


@pytest.mark.parametrize("version,coeffs", sorted(APPENDIX_COEFFS.items()))
def test_version_reproduces_appendix_coefficients(version, coeffs):
    got = Youngs2003PrimarySR(version=version).get_prob(MAGS)
    np.testing.assert_allclose(got, logistic(*coeffs, MAGS), rtol=1e-12)


def test_wc93_version_is_the_wells_coppersmith_model():
    """The two are the same regression, not independent models."""
    np.testing.assert_array_equal(
        Youngs2003PrimarySR(version="WC93").get_prob(MAGS),
        WC1993PrimarySR().get_prob(MAGS))


@pytest.mark.parametrize("spelling", [
    "GreatBasin", "greatbasin", "Great Basin", "GREAT_BASIN", "great-basin",
    "gb",
])
def test_version_spelling_is_forgiving(spelling):
    np.testing.assert_allclose(
        Youngs2003PrimarySR(version=spelling).get_prob(MAGS),
        logistic(*APPENDIX_COEFFS["GreatBasin"], MAGS), rtol=1e-12)


@pytest.mark.parametrize("style,version", [
    ("all", "WC93"),
    ("normal", "GreatBasin"),
])
def test_legacy_style_still_selects_the_same_dataset(style, version):
    """Inputs written before the rename must be unchanged."""
    np.testing.assert_allclose(
        Youngs2003PrimarySR(style=style).get_prob(MAGS),
        logistic(*APPENDIX_COEFFS[version], MAGS), rtol=1e-12)


def test_legacy_default_is_wc93():
    """No selector anywhere: keep the historical default."""
    np.testing.assert_allclose(
        Youngs2003PrimarySR().get_prob(MAGS),
        logistic(*APPENDIX_COEFFS["WC93"], MAGS), rtol=1e-12)


def test_call_time_style_selects_dataset_when_nothing_is_pinned():
    """A bare branch takes the rake-derived style the adapter injects."""
    np.testing.assert_allclose(
        Youngs2003PrimarySR().get_prob(MAGS, style="normal"),
        logistic(*APPENDIX_COEFFS["GreatBasin"], MAGS), rtol=1e-12)


def test_pinned_version_outranks_injected_style():
    """
    The adapter injects a rake-derived style on every call. A branch that
    pins ``version`` must keep it, otherwise a normal-rake source would
    silently serve Great Basin coefficients to a WC93 branch.
    """
    np.testing.assert_allclose(
        Youngs2003PrimarySR(version="WC93").get_prob(MAGS, style="normal"),
        logistic(*APPENDIX_COEFFS["WC93"], MAGS), rtol=1e-12)


def test_call_time_version_outranks_everything():
    np.testing.assert_allclose(
        Youngs2003PrimarySR(version="WC93").get_prob(
            MAGS, version="ExtensionalCordillera", style="normal"),
        logistic(*APPENDIX_COEFFS["ExtensionalCordillera"], MAGS),
        rtol=1e-12)


@pytest.mark.parametrize("kwargs", [
    {"version": "NotADataset"},
    {"style": "reverse"},
])
def test_unknown_selector_is_rejected(kwargs):
    with pytest.raises(ValueError, match="Youngs2003PrimarySR"):
        Youngs2003PrimarySR(**kwargs)


def test_scalar_input_returns_scalar():
    assert isinstance(
        Youngs2003PrimarySR(version="GreatBasin").get_prob(6.5), float)


def test_normal_faulting_sets_exceed_worldwide_above_magnitude_6():
    """
    Sanity check on Figure 4: the Great Basin fit is steeper than the
    worldwide all-slip-type one, so it predicts more surface rupture at
    moderate-to-large magnitude.
    """
    mags = np.array([6.0, 6.5, 7.0])
    gb = Youngs2003PrimarySR(version="GreatBasin").get_prob(mags)
    wc = Youngs2003PrimarySR(version="WC93").get_prob(mags)
    assert np.all(gb > wc)
