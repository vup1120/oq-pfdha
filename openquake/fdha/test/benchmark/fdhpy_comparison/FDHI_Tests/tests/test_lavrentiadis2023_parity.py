"""
Lavrentiadis & Abrahamson (2023) parity tests.

Three references, none of which needs fdhpy installed:

- golden fixtures frozen from fdhpy 1.0.3 (``fixtures/la23``);
- the authors' original implementation (``reference/la23_pfdha_model``,
  NHR3-UCLA/LA23_PFDHA_model at c3ba1136), run live;
- the paper's equations, where a decision departs from fdhpy: the mixture
  percentile with zero terms (G3), the zero term of the individual-segment
  principal model (B2) and the aggregate-only maximum displacement (G5).
"""

import functools
import importlib.util
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from openquake.fdha.primary_surf_displ import lavrentiadis2023 as la23
from openquake.fdha.primary_surf_displ.lavrentiadis2023 import (
    Lavrentiadis2023PrimaryFD_aggregate as Aggregate,
    Lavrentiadis2023PrimaryFD_principal as Principal,
)

pytestmark = pytest.mark.lavrentiadis2023

HERE = Path(__file__).resolve().parent
FIXTURES = HERE.parent / "fixtures" / "la23"
AUTHORS_CODE = HERE.parent / "reference" / "la23_pfdha_model" / "LA23_flt_rup_model.py"

RTOL, ATOL = 1e-6, 1e-10   # fdhpy fixtures, stored with 10 significant digits
RTOL_MEAN = 1e-8           # power-normal mean (quadrature in both codes)
RTOL_AUTHORS = 1e-12       # same arithmetic, same machine

COMBOS = {
    ("aggregate", "full rupture"): (Aggregate, "disp_agg_prime"),
    ("aggregate", "individual segment"): (Aggregate, "disp_agg_seg"),
    ("sum-of-principal", "full rupture"): (Principal, "disp_prnc_prime"),
    ("sum-of-principal", "individual segment"): (Principal, "disp_prnc_seg"),
}
FULL_RUPTURE = [c for c in COMBOS if c[1] == "full rupture"]
# Table 5, typed from the paper (independent of the module constants)
B0 = {"normal": -3.24, "strike-slip": 0.867, "reverse": 1.65}
B1 = {"normal": 9.105, "strike-slip": 3.767, "reverse": 1.349}


@functools.lru_cache(maxsize=None)
def fixture(name):
    return pd.read_csv(FIXTURES / name, comment="#", float_precision="round_trip")


def displacements():
    return fixture("displacements.csv").displ_m.to_numpy()


def prob_columns(df):
    return [c for c in df.columns if c[0] == "p" and c[1:].isdigit()]


def select(df, metric, version, zero=None):
    mask = (df.metric == metric) & (df.version == version)
    if zero is not None:
        mask &= df.include_prob_zero == zero
    return df[mask]


@functools.lru_cache(maxsize=None)
def authors():
    spec = importlib.util.spec_from_file_location("la23_authors_c3ba1136", AUTHORS_CODE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(autouse=True)
def fresh_warnings(monkeypatch):
    """Warnings are logged once per process; reset for each test."""
    monkeypatch.setattr(la23, "_WARNED", set())


# ----------------------------------------------------------------- F1-F5, F13
PROBABILITY_CASES = [
    (name, method, combo, zero)
    for name, method in [("prob_exceed.csv", "get_prob"), ("cdf.csv", "get_cdf"),
                         ("prob_exceed_out_of_range.csv", "get_prob")]
    for combo in COMBOS
    for zero in (False, True)
    # decision B2 departs from fdhpy here: test_principal_segment_zero_term
    if not (combo == ("sum-of-principal", "individual segment") and zero)
]


@pytest.mark.parametrize("name,method,combo,zero", PROBABILITY_CASES)
def test_probabilities_match_fdhpy(name, method, combo, zero):
    cls, output_type = COMBOS[combo]
    df = fixture(name)
    cols = prob_columns(df)
    model = cls()
    d = displacements()
    for row in select(df, *combo, zero).itertuples(index=False):
        got = getattr(model, method)(
            d, np.array([row.xl]), row.magnitude, style=row.style,
            output_type=output_type, include_zero_slip=zero)
        np.testing.assert_allclose(
            got, [getattr(row, c) for c in cols], rtol=RTOL, atol=ATOL,
            err_msg=f"{name} {combo} zero={zero}: {row.style} M{row.magnitude} x/L {row.xl}")


def test_aggregate_segment_ignores_zero_terms(caplog):
    """F3: a single segment has no gap (p. 21) and P(D_P = 0) is the
    zero-PRINCIPAL probability, so the flag is ignored with a warning.
    Non-vacuous: P(D_P = 0) exceeds 0.5 at the probe."""
    caplog.set_level(logging.WARNING, logger=la23.__name__)
    model, d, xl = Aggregate(), displacements(), np.array([0.05])
    assert model.get_p_zero_slip(xl, 5.5, style="normal")[0] > 0.5
    with_zero = model.get_prob(d, xl, 5.5, style="normal",
                               output_type="disp_agg_seg", include_zero_slip=True)
    without = model.get_prob(d, xl, 5.5, style="normal", output_type="disp_agg_seg")
    np.testing.assert_array_equal(with_zero, without)
    assert "include_zero_slip is ignored for disp_agg_seg" in caplog.text


@pytest.mark.parametrize("mag", [4.99, 8.51])
def test_out_of_range_magnitude_warns_and_extrapolates(mag, caplog):
    """F13 (decision Q3): warn, then compute as fdhpy does; the values are
    checked against prob_exceed_out_of_range.csv above."""
    caplog.set_level(logging.WARNING, logger=la23.__name__)
    for _ in range(3):
        out = Principal().get_prob(displacements(), np.array([0.3]), mag,
                                   style="strike-slip", include_zero_slip=True)
    assert np.all(np.isfinite(out))
    records = [r for r in caplog.records if "outside the model range" in r.message]
    assert len(records) == 1


@pytest.mark.parametrize("mag", [5.0, 8.5])
def test_range_edges_do_not_warn(mag, caplog):
    caplog.set_level(logging.WARNING, logger=la23.__name__)
    Aggregate().get_prob(displacements(), np.array([0.3]), mag, style="reverse")
    assert "outside the model range" not in caplog.text


@pytest.mark.parametrize("style", ["all", "strike_slip", "thrust"])
def test_invalid_style_raises(style):
    with pytest.raises(ValueError, match="style"):
        Aggregate(style=style)
    with pytest.raises(ValueError, match="style"):
        Aggregate().get_prob(np.array([0.1]), np.array([0.3]), 7.0, style=style)


def test_missing_style_warns_and_uses_normal(caplog):
    caplog.set_level(logging.WARNING, logger=la23.__name__)
    d, xl = displacements(), np.array([0.3])
    np.testing.assert_array_equal(
        Aggregate().get_prob(d, xl, 7.0),
        Aggregate().get_prob(d, xl, 7.0, style="normal"))
    assert "no style given" in caplog.text


def test_x_over_l_is_folded():
    d, model = displacements(), Principal()
    for a, b in [(0.77, 0.23), (1.2, 0.2), (1.0, 0.0)]:
        np.testing.assert_allclose(
            model.get_prob(d, np.array([a]), 6.5, style="reverse", include_zero_slip=True),
            model.get_prob(d, np.array([b]), 6.5, style="reverse", include_zero_slip=True),
            rtol=1e-13, atol=0)


# ------------------------------------------------------------------- F4, B2
def test_principal_segment_zero_term():
    """B2: P(D_P > d) = (1 - P(D_P = 0)) * fdhpy's segment exceedance, with
    P(D_P = 0) from Eq. 32 and the Eq. 14 predictor (fdhpy's segment mu) and
    P(Gap) = 0 for a single segment (p. 21)."""
    pe, params = fixture("prob_exceed.csv"), fixture("params.csv")
    cols, d = prob_columns(pe), displacements()
    mu_seg = select(params, "aggregate", "individual segment").set_index(
        ["style", "magnitude", "xl"]).mu
    largest = 0.0
    for row in select(pe, "sum-of-principal", "individual segment", False).itertuples(index=False):
        mu = mu_seg[(row.style, row.magnitude, row.xl)]
        p0 = 1.0 / (1.0 + np.exp(B0[row.style] + B1[row.style] * mu))
        largest = max(largest, p0)
        got = Principal().get_prob(d, np.array([row.xl]), row.magnitude, style=row.style,
                                   output_type="disp_prnc_seg", include_zero_slip=True)
        np.testing.assert_allclose(
            got, (1.0 - p0) * np.array([getattr(row, c) for c in cols]),
            rtol=RTOL, atol=ATOL, err_msg=f"{row.style} M{row.magnitude} x/L {row.xl}")
    assert largest > 0.5  # the term is not vacuous on the grid


@pytest.mark.parametrize("combo", list(COMBOS))
def test_stat_params_match_fdhpy(combo):
    """F16, including the F4 median and sigma of disp_prnc_seg."""
    cls, output_type = COMBOS[combo]
    for row in select(fixture("params.csv"), *combo).itertuples(index=False):
        got = cls().get_stat_params(np.array([row.xl]), row.magnitude, style=row.style,
                                    output_type=output_type)
        np.testing.assert_allclose([got["mu"][0], got["sigma"][0]], [row.mu, row.sigma],
                                   rtol=RTOL, atol=ATOL)


# ---------------------------------------------------------------- F11, F12
def test_zero_probabilities_and_epistemic_sigma_match_fdhpy():
    rows = select(fixture("params.csv"), "sum-of-principal", "full rupture")
    model = Aggregate()
    for row in rows.itertuples(index=False):
        xl = np.array([row.xl])
        np.testing.assert_allclose(
            [model.get_p_gap(xl, row.magnitude, style=row.style)[0],
             model.get_p_zero_slip(xl, row.magnitude, style=row.style)[0],
             model.get_sigma_mu_agg(row.magnitude, style=row.style)],
            [row.p_gap, row.p_zero_slip, row.sigma_mu_agg], rtol=RTOL, atol=ATOL)


@pytest.mark.parametrize("mag,style,expected", [
    (6.0, "normal", 0.035 + 0.064 * 1.1),
    (6.0, "reverse", 0.035 + 0.036 * 1.1),
    (8.0, "strike-slip", 0.035 + 0.025 * 0.9),
])
def test_sigma_mu_agg_follows_eq29(mag, style, expected):
    assert Aggregate().get_sigma_mu_agg(mag, style=style) == pytest.approx(expected, rel=1e-12)


# ------------------------------------------------------------------- F6-F8
def _site_cases(df, combo, zero, mean_only=False):
    rows = select(df, *combo, zero)
    if mean_only:
        rows = rows[rows.percentile == -1]
    return rows.groupby(["style", "magnitude", "percentile"], sort=False)


@pytest.mark.parametrize("combo", list(COMBOS))
def test_displ_site_without_zero_terms_matches_fdhpy(combo):
    cls, output_type = COMBOS[combo]
    for (style, mag, prc), g in _site_cases(fixture("displ_site.csv"), combo, False):
        got = cls().get_displ_site(g.xl.to_numpy(), mag, percentile=prc, style=style,
                                   output_type=output_type)
        np.testing.assert_allclose(
            got, g.displ_m.to_numpy(), rtol=RTOL_MEAN if prc == -1 else RTOL, atol=1e-12,
            err_msg=f"{combo} {style} M{mag} percentile {prc}")


@pytest.mark.parametrize("combo", FULL_RUPTURE)
def test_mean_with_zero_terms_matches_fdhpy(combo):
    """G3: with the zero terms the mean is (1 - P0) times the mean of the
    non-zero part, as in fdhpy (and Eq. 35 scaled by 1 - P(Gap), p. 24)."""
    cls, output_type = COMBOS[combo]
    for (style, mag, prc), g in _site_cases(fixture("displ_site.csv"), combo, True, True):
        got = cls().get_displ_site(g.xl.to_numpy(), mag, percentile=-1, style=style,
                                   output_type=output_type, include_zero_slip=True)
        np.testing.assert_allclose(got, g.displ_m.to_numpy(), rtol=RTOL_MEAN, atol=1e-12)


@pytest.mark.parametrize("combo", list(COMBOS))
def test_percentile_with_zero_terms_is_the_mixture_percentile(combo):
    """G3: d_p is the root of P(D > d) = 1 - p, the zero-scaled exceedance
    of Eqs. 25 and 31 (B2 for the principal segment); d_p = 0 where the
    zero mass alone exceeds p."""
    cls, output_type = COMBOS[combo]
    model, n_zero = cls(), 0
    for style in la23.STYLES:
        for mag in (5.0, 5.5, 6.5, 7.5, 8.5):
            xl = np.array([0.01, 0.05, 0.15, 0.3, 0.5])
            for p in (0.16, 0.5, 0.84):
                d_p = model.get_displ_site(xl, mag, percentile=p, style=style,
                                           output_type=output_type, include_zero_slip=True)
                exceed = np.array([
                    model.get_prob(np.array([max(d, 1e-300)]), np.array([x]), mag, style=style,
                                   output_type=output_type, include_zero_slip=True)
                    for d, x in zip(d_p, xl)]).ravel()
                pos = d_p > 0
                np.testing.assert_allclose(exceed[pos], 1 - p, rtol=1e-9)
                assert np.all(exceed[~pos] <= 1 - p + 1e-12)
                n_zero += int((~pos).sum())
    if combo[0] == "sum-of-principal":
        assert n_zero > 0  # the zero mass decides some percentiles


def test_mixture_percentile_differs_from_fdhpy_as_documented():
    """fdhpy multiplies the percentile displacement by 1 - P0 (frozen in the
    fixture); pfdha returns the mixture percentile. LA23 M6.5 strike-slip
    x/L 0.23, principal median: 0.1258 m (pfdha) vs 0.1349 m (fdhpy)."""
    df = select(fixture("displ_site.csv"), "sum-of-principal", "full rupture")
    key = (df["style"] == "strike-slip") & (df.magnitude == 6.5) & (df.xl == 0.23) & (df.percentile == 0.5)
    fdhpy_zero = df[key & df.include_prob_zero].displ_m.item()
    fdhpy_plain = df[key & ~df.include_prob_zero].displ_m.item()
    model = Principal()
    w = (1 - model.get_p_gap(np.array([0.23]), 6.5, style="strike-slip")[0]) * (
        1 - model.get_p_zero_slip(np.array([0.23]), 6.5, style="strike-slip")[0])
    assert fdhpy_zero == pytest.approx(w * fdhpy_plain, rel=1e-6)
    ours = model.get_displ_site(np.array([0.23]), 6.5, percentile=0.5, style="strike-slip",
                                include_zero_slip=True)[0]
    assert ours == pytest.approx(0.1258, abs=5e-5)
    assert fdhpy_zero == pytest.approx(0.1349, abs=5e-5)


@pytest.mark.parametrize("combo", list(COMBOS))
def test_displ_profile_matches_fdhpy(combo):
    cls, output_type = COMBOS[combo]
    df = fixture("displ_profile.csv")
    for zero in (False, True):
        rows = select(df, *combo, zero)
        if zero:  # targets with zero terms: the full-rupture mean only (G3, B2)
            if combo[1] != "full rupture":
                continue
            rows = rows[rows.percentile == -1]
        for (style, mag, prc), g in rows.groupby(["style", "magnitude", "percentile"]):
            xl, got = cls().get_displ_profile(mag, percentile=prc, xl_step=0.1, style=style,
                                              output_type=output_type, include_zero_slip=zero)
            np.testing.assert_allclose(xl, g.xl.to_numpy(), rtol=0, atol=1e-12)
            np.testing.assert_allclose(got, g.displ_m.to_numpy(),
                                       rtol=RTOL_MEAN if prc == -1 else RTOL, atol=1e-12)


def test_profile_grid_closes_at_one():
    """fdhpy returns no profile when xl_step does not divide 1."""
    xl, displ = Aggregate().get_displ_profile(7.0, xl_step=0.03, style="normal")
    assert xl[0] == 0.0 and xl[-1] == 1.0 and np.all(np.diff(xl) <= 0.03 + 1e-12)
    assert displ.shape == xl.shape


@pytest.mark.parametrize("percentile", [0.0, 1.0, 1.5, -2])
def test_invalid_percentile_raises(percentile):
    with pytest.raises(ValueError, match="percentile"):
        Aggregate().get_displ_site(np.array([0.3]), 7.0, percentile=percentile, style="normal")


# ------------------------------------------------------------------ F9, F10
def test_average_displacement_matches_fdhpy():
    for row in fixture("displ_avg.csv").itertuples(index=False):
        assert Principal().get_displ_avg(row.magnitude, style=row.style) == pytest.approx(
            row.displ_m, rel=RTOL)


def test_average_displacement_is_principal_full_rupture_only():
    with pytest.raises(ValueError, match="principal displacement only"):
        Aggregate().get_displ_avg(7.0, style="normal")
    with pytest.raises(ValueError, match="full rupture only"):
        Principal(output_type="disp_prnc_seg").get_displ_avg(7.0, style="normal")


def test_maximum_displacement_matches_fdhpy():
    for row in fixture("displ_max.csv").itertuples(index=False):
        got = Aggregate().get_displ_max(row.magnitude, percentile=row.percentile, style=row.style)
        assert got == pytest.approx(row.displ_m, rel=RTOL_MEAN if row.percentile == -1 else RTOL)


def test_maximum_displacement_is_aggregate_full_rupture_only():
    """G5: Eqs. 38-39 and Table 6 give the maximum AGGREGATE displacement;
    fdhpy's sum-of-principal variant is not in the paper."""
    with pytest.raises(ValueError, match="aggregate displacement only"):
        Principal().get_displ_max(7.0, style="normal")
    with pytest.raises(ValueError, match="full rupture only"):
        Aggregate(output_type="disp_agg_seg").get_displ_max(7.0, style="normal")


# ------------------------------------------------- authors' implementation
MAGS = (5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5)
XLS = np.array([0.01, 0.05, 0.10, 0.15, 0.23, 0.30, 0.40, 0.50, 0.77])
SRL = 30_000.0


@pytest.mark.parametrize("sof", ["Normal", "Strike-Slip", "Reverse"])
def test_slip_profile_matches_authors_code(sof):
    """The 11 terms of LavrentiadisAbrahamson2023SlipProfile (medians,
    sigmas, P(Gap), P(D_P = 0))."""
    for mag in MAGS:
        ref = authors().LavrentiadisAbrahamson2023SlipProfile(XLS * SRL, mag, SRL, sof)
        got = Aggregate().LavrentiadisAbrahamson2023SlipProfile(XLS * SRL, mag, SRL, sof)
        for i, (g, r) in enumerate(zip(got, ref)):
            r = np.ravel(r)
            np.testing.assert_allclose(np.broadcast_to(np.ravel(g), r.shape), r,
                                       rtol=RTOL_AUTHORS, atol=0, err_msg=f"term {i} M{mag}")


@pytest.mark.parametrize("sof", ["Normal", "Strike-Slip", "Reverse"])
def test_percentiles_match_authors_code(sof):
    """Full-rupture aggregate and principal percentiles (no zero terms)."""
    style = sof.lower()
    for mag in MAGS:
        for prc in (0.16, 0.5, 0.84):
            ref_agg, ref_prnc = authors().LavrentiadisAbrahamson2023SlipProfilePrc(
                XLS * SRL, mag, SRL, sof, prc)
            agg = Aggregate().get_displ_site(XLS, mag, percentile=prc, style=style)
            prnc = Principal().get_displ_site(XLS, mag, percentile=prc, style=style)
            np.testing.assert_allclose(agg, ref_agg, rtol=RTOL_AUTHORS, atol=0)
            np.testing.assert_allclose(prnc, ref_prnc, rtol=RTOL_AUTHORS, atol=0)
            legacy = Aggregate().LavrentiadisAbrahamson2023SlipProfilePrc(
                XLS * SRL, mag, SRL, sof, prc)
            np.testing.assert_allclose(legacy[0], ref_agg, rtol=RTOL_AUTHORS, atol=0)
            np.testing.assert_allclose(legacy[1], ref_prnc, rtol=RTOL_AUTHORS, atol=0)


@pytest.mark.parametrize("sof", ["Normal", "Strike-Slip", "Reverse"])
def test_average_and_maximum_displacement_match_authors_code(sof):
    style = sof.lower()
    for mag in MAGS:
        ad, ratio = authors().LavrentiadisAbrahamson2023AvgDisp(mag, SRL, sof)
        md, sig_md = authors().LavrentiadisAbrahamson2023MaxDisp(mag, SRL, sof)
        assert Principal().get_displ_avg(mag, style=style) == pytest.approx(
            np.ravel(ad)[0], rel=RTOL_AUTHORS)
        assert Aggregate().get_displ_max(mag, style=style) == pytest.approx(
            np.ravel(md)[0], rel=RTOL_AUTHORS)
        legacy_ad = Aggregate().LavrentiadisAbrahamson2023AvgDisp(mag, SRL, sof)
        legacy_md = Aggregate().LavrentiadisAbrahamson2023MaxDisp(mag, SRL, sof)
        np.testing.assert_allclose(legacy_ad, [np.ravel(ad)[0], np.ravel(ratio)[0]],
                                   rtol=RTOL_AUTHORS)
        np.testing.assert_allclose(legacy_md, [np.ravel(md)[0], np.ravel(sig_md)[0]],
                                   rtol=RTOL_AUTHORS)
