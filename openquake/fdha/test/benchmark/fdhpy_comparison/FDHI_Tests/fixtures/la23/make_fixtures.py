# -*- coding: utf-8 -*-
"""
Freeze the fdhpy ``LavrentiadisAbrahamson2023`` outputs as golden fixtures.

The CSV files next to this script are the regression baseline for the LA23
parity work (features F1-F12 and F16). They are generated once from the
pinned fdhpy release and read by the tests without fdhpy installed.

Run with the pinned reference installed (``pip install -e .[test]``)::

    python make_fixtures.py

The script refuses to run against any other fdhpy version: fdhpy <= 1.0.2
returns a wrong CDF when the zero-displacement terms are included.
"""

import datetime
import importlib.metadata
import logging
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

FDHPY_VERSION = "1.0.3"

HERE = Path(__file__).resolve().parent
MAGNITUDES = (5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5)
XLS = (0.01, 0.05, 0.10, 0.15, 0.23, 0.30, 0.40, 0.50, 0.77)
STYLES = ("normal", "strike-slip", "reverse")
OUT_OF_RANGE_MAGNITUDES = (4.99, 8.51)
OUT_OF_RANGE_XLS = (0.05, 0.30, 0.50)
PERCENTILES = (0.16, 0.5, 0.84, -1)  # -1 = mean
COMBOS = (
    ("aggregate", "full rupture"),
    ("aggregate", "individual segment"),
    ("sum-of-principal", "full rupture"),
    ("sum-of-principal", "individual segment"),
)
INCLUDE_PROB_ZERO = (False, True)
DISPLACEMENTS = np.logspace(-3, np.log10(20.0), 60)
PROFILE_XL_STEP = 0.1  # must divide 1: fdhpy returns None otherwise
FLOAT_FORMAT = "%.10g"
PROB_COLUMNS = [f"p{i:02d}" for i in range(DISPLACEMENTS.size)]


def _check_version():
    installed = importlib.metadata.version("fdhpy")
    if installed != FDHPY_VERSION:
        raise SystemExit(
            f"fdhpy {installed} is installed; the fixtures must be "
            f"generated with fdhpy=={FDHPY_VERSION}")


def _header(description, columns):
    lines = [
        f"generator: {Path(__file__).name}",
        f"fdhpy: {FDHPY_VERSION} (LavrentiadisAbrahamson2023)",
        f"numpy: {np.__version__}; scipy: {scipy.__version__}; "
        f"pandas: {pd.__version__}; python: {platform.python_version()}",
        f"generated: {datetime.date.today().isoformat()}",
        f"content: {description}",
    ]
    lines += [f"column {name}: {meaning}" for name, meaning in columns]
    return "".join(f"# {line}\n" for line in lines)


def _write(name, df, description, columns, float_format=FLOAT_FORMAT):
    path = HERE / name
    with path.open("w", newline="") as f:
        f.write(_header(description, columns))
        df.to_csv(f, index=False, float_format=float_format)
    print(f"{name}: {len(df)} rows, {path.stat().st_size / 1e3:.0f} kB")


def _nan_if_none(value):
    return np.nan if value is None else float(value)


SCENARIO_COLUMNS = [
    ("style", "faulting style (fdhpy spelling)"),
    ("magnitude", "moment magnitude"),
    ("xl", "x/L passed to fdhpy (unfolded)"),
    ("metric", "fdhpy metric"),
    ("version", "fdhpy version"),
    ("include_prob_zero", "fdhpy include_prob_zero flag"),
]


def make_displacements():
    df = pd.DataFrame({"column": PROB_COLUMNS, "displ_m": DISPLACEMENTS})
    _write("displacements.csv", df,
           "the 60 log-spaced test displacements (1e-3 to 20 m) used by "
           "prob_exceed.csv and cdf.csv",
           [("column", "column name in prob_exceed.csv / cdf.csv"),
            ("displ_m", "displacement in metres")],
           float_format="%.17g")  # exact round trip: the grid is an input


def make_prob_exceed_out_of_range(model):
    # fdhpy warns outside its recommended magnitude range and still computes
    rows = []
    for style in STYLES:
        for mag in OUT_OF_RANGE_MAGNITUDES:
            for xl in OUT_OF_RANGE_XLS:
                for metric, version in COMBOS:
                    for zero in INCLUDE_PROB_ZERO:
                        m = model(style=style, magnitude=mag, xl=xl,
                                  metric=metric, version=version,
                                  displ_array=DISPLACEMENTS,
                                  include_prob_zero=zero)
                        rows.append([style, mag, xl, metric, version, zero]
                                    + list(m.prob_exceed))
    names = [c for c, _ in SCENARIO_COLUMNS] + PROB_COLUMNS
    _write("prob_exceed_out_of_range.csv", pd.DataFrame(rows, columns=names),
           "P(D > d) just outside the recommended range M 5.0-8.5 (F13)",
           SCENARIO_COLUMNS + [
               ("p00..p59",
                "value at the displacement listed in displacements.csv")])


def make_prob_exceed_and_cdf(model):
    pe_rows, cdf_rows = [], []
    for style in STYLES:
        for mag in MAGNITUDES:
            for xl in XLS:
                for metric, version in COMBOS:
                    for zero in INCLUDE_PROB_ZERO:
                        m = model(style=style, magnitude=mag, xl=xl,
                                  metric=metric, version=version,
                                  displ_array=DISPLACEMENTS,
                                  include_prob_zero=zero)
                        key = [style, mag, xl, metric, version, zero]
                        pe_rows.append(key + list(m.prob_exceed))
                        cdf_rows.append(key + list(m.cdf))
    names = [c for c, _ in SCENARIO_COLUMNS] + PROB_COLUMNS
    columns = SCENARIO_COLUMNS + [
        ("p00..p59", "value at the displacement listed in displacements.csv")]
    _write("prob_exceed.csv", pd.DataFrame(pe_rows, columns=names),
           "P(D > d), fdhpy prob_exceed (F1-F4)", columns)
    _write("cdf.csv", pd.DataFrame(cdf_rows, columns=names),
           "P(D <= d), fdhpy cdf (F5)", columns)


def make_displ_site(model):
    rows = []
    for style in STYLES:
        for mag in MAGNITUDES:
            for xl in XLS:
                for metric, version in COMBOS:
                    for zero in INCLUDE_PROB_ZERO:
                        for prc in PERCENTILES:
                            m = model(style=style, magnitude=mag, xl=xl,
                                      metric=metric, version=version,
                                      percentile=prc, include_prob_zero=zero)
                            rows.append([style, mag, xl, metric, version, zero,
                                         prc, _nan_if_none(m.displ_site)])
    names = [c for c, _ in SCENARIO_COLUMNS] + ["percentile", "displ_m"]
    _write("displ_site.csv", pd.DataFrame(rows, columns=names),
           "fdhpy displ_site (F6, F7); NaN where fdhpy returns None "
           "(individual segment with include_prob_zero=True)",
           SCENARIO_COLUMNS + [("percentile", "aleatory quantile; -1 = mean"),
                               ("displ_m", "displacement in metres")])


def make_displ_profile(model):
    rows = []
    for style in STYLES:
        for mag in MAGNITUDES:
            for metric, version in COMBOS:
                for zero in INCLUDE_PROB_ZERO:
                    for prc in PERCENTILES:
                        m = model(style=style, magnitude=mag, metric=metric,
                                  version=version, percentile=prc,
                                  xl_step=PROFILE_XL_STEP,
                                  include_prob_zero=zero)
                        xl_array, displ = m.displ_profile
                        for xl, d in zip(xl_array, displ):
                            rows.append([style, mag, metric, version, zero,
                                         prc, float(xl), _nan_if_none(d)])
    names = ["style", "magnitude", "metric", "version", "include_prob_zero",
             "percentile", "xl", "displ_m"]
    _write("displ_profile.csv", pd.DataFrame(rows, columns=names),
           f"fdhpy displ_profile with xl_step={PROFILE_XL_STEP} (F8)",
           [c for c in SCENARIO_COLUMNS if c[0] != "xl"]
           + [("percentile", "aleatory quantile; -1 = mean"),
              ("xl", "profile x/L"), ("displ_m", "displacement in metres")])


def make_displ_avg(model):
    rows = []
    for style in STYLES:
        for mag in MAGNITUDES:
            m = model(style=style, magnitude=mag, metric="sum-of-principal",
                      version="full rupture", percentile=0.5,
                      include_prob_zero=False)
            rows.append([style, mag, float(m.displ_avg)])
    _write("displ_avg.csv",
           pd.DataFrame(rows, columns=["style", "magnitude", "displ_m"]),
           "fdhpy displ_avg, sum-of-principal, full rupture, median (F9)",
           [("style", "faulting style"), ("magnitude", "moment magnitude"),
            ("displ_m", "average displacement in metres")])


def make_displ_max(model):
    rows = []
    for style in STYLES:
        for mag in MAGNITUDES:
            for prc in PERCENTILES:
                m = model(style=style, magnitude=mag, metric="aggregate",
                          version="full rupture", percentile=prc,
                          include_prob_zero=False)
                rows.append([style, mag, prc, float(m.displ_max)])
    _write("displ_max.csv",
           pd.DataFrame(rows, columns=["style", "magnitude", "percentile",
                                       "displ_m"]),
           "fdhpy displ_max, aggregate, full rupture (F10); the "
           "sum-of-principal variant is deliberately not frozen",
           [("style", "faulting style"), ("magnitude", "moment magnitude"),
            ("percentile", "aleatory quantile; -1 = mean"),
            ("displ_m", "maximum displacement in metres")])


def make_params(model):
    rows = []
    for style in STYLES:
        for mag in MAGNITUDES:
            for xl in XLS:
                for metric, version in COMBOS:
                    m = model(style=style, magnitude=mag, xl=xl, metric=metric,
                              version=version, include_prob_zero=True)
                    params = m.stat_params_info["params"]
                    rows.append([style, mag, xl, metric, version,
                                 params["mu"], params["sigma"], m.p_gap,
                                 m.p_zero_slip, m.sigma_mu_agg])
    names = ["style", "magnitude", "xl", "metric", "version", "mu", "sigma",
             "p_gap", "p_zero_slip", "sigma_mu_agg"]
    _write("params.csv", pd.DataFrame(rows, columns=names),
           "power-normal parameters, zero-displacement probabilities and "
           "epistemic sigma, include_prob_zero=True (F11, F12, F16)",
           [c for c in SCENARIO_COLUMNS if c[0] != "include_prob_zero"]
           + [("mu", "median of D^0.3"), ("sigma", "aleatory sigma of D^0.3"),
              ("p_gap", "fdhpy p_gap (0 for individual segment)"),
              ("p_zero_slip", "fdhpy p_zero_slip"),
              ("sigma_mu_agg", "fdhpy sigma_mu_agg (Eq. 29)")])


def main():
    _check_version()
    from fdhpy import LavrentiadisAbrahamson2023

    # fdhpy logs a warning for every flag it ignores; those cases are
    # documented in README.md
    logging.disable(logging.WARNING)
    make_displacements()
    make_prob_exceed_and_cdf(LavrentiadisAbrahamson2023)
    make_prob_exceed_out_of_range(LavrentiadisAbrahamson2023)
    make_displ_site(LavrentiadisAbrahamson2023)
    make_displ_profile(LavrentiadisAbrahamson2023)
    make_displ_avg(LavrentiadisAbrahamson2023)
    make_displ_max(LavrentiadisAbrahamson2023)
    make_params(LavrentiadisAbrahamson2023)


if __name__ == "__main__":
    main()
