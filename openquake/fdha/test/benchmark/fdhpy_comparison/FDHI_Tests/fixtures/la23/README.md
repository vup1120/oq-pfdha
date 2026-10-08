# LA23 golden fixtures (fdhpy 1.0.3)

Frozen outputs of fdhpy 1.0.3 `LavrentiadisAbrahamson2023`. They are the regression baseline for the
Lavrentiadis & Abrahamson (2023) parity tests, and they let those tests run without fdhpy installed.

Regenerate with `python make_fixtures.py`. The script refuses any fdhpy version other than 1.0.3:
fdhpy ≤ 1.0.2 returns a wrong `cdf` when the zero-displacement terms are included.

## Grid

M {5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5} × x/L {0.01, 0.05, 0.10, 0.15, 0.23, 0.30, 0.40, 0.50, 0.77}
× {normal, strike-slip, reverse} × the four metric/version combinations × `include_prob_zero` {False, True}.
Percentiles {0.16, 0.5, 0.84, −1 (mean)}. Displacements: 60 log-spaced values from 1e-3 to 20 m.

## Files

| File | Content | Feature |
|---|---|---|
| `displacements.csv` | the 60 displacements, full precision | — |
| `prob_exceed.csv` | P(D > d), one row per scenario, columns `p00`…`p59` | F1–F4 |
| `cdf.csv` | P(D ≤ d), same layout | F5 |
| `prob_exceed_out_of_range.csv` | P(D > d) at M 4.99 and 8.51 (fdhpy warns and extrapolates) | F13 |
| `displ_site.csv` | displacement at a percentile or the mean | F6, F7 |
| `displ_profile.csv` | profile with `xl_step = 0.1` | F8 |
| `displ_avg.csv` | average displacement, sum-of-principal, full rupture, median | F9 |
| `displ_max.csv` | maximum displacement, aggregate, full rupture | F10 |
| `params.csv` | `mu`, `sigma`, `p_gap`, `p_zero_slip`, `sigma_mu_agg` (`include_prob_zero=True`) | F11, F12, F16 |

Each file starts with `#` header lines that record the fdhpy, numpy, scipy and pandas versions. Read the files with

```python
pd.read_csv(path, comment="#", float_precision="round_trip")
```

pandas' default float parser does not round-trip the last bits of `displacements.csv`. Outputs are stored with
10 significant digits (relative rounding ≤ 5e-10).

## fdhpy values that are not test targets

These rows are frozen as fdhpy returns them, but by decision the pfdha implementation does not reproduce them:

| Rows | fdhpy behaviour | pfdha target |
|---|---|---|
| `displ_site`, `displ_profile`: full rupture, `include_prob_zero=True`, percentile ≠ −1 | the percentile displacement is multiplied by (1 − P) | the mixture percentile: the root of Eq. 31 exceedance = 1 − p (decision G3) |
| `displ_site`, `displ_profile`: individual segment, `include_prob_zero=True` | returns `None`, stored as NaN | see decision G2 |
| sum-of-principal, individual segment, `include_prob_zero=True` (all files) | the flag is ignored | see decision G2 |
| `params.csv`: `p_gap` for individual segment | 0 | — |

The mean rows (percentile −1) are targets: fdhpy integrates E[max(Y, 0)^(10/3)], the censored mean that Eq. 35 also
uses. The sum-of-principal maximum displacement is not frozen; it is not in LA23 (decision G5).
`displ_profile.csv` uses `xl_step = 0.1` because fdhpy returns `None` when the step does not divide 1.
