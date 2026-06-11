# Figure 13 benchmark — Visini et al. (2025)

Reproduces the decision-tree example of Visini et al. (2025), *Empirical
regressions for distributed faulting of dip-slip earthquakes*, Earthquake
Spectra 41(4), doi:10.1177/87552930241308860 — Figure 13 — and compares the
computed conditional probabilities of exceedance against the published curves
digitized into `reference_data/visini2025_case{1,2,3}.csv`.

## Scenario (paper, "Suggestions for the application" section)

| Quantity | Value | Where it lives here |
| --- | --- | --- |
| Earthquake | Mw 7.0, normal, characteristic | `source_model_fig13.xml` |
| PF rupture length | 40 km (trace E–W at lat 39.633, dip 36° N) | `source_model_fig13.xml` |
| Site footprint | 100 m × 100 m | `corner_points` in `job_case*.ini` |
| Site position | hanging wall, r = 2000 m from the trace, opposite the PF midpoint (x/L = 0.5) | site centered at (16.2275, 39.650986) |
| TPFm | 1.42 m (value calculated in the paper for this scenario) | `tpfm = 1.42` in the FD logic-tree branches |
| Sigma truncation | ±3σ with truncation-range normalization (FDHLab `eps = 3`) | `n_sigma = 3` |
| Case 1 | Combinations A + B + C (Rank 1.5 fault beneath the site, plus nearby fault) | `case = case1`, traces `R1p5_local_A`, `R1p5_local_B` |
| Case 2 | Combinations A + B (nearby Rank 1.5 fault only) | `case = case2`, trace `R1p5_500m` |
| Case 3 | Combination A only (no Rank 1.5 fault within 1 km) | `case = case3`, trace `R1p5_far1` (~2.6 km) |

The geometry was verified against the calculator: the site context resolves to
r = 1.9995 km, x/L = 0.498, rx > 0 (hanging wall), L = 40.0 km.

## Two findings worth knowing about (established while reconstructing this benchmark)

1. **The Case 2 distance in the paper text (200 m) does not match the
   published figure.** Fitting the digitized Combination B curve with the
   paper's Equation 3 (σ = 1.0271, ±3σ truncation) gives a median normalized
   displacement of ≈0.17 m and an occurrence probability of ≈0.0084 —
   consistent with the Rank 1.5 fault at **500 m**
   (lnY = −1.80, P ≈ 0.149 × P_along), and inconsistent with 200 m
   (lnY = −1.61, P ≈ 0.266 × P_along). The reference implementation
   ([FDHLab](https://github.com/fault2shaESCWG/FDHLab),
   `Model/script_combB_pfdhcurves.m`) indeed uses `dist = [500]` and the
   `LOGISTIC_CombB_Mw7_SiteDim100_SiteDist500_HW` table. This benchmark
   therefore places the Case 2 trace at 500 m to match the figure.

2. **The reference occurrence factor is logistic × Monte Carlo.** In FDHLab
   the curves are `P_logistic(s) × P_montecarlo × exceedance`, where
   `P_montecarlo` is the along-strike intersection probability obtained by
   simulating DR segments (lognormal lengths; uniform and clustered
   placements averaged) over a total DR length = F-ratio × fault length.
   This pipeline reproduces the implied reference values
   (e.g. Combination A at 2 km: P_occ ≈ 0.0130 reference vs ≈ 0.0134 ± 0.0006
   here); the MC factor is not seeded, so each run carries ≈ ±5% noise.

## Running

```bash
# the three cases + comparison plot
python run_all_cases_and_plot.py

# pytest benchmark (slow; runs the CLI three times)
pytest test_fig13_reproduction.py -m slow
```

Current agreement against the digitized curves (model/reference ratio at the
calculation displacement levels, within the digitized range and below the
+3σ truncation cliff at ~2.5 m):

| Case | min | median | max |
| --- | --- | --- | --- |
| 1 (A+B+C) | 0.97 | 1.02 | 1.07 |
| 2 (A+B) | 0.81 | 0.87 | 0.97 |
| 3 (A) | 0.94 | 0.96 | 1.03 |

`case{1,2,3}_results.json` are snapshots of these runs (`imls` / `poes`).
