# Figure 13 benchmark - Visini et al. (2025)

Reproduces the decision-tree example of Visini et al. (2025), *Empirical
regressions for distributed faulting of dip-slip earthquakes*, Earthquake
Spectra 41(4), doi:10.1177/87552930241308860 - Figure 13 - and compares the
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
| P(SR_primary) | **1.0 (pinned)** - Fig 13 shows *conditional* probabilities | `[FixedPrimarySR] value = 1.0` in the SR logic-tree branches |

The geometry was verified against the calculator: the site context resolves to
r = 1.9995 km, x/L = 0.498, rx > 0 (hanging wall), L = 40.0 km.

**Why P(SR_primary) is pinned to 1.** The paper's Figure 13 curves are
conditional probabilities of exceedance: Visini et al. (2025) explicitly
exclude both the earthquake rate and the primary surface-rupture probability
from the worked example. The hazard pipeline, however, deliberately gates the
distributed contribution with P(SR_primary) - the DR occurrence regressions
are fit on the SURE database, which contains only earthquakes with a mapped
Rank-1 surface rupture, so P_dist is conditional on the principal fault
reaching the surface and the gate converts it into a per-rupture rate
contribution (see `calc/hazard.py`). To compare like-with-like against the
published conditional curves, this benchmark therefore fixes the gate to 1
via `FixedPrimarySR`. (An earlier revision of these configs used
`Moss2013PrimarySR`, chosen before the gate existed; once the gate was
introduced it scaled all three curves by P_sr(Mw 7) ≈ 0.39 - a uniform
×0.4 offset against the figure.)

## Reading the truncation cliff in the comparison figures

Each case's computed curve ends in a sharp, smooth roll-off to exactly zero
rather than an extended power-law tail. This is by design:
`Visini2025SecondaryFD` implements the along-strike displacement as a
**truncated** log-normal (`n_sigma = 3`, matching the FDHLab MATLAB
reference's `eps = 3`), and a truncated distribution has zero density beyond
its truncation bound. The bound scales with each case's TPFm/combination and
falls at different displacements: ≈ 2.5 m for cases 2 and 3 (combination
A/A+B), ≈ 8-9 m for case 1 (A+B+C, a larger combined median). Evaluated at
the reference's own last digitized point (case 3, d = 2.2175 m, published
P ≈ 1.00e-5) the model gives ≈ 8.8e-6 - matching to within the same ~5-10%
agreement seen everywhere else on the curve.

`job_case*.ini`'s `displacement_measure_levels` jumps straight across each
of these gaps (`..., 1.0, 3.0, ...` and `..., 3.0, 5.0, 7.5, 10.0`), so a
naive plot of the raw calculation grid draws a straight line from the last
nonzero point to the hard zero - which on a log-y axis looks like a
near-vertical cliff and can read as a modelling error. It isn't: it's a
sampling artifact of the coarse grid, not a divergence from the reference.
`plot_fig13_cases_vs_reference.py` fixes this by densifying the
*plotted* curve across both gaps (see its module docstring) while leaving
the coarse grid used by `test_fig13_reproduction.py` and the `agreement`
stats below untouched.

## Two findings worth knowing about (established while reconstructing this benchmark)

1. **The Case 2 distance in the paper text (200 m) does not match the
   published figure.** Fitting the digitized Combination B curve with the
   paper's Equation 3 (σ = 1.0271, ±3σ truncation) gives a median normalized
   displacement of ≈0.17 m and an occurrence probability of ≈0.0084 -
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
# the three cases + comparison figures + ratio stats
# (writes ../Figures/visini2025_fig13_*.png and fig13_agreement.json)
PYTHONPATH=. python plot_fig13_cases_vs_reference.py

# legacy runner (three cases + single combined plot, results JSONs)
python run_all_cases_and_plot.py

# pytest benchmark (slow; runs the CLI three times)
pytest test_fig13_reproduction.py -m slow
```

The headline comparison figure is
`../Figures/visini2025_fig13_all_cases_ref_vs_impl.png` (published Fig. 13
curves dashed, oq-pfdha solid); per-case overlay and relative-error panels
sit next to it, and `fig13_agreement.json` records the computed/reference
ratio statistics of the latest run.

Current agreement against the digitized curves (model/reference ratio at the
calculation displacement levels, within the digitized range and below the
+3σ truncation cliff, d ≤ 2 m; run of 2026-07-10, see
`fig13_agreement.json` - the unseeded Monte Carlo along-strike factor moves
individual points by ≈ ±5% between runs):

| Case | min | median | max |
| --- | --- | --- | --- |
| 1 (A+B+C) | 0.97 | 1.03 | 1.07 |
| 2 (A+B) | 0.79 | 0.85 | 0.94 |
| 3 (A) | 0.89 | 0.91 | 0.99 |

`case{1,2,3}_results.json` are snapshots of earlier runs (`imls` / `poes`).
