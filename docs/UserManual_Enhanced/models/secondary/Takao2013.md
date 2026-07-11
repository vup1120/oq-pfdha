# Configuration guide for Takao et al. (2013) — Secondary Displacement

Reference: Takao, M., Tsuchiyama, J., Annaka, T., & Kurita, T. (2013). Application of probabilistic fault displacement hazard analysis in Japan. Journal of Japan Association for Earthquake Engineering, 13(1), 17-36. https://doi.org/10.5610/jaee.13.17

The Takao et al. (2013) secondary (distributed) surface displacement model predicts displacement exceedance probabilities at off-fault locations, regressed on Japanese reverse- and strike-slip-faulting earthquakes. The distributed displacement is normalized by the maximum (PMD) or average (PAD) displacement of the principal fault; the 90% non-exceedance level decays exponentially with the closest distance from the principal fault (their Eqs. 15-16), and the conditional distribution is a Gamma distribution with shape a = 2.5 anchored at that level (their Eq. 17).

## Model selection keys

| Path | Required? | Allowed values | Purpose |
| --- | --- | --- | --- |
| `models.secondary_surf_displ.type` | Yes | `Takao2013SecondaryFD` | Chooses the Takao et al. (2013) secondary displacement exceedance model. |

## Scope and behavior

- **Supported fault styles:** Reverse and strike-slip faulting (Japanese data)
- **Displacement metric:** Distributed (secondary) displacement
- **Normalization:** DD/PMD or DD/PAD (maximum or average displacement on the principal fault)
- **Statistical distribution:** Gamma distribution (shape a = 2.5), scale anchored so the 90th percentile equals the Eq. 15/16 regression level
- **Slip component:** Net displacement
- **Distance range:** calibrated on data within ~0-20 km of the principal fault
- **Classification:** Single Principal & Distributed

## Parameters (`models.secondary_surf_displ.parameters`)

| Name | Type | Units | Default | Allowed | Required? | Description |
| --- | --- | --- | --- | --- | --- | --- |
| `norm_disp_type` | string | – | `"AD"` | `"AD"`, `"MD"` | No | Normalization displacement type. `"AD"` uses the average-displacement regression (Eq. 16, `DD/PAD = 1.9 exp(-0.17 r)`), `"MD"` the maximum-displacement regression (Eq. 15, `DD/PMD = 0.55 exp(-0.17 r)`). The paper's own case study uses `"AD"`. |
| `n_sigma` | float | – | `3.0` | > 0 | No | Truncation of the PMD/PAD lognormal at `mean ± n_sigma · sigma` in log10 space, mirroring `Takao2013PrimaryFD`. |

## Notes and cautions

- **Companion models:** Pair with `Takao2013PrimarySR` (P1p) and
  `Takao2013SecondarySR` (P2d, their Eq. 14) or `Takao2014SecondarySR`
  (pixel-size-dependent P2d refit) for the full Takao chain of their
  Eq. (13).
- **PMD/PAD scaling:** PMD uses the paper's own refit
  `log10(PMD) = -5.16 + 0.82 Mw` (their Eq. 9, constant 0.3 above Wells &
  Coppersmith 1994), PAD the Wells & Coppersmith relation
  `log10(PAD) = -4.80 + 0.69 Mw` (their Eq. 10); lognormal sigmas 0.42 and
  0.36 respectively.
- **Truncation and tails:** the deep hazard-curve tail (rates below
  ~1e-9/yr) is controlled by `n_sigma`; the benchmark chains use 5.0.
- **Validation:** `openquake/fdha/test/benchmark/takao_2013/` reproduces the
  paper's Fig. 11 case (b) distributed example (published anchor 7.0e-7/yr
  at 0.01 m reproduced to 0.2%).

## References

- Takao, M., Tsuchiyama, J., Annaka, T., & Kurita, T. (2013). Application of probabilistic fault displacement hazard analysis in Japan. Journal of Japan Association for Earthquake Engineering, 13(1), 17-36.
- Wells, D. L., & Coppersmith, K. J. (1994). New empirical relationships among magnitude, rupture length, rupture width, rupture area, and surface displacement. Bulletin of the Seismological Society of America, 84(4), 974-1002.
