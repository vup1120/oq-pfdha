# Configuration guide for Lavrentiadis & Abrahamson (2023)

Reference: Lavrentiadis, G., & Abrahamson, N. A. (2023). Fault-displacement models for aggregate and principal
displacements. *Earthquake Spectra*, 41(4), 2806–2837. https://doi.org/10.1177/87552930231201531

The Lavrentiadis & Abrahamson (2023) model (LA23) predicts the probability of exceeding a fault displacement on normal,
strike-slip and reverse faults. It works with the power-transformed displacement D^0.3, which is normally
distributed, and it can include the probability of zero displacement (a site in a gap between segments, or zero
principal displacement where only distributed ruptures occur).

## Model selection keys

The model is exposed as **two classes, one per displacement definition** (the class choice is the definition, like
the `Petersen2011PrimaryFD_*` shape variants). Within a class, `output_type` selects the version: the full rupture
(the paper's "simplified FDM without segmentation") or an individual segment.

| Path | Required? | Allowed values | Purpose |
| --- | --- | --- | --- |
| `models.primary_surf_displ.type` | Yes | `Lavrentiadis2023PrimaryFD_aggregate` | The **aggregate** displacement (principal plus distributed): `output_type` `disp_agg_prime` (default) or `disp_agg_seg`. Aggregate chains run single-bucket; secondary-slot models are forbidden (FDLT-013). |
| `models.primary_surf_displ.type` | Yes | `Lavrentiadis2023PrimaryFD_principal` | The **sum-of-principal** displacement (no distributed ruptures): `output_type` `disp_prnc_prime` (default) or `disp_prnc_seg`. Not aggregate: secondary models remain legitimate. |

| `output_type` | Class | Version | Median of D^0.3 | Sigma of D^0.3 | Zero term with `include_zero_slip = true` |
| --- | --- | --- | --- | --- | --- |
| `disp_agg_prime` | `_aggregate` | full rupture | Eq. 22 | Eq. 24 | 1 − P(Gap), Eq. 25 |
| `disp_agg_seg` | `_aggregate` | individual segment | Eq. 14 | √(φ² + τ²), p. 16 | none: the flag is ignored with a warning |
| `disp_prnc_prime` | `_principal` | full rupture | Eq. 33 | Eq. 34 with φ_add | (1 − P(D_P = 0)) (1 − P(Gap)), Eq. 31 |
| `disp_prnc_seg` | `_principal` | individual segment | max(Eq. 14 + b₂, 0) | √(φ_P² + τ²) | 1 − P(D_P = 0), Eq. 32 with the Eq. 14 predictor |

## Scope and behavior

- **Supported fault styles:** normal, strike-slip, reverse.
- **Magnitude range:** M 5.0–8.5 (p. 24). Outside it the model logs one warning per magnitude and extrapolates.
- **Distribution:** D^0.3 is normal, left-truncated at zero (Eq. 17).
- **Slip component:** net displacement.
- **Classification:** aggregate (`Lavrentiadis2023PrimaryFD_aggregate`) / sum-of-principal
  (`Lavrentiadis2023PrimaryFD_principal`).

## Parameters (`models.primary_surf_displ.parameters`)

| Name | Type | Units | Default | Allowed | Required? | Description |
| --- | --- | --- | --- | --- | --- | --- |
| `style` | string | - | derived from the rupture rake | `"normal"`, `"strike-slip"`, `"reverse"` (case-insensitive) | No | Style of faulting. Any other value (including `"all"`) is an error. |
| `output_type` | string | - | `"disp_agg_prime"` / `"disp_prnc_prime"` | `_aggregate`: `"disp_agg_prime"`, `"disp_agg_seg"`. `_principal`: `"disp_prnc_prime"`, `"disp_prnc_seg"` | No | Version of the class's displacement definition: full rupture (`*_prime`) or individual segment (`*_seg`). A value of the other class is rejected (FDLT-015). |
| `include_zero_slip` | boolean | - | `false` | `true`, `false` | No | If `true`, the exceedance probability is scaled by the probability of a non-zero displacement (see the `output_type` table). fdhpy defaults to `true`. |

## Notes and cautions

- **Zero displacement.** A single segment has no gap (p. 21), so the individual-segment versions never use P(Gap).
  P(D_P = 0) is the probability of zero **principal** displacement where only distributed ruptures occur
  (Eqs. 31–32). It applies to the principal versions only: the aggregate already contains the distributed displacement.
- **Individual-segment versions on multi-section ruptures.** `disp_agg_seg` and `disp_prnc_seg` are functions of the
  position along one segment, X_seg/L_seg. For a `multiFaultSource` rupture the calculator measures x/L along the
  whole rupture, so a run that combines the two logs a warning. Use the full-rupture version there, or model the
  segments as separate sources.
- **Folding.** x/L is folded to [0, 0.5]; x/L and 1 − x/L give the same result.
- **Paper and code.** The regression follows the authors' reference implementation
  ([NHR3-UCLA/LA23_PFDHA_model](https://github.com/NHR3-UCLA/LA23_PFDHA_model), 2023-10-04). It differs from the printed
  paper in two places:
  - Eq. 34. The printed cross term 2ρ·φ_agg·φ_b2² is dimensionally inconsistent; the code uses 2ρ·φ_agg·φ_b2. Only the
    code form reproduces the IAEA TECDOC-2092 team L23 curves (the printed form fails Kumamoto, Le Teil and Norcia).
  - Eq. 32. For the full rupture the code predicts P(D_P = 0) from μ′_agg (Eq. 22), where the paper's text names μ_agg
    (Eq. 14). The differences are ≤ 0.004 for strike-slip and reverse faults, and up to 0.07 for normal faults near the
    rupture ends. The individual-segment version `disp_prnc_seg` uses μ_agg, as printed.
- **Differences from fdhpy 1.0.3.** `include_zero_slip` defaults to `false`. With the zero terms, a percentile
  displacement is the percentile of the mixture of zero and the power-normal distribution (fdhpy multiplies the
  percentile displacement by the non-zero probability). `disp_prnc_seg` applies P(D_P = 0) (fdhpy ignores the flag
  for segments). The maximum displacement is given for the aggregate full rupture only (Eqs. 38–39).

## Scenario (deterministic) interface

Both classes also provide the paper's deterministic quantities, for use outside the hazard calculation. Arguments
follow `get_prob` (`style`, `output_type`, `include_zero_slip`); `percentile = -1` gives the mean.

| Method | Returns | Source |
| --- | --- | --- |
| `get_prob(d, X_L_ratio, mag)` / `get_cdf(...)` | P(D > d) / P(D ≤ d) | Eqs. 17, 22–25, 31–34 |
| `get_displ_site(X_L_ratio, mag, percentile)` | displacement (m) at a percentile, or the mean | Eq. 17; mean as in Eq. 35 |
| `get_displ_profile(mag, percentile, xl_step)` | `(xl, displacement)` along the rupture, x/L from 0 to 1 | — |
| `get_displ_avg(mag)` | median average principal displacement (m); `_principal`, full rupture only | Eqs. 36–37, Table 5 |
| `get_displ_max(mag, percentile)` | maximum aggregate displacement (m); `_aggregate`, full rupture only | Eqs. 38–39, Table 6 |
| `get_stat_params(X_L_ratio, mag)` | median and sigma of D^0.3 | — |
| `get_p_gap(X_L_ratio, mag)`, `get_p_zero_slip(X_L_ratio, mag, segment=False)` | P(Gap); P(D_P = 0 \| not in a gap) | Eqs. 25, 32 |
| `get_sigma_mu_agg(mag)` | epistemic standard deviation of the median D_agg^0.3 | Eq. 29 |

## Verification

The classes are tested in `openquake/fdha/test/benchmark/fdhpy_comparison/FDHI_Tests/tests/test_lavrentiadis2023_parity.py`
against three references: the authors' implementation (copied unmodified to `FDHI_Tests/reference/la23_pfdha_model`,
agreement to 1e-12), golden fixtures frozen from fdhpy 1.0.3 (`FDHI_Tests/fixtures/la23`, 1e-6; the mean 1e-8), and
the paper's equations where pfdha departs from fdhpy. The IAEA benchmark (`test/benchmark/IAEA`, L23 entries) compares
the hazard curves with the TECDOC-2092 team curves.

## References

- Lavrentiadis, G., & Abrahamson, N. A. (2023). Fault-displacement models for aggregate and principal displacements.
  *Earthquake Spectra*, 41(4), 2806–2837. https://doi.org/10.1177/87552930231201531
- Sarmiento, A., et al. (2025). Comparisons of FDHI fault displacement models for principal and aggregate displacement.
  *Earthquake Spectra*, 41(4), 2691–2720. https://doi.org/10.1177/87552930251327894
