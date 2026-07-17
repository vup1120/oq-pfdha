# Configuration guide for Lavrentiadis & Abrahamson (2023)

Reference: Lavrentiadis, G., & Abrahamson, N. A. (2023). Fault displacement model for aggregate and principal displacement. Earthquake Spectra.

The Lavrentiadis & Abrahamson (2023) primary surface displacement model predicts displacement exceedance probabilities for normal, strike-slip, and reverse faults. The model works in power-normal space (m^0.3) and can account for zero-slip probability and rupture gaps.

## Model selection keys

The model is exposed as **two classes, one per displacement definition**
(the class choice IS the definition, like the `Petersen2011PrimaryFD_*`
shape variants):

| Path | Required? | Allowed values | Purpose |
| --- | --- | --- | --- |
| `models.primary_surf_displ.type` | Yes | `Lavrentiadis2023PrimaryFD_aggregate` | The **aggregate**-definition variants (`output_type` `disp_agg_prime`, default, or `disp_agg_seg`). Aggregate chains run single-bucket; secondary-slot models are forbidden (FDLT-013). |
| `models.primary_surf_displ.type` | Yes | `Lavrentiadis2023PrimaryFD_principal` | The **sum-of-principal** `disp_prnc_prime` variant; `output_type` is pinned by the class. Not aggregate: secondary models remain legitimate. |

## Scope and behavior

- **Supported fault styles:** Normal, strike-slip, reverse
- **Magnitude range:** 5.0–8.5
- **Displacement metrics:** 
  - Aggregate displacement for entire event rupture (`disp_agg_prime`) — `Lavrentiadis2023PrimaryFD_aggregate`
  - Aggregate displacement for single segment (`disp_agg_seg`) — `Lavrentiadis2023PrimaryFD_aggregate`
  - Principal displacement for entire event rupture (`disp_prnc_prime`) — `Lavrentiadis2023PrimaryFD_principal`
- **Statistical distribution:** Normal distribution in power-normal space (m^0.3)
- **Slip component:** Net displacement
- **Classification:** Aggregate (`Lavrentiadis2023PrimaryFD_aggregate`) / Sum of principal (`Lavrentiadis2023PrimaryFD_principal`)

## Parameters (`models.primary_surf_displ.parameters`)

| Name | Type | Units | Default | Allowed | Required? | Description |
| --- | --- | --- | --- | --- | --- | --- |
| `style` | string | – | `"normal"` | `"normal"`, `"strike-slip"`, `"reverse"` (case-insensitive) | No | Style of faulting. |
| `output_type` | string | – | `"disp_agg_prime"` | `"disp_agg_prime"`, `"disp_agg_seg"` (`Lavrentiadis2023PrimaryFD_aggregate` only) | No | Aggregate metric to evaluate. `disp_prnc_prime` is **rejected** on this class (FDLT-015): select `Lavrentiadis2023PrimaryFD_principal` instead. On `Lavrentiadis2023PrimaryFD_principal` the metric is pinned by the class and passing any explicit `output_type` is an error. |
| `include_zero_slip` | boolean | – | `false` | `true`, `false` | No | If `true`, the probability accounts for zero slip and gap probabilities. If `false`, uses only the displacement distribution. Accepted by both classes. |

## Notes and cautions

- **Output type validation:** Invalid `output_type` values raise `ValueError` with allowed options; `disp_prnc_prime` on the aggregate class (or any explicit `output_type` on the `_principal` class) raises a `ValueError` naming the correct class, and is also rejected at logic-tree validation time (FDLT-015).
- **Style handling:** Style is converted to lowercase for matching. Accepted values are "normal", "strike-slip", "reverse".
- **Surface rupture length:** The model internally sets `srl=1.0` when calling the slip profile function, as the model normalizes by rupture length.
- **Power-normal space:** All internal calculations use the 0.3 power transformation. This is a key feature of the model.
- **Vectorization:** The model supports vectorized inputs for `d` and `X_L_ratio`, with proper broadcasting and output shape normalization.
- **Zero-slip inclusion:** When `include_zero_slip=true`, the probabilities are reduced to account for the possibility of no displacement. This is important for hazard calculations.
- **Gap probability:** The gap probability accounts for along-strike variations in rupture continuity and is magnitude and position dependent.

## References

- Lavrentiadis, G., & Abrahamson, N. A. (2023). Fault displacement model for aggregate and principal displacement. Earthquake Spectra.
- Sarmiento, A., et al. (2025). Comparisons of FDHI fault displacement models for principal and aggregate displacement. Earthquake Spectra. https://doi.org/10.1177/87552930251327894




