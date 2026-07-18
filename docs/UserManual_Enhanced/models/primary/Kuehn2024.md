# Configuration guide for Kuehn et al. (2024)

Reference: Kuehn, N. M., et al. (2024). Primary surface fault displacement models for probabilistic fault displacement hazard analysis. Earthquake Spectra.

The Kuehn et al. (2024) primary surface displacement model predicts aggregate displacement exceedance probabilities for normal, reverse, and strike-slip faults using a Box-Cox transformed normal distribution with optional epistemic uncertainty via posterior coefficient sampling.

## Model selection keys

| Path | Required? | Allowed values | Purpose |
| --- | --- | --- | --- |
| `models.primary_surf_displ.type` | Yes | `Kuehn2024PrimaryFD` | Chooses the Kuehn (2024) displacement exceedance model. |

## Scope and behavior

- **Supported fault styles:** Normal, reverse, strike-slip
- **Magnitude ranges:**
  - Strike-slip: 6.0–8.0
  - Reverse: 5.0–8.0
  - Normal: (see coefficient data)
- **Displacement metric:** Aggregate displacement
- **Statistical distribution:** Normal distribution in transformed space (Box-Cox transformation)
- **Classification:** Aggregate
- **Slip component:** Net displacement

## Parameters (`models.primary_surf_displ.parameters`)

| Name | Type | Units | Default | Allowed | Required? | Description |
| --- | --- | --- | --- | --- | --- | --- |
| `style` | string | - | - | `"normal"`, `"reverse"`, `"strike-slip"` (case-insensitive) | Yes | Fault style. Invalid values raise ValueError. |
| `folded` | boolean | - | `true` | `true`, `false` | No | If true, averages probabilities from x/L and 1−x/L positions. If false, uses only the site position. |
| `epistemic_uncertainty` | boolean | - | `true` | `true`, `false` | No | If true, uses ensemble of coefficient realizations (posterior sampling). If false, uses mean coefficients only. |
| `coefficient_type` | string | - | - | `"full"`, `"mean"` | No | Alternative to `epistemic_uncertainty`. `"full"` enables epistemic uncertainty, `"mean"` uses mean coefficients. If provided, overrides `epistemic_uncertainty`. |

## Notes and cautions

- **Style validation:** Only `"normal"`, `"reverse"`, and `"strike-slip"` are accepted (case-insensitive). Invalid values raise `ValueError`.
- **Magnitude:** Only single scalar magnitude values are allowed. Array inputs raise `ValueError`.
- **Epistemic uncertainty output shape:** When `epistemic_uncertainty=true`, the output includes an additional dimension for coefficient realizations. Shape is `(n_models, n_displ, n_sites)` or reduced appropriately.
- **Coefficient data:** The model requires coefficient CSV files to be present in the data directory. Missing files will raise errors during coefficient loading.
- **Vectorization:** The model supports vectorized inputs for `d` and `X_L_ratio`, with proper broadcasting.
- **Output shape normalization:** For single sites, the model may return reduced dimensions for backward compatibility with tests.

## References

- Kuehn, N. M., et al. (2024). Primary surface fault displacement models for probabilistic fault displacement hazard analysis. Earthquake Spectra.
- Sarmiento, A., et al. (2025). Comparisons of FDHI fault displacement models for principal and aggregate displacement. Earthquake Spectra. https://doi.org/10.1177/87552930251327894




