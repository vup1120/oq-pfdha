# Configuration guide for Youngs et al. (2003)

Reference: Youngs, R. R., et al. (2003). A methodology for probabilistic fault displacement hazard analysis (PFDHA). Earthquake Spectra, 19(1), 191-219.

The Youngs et al. (2003) primary surface displacement model predicts displacement exceedance probabilities for normal faults using normalized displacement (d/AD or d/MD) with Gamma or Beta distributions. The model uses Wells & Coppersmith (1994) magnitude-displacement scaling relationships and numerical integration over epsilon space.

## Model selection keys

| Path | Required? | Allowed values | Purpose |
| --- | --- | --- | --- |
| `models.primary_surf_displ.type` | Yes | `Youngs2003PrimaryFD` | Chooses the Youngs (2003) displacement exceedance model. |

## Scope and behavior

- **Supported fault styles:** Normal faulting
- **Magnitude range:** (varies by database, typically 4.5–7.6)
- **Displacement metric:** Principal displacement
- **Normalization:** d/AD (average displacement) or d/MD (maximum displacement)
- **Statistical distribution:** 
  - Gamma distribution for d/AD
  - Beta distribution for d/MD
- **Slip component:** Vertical displacement
- **Classification:** Single Principal

## Parameters (`models.primary_surf_displ.parameters`)

| Name | Type | Units | Default | Allowed | Required? | Description |
| --- | --- | --- | --- | --- | --- | --- |
| `style` | string | - | - | `"all"`, `"normal"` (case-insensitive) | Yes | Faulting style for Wells & Coppersmith (1994) coefficients. `"all"` uses "All styles" coefficients (recommended, consistent with paper and fdhpy). `"normal"` uses "Normal faulting" specific coefficients. |
| `norm_disp_type` | string | - | - | `"AD"`, `"MD"` | Yes | Normalization displacement type. `"AD"` uses average displacement normalization (Gamma distribution), `"MD"` uses maximum displacement normalization (Beta distribution). |

## Notes and cautions

- **Style limitation:** This model is designed specifically for normal faulting. The `style` parameter selects which Wells & Coppersmith (1994) coefficient set to use, not the fault mechanism.
- **Style validation:** Only `"all"` and `"normal"` are accepted (case-insensitive). Invalid values raise `ValueError`.
- **Normalization type validation:** Only `"AD"` and `"MD"` are accepted. Invalid values raise `ValueError`.
- **Magnitude validation:** Magnitude must be a scalar value. Array inputs raise `ValueError`.
- **Numerical integration:** The model uses epsilon-based numerical integration (matching fdhpy implementation) to capture total aleatory variability. The integration uses ±3 sigma truncation with 0.1 step size.
- **Truncation correction:** For `norm_disp_type="MD"`, the model applies a truncation correction to account for the constraint that D/MD ≤ 1. This correction divides the CDF by `CDF(1; α, β)` for values where `y ≤ 1`.
- **Vectorization:** The model supports vectorized inputs for `d` and `X_L_ratio`, returning output of shape `(n_displacements, n_sites)`.
- **Wells & Coppersmith scaling:** The model uses Wells & Coppersmith (1994) relationships for magnitude-displacement scaling. The "all" style coefficients are recommended for consistency with the original paper and fdhpy reference implementation.

## References

- Youngs, R. R., et al. (2003). A methodology for probabilistic fault displacement hazard analysis (PFDHA). Earthquake Spectra, 19(1), 191-219.
- Wells, D. L., & Coppersmith, K. J. (1994). New empirical relationships among magnitude, rupture length, rupture width, rupture area, and surface displacement. Bulletin of the Seismological Society of America, 84(4), 974-1002.
- Sarmiento, A., et al. (2025). Comparisons of FDHI fault displacement models for principal and aggregate displacement. Earthquake Spectra. https://doi.org/10.1177/87552930251327894




