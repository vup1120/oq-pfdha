# Configuration guide for Petersen et al. (2011)

Reference: Petersen, M. D., et al. (2011). Fault displacement hazard for strike-slip faults. Bulletin of the Seismological Society of America, 101(2), 805-825. https://doi.org/10.1785/0120100035

The Petersen et al. (2011) primary surface displacement model predicts displacement exceedance probabilities for strike-slip faults using a log-normal distribution in ln(cm) space. The model provides three functional forms for the along-strike displacement profile: bilinear, elliptical, and quadratic.

## Model selection keys

| Path | Required? | Allowed values | Purpose |
| --- | --- | --- | --- |
| `models.primary_surf_displ.type` | Yes | `Petersen2011PrimaryFD` | Chooses the Petersen (2011) displacement exceedance model. |

## Scope and behavior

- **Supported fault styles:** Strike-slip only
- **Magnitude range:** 6.0–8.0 (principal FDM)
- **Displacement metric:** Principal displacement
- **Statistical distribution:** Normal distribution in ln(cm) space
- **Slip component:** Lateral displacement
- **Classification:** Single Principal

## Parameters (`models.primary_surf_displ.parameters`)

| Name | Type | Units | Default | Allowed | Required? | Description |
| --- | --- | --- | --- | --- | --- | --- |
| `version` | string | – | `"quadratic"` | `"bilinear"`, `"elliptical"`, `"quadratic"` (case-insensitive) | No | Functional form for the along-strike displacement profile. `"bilinear"` uses piecewise linear (Eqns 7–9), `"elliptical"` uses elliptical profile (Eqn 13), `"quadratic"` uses quadratic profile with folding (Eqn 10). |

## Notes and cautions

- **Style limitation:** This model is designed specifically for strike-slip faults only. Using it for other fault styles may produce invalid results.
- **X_L_ratio validation:** X_L_ratio must be between 0 and 1. Values outside this range raise `ValueError`.
- **Version validation:** Invalid version strings raise `ValueError` with allowed options: "bilinear", "elliptical", "quadratic".
- **Units:** The model internally works in centimeters and log space (ln(cm)). All input displacements are converted from meters to centimeters.
- **Magnitude handling:** Magnitude can be scalar or array, but will be broadcast to match X_L_ratio shape if needed.
- **Vectorization:** The model fully supports vectorized inputs for `d` and `X_L_ratio`, with proper broadcasting to produce output of shape `(n_displacements, n_sites)`.
- **Bilinear intersection:** The bilinear version clips the intersection point X_L_prime to [0.25, 0.26] to ensure numerical stability.

## References

- Petersen, M. D., et al. (2011). Fault displacement hazard for strike-slip faults. Bulletin of the Seismological Society of America, 101(2), 805-825. https://doi.org/10.1785/0120100035
- Sarmiento, A., et al. (2025). Comparisons of FDHI fault displacement models for principal and aggregate displacement. Earthquake Spectra. https://doi.org/10.1177/87552930251327894




