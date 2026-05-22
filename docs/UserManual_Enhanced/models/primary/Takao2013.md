# Configuration guide for Takao et al. (2013)

Reference: Takao, M., et al. (2013). Development of probabilistic fault displacement hazard analysis method for reverse and strike-slip faults. Journal of Japan Association for Earthquake Engineering, 13(4), 1-20.

The Takao et al. (2013) primary surface displacement model predicts displacement exceedance probabilities for reverse and strike-slip faults using normalized displacement (d/AD or d/MD) with Gamma or Beta distributions. The model parameters depend on surface rupture length.

## Model selection keys

| Path | Required? | Allowed values | Purpose |
| --- | --- | --- | --- |
| `models.primary_surf_displ.type` | Yes | `Takao2013PrimaryFD` | Chooses the Takao (2013) displacement exceedance model. |

## Scope and behavior

- **Supported fault styles:** Reverse, strike-slip
- **Magnitude range:** 5.8–7.4
- **Displacement metric:** Principal displacement
- **Normalization:** d/AD (average displacement) or d/MD (maximum displacement)
- **Statistical distribution:** 
  - Gamma distribution for d/AD
  - Beta distribution for d/MD
- **Slip component:** Net displacement
- **Classification:** Single Principal

## Parameters (`models.primary_surf_displ.parameters`)

| Name | Type | Units | Default | Allowed | Required? | Description |
| --- | --- | --- | --- | --- | --- | --- |
| `norm_disp_type` | string | – | – | `"AD"`, `"MD"` | Yes | Normalization displacement type. `"AD"` uses average displacement normalization (Gamma distribution), `"MD"` uses maximum displacement normalization (Beta distribution). |

## Notes and cautions

- **Style limitation:** This model is designed for reverse and strike-slip faults only. Using it for normal faults may produce invalid results.
- **Magnitude range:** The model is calibrated for magnitudes 5.8–7.4. Results outside this range may be less reliable.
- **Surface rupture length:** The model automatically estimates SRL from magnitude using Wells & Coppersmith (1994). The distribution parameters (Gamma or Beta) depend on whether SRL < 10 km or SRL ≥ 10 km.
- **Normalization type validation:** Only `"AD"` and `"MD"` are accepted. Invalid values raise `ValueError`.
- **Numerical integration:** The model uses 1000-point discretization for numerical integration. This is an internal implementation detail and not user-configurable.
- **Truncation:** The AD/MD distribution is truncated at ±3 sigma to ensure numerical stability in the integration.
- **Vectorization:** The model supports vectorized inputs for `d` and `X_L_ratio`, returning output of shape `(n_displacements, n_sites)`.
- **Wells & Coppersmith scaling:** The model uses Wells & Coppersmith (1994) "All styles" coefficients for both AD and MD scaling relationships.

## References

- Takao, M., et al. (2013). Development of probabilistic fault displacement hazard analysis method for reverse and strike-slip faults. Journal of Japan Association for Earthquake Engineering, 13(4), 1-20.
- Wells, D. L., & Coppersmith, K. J. (1994). New empirical relationships among magnitude, rupture length, rupture width, rupture area, and surface displacement. Bulletin of the Seismological Society of America, 84(4), 974-1002.
- Sarmiento, A., et al. (2025). Comparisons of FDHI fault displacement models for principal and aggregate displacement. Earthquake Spectra. https://doi.org/10.1177/87552930251327894




