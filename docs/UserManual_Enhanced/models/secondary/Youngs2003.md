# Configuration guide for Youngs et al. (2003) — Secondary Displacement

Reference: Youngs, R. R., et al. (2003). A methodology for probabilistic fault displacement hazard analysis (PFDHA). Earthquake Spectra, 19(1), 191-219.

The Youngs et al. (2003) secondary (distributed) surface displacement model predicts displacement exceedance probabilities for normal faults at off-fault locations. The model uses normalized displacement (D/MD) with a Gamma distribution and Wells & Coppersmith (1994) scaling for maximum displacement.

## Model selection keys

| Path | Required? | Allowed values | Purpose |
| --- | --- | --- | --- |
| `models.secondary_surf_displ.type` | Yes | `Youngs2003SecondaryFD` | Chooses the Youngs (2003) secondary displacement exceedance model. |

## Scope and behavior

- **Supported fault styles:** Normal faulting
- **Magnitude range:** (varies, typically 5.5–7.4)
- **Displacement metric:** Distributed (secondary) displacement
- **Normalization:** D/MD (maximum displacement on principal fault)
- **Statistical distribution:** Gamma distribution
- **Slip component:** Vertical displacement
- **Distance range:** 0 ≤ r ≤ 15 km
- **Classification:** Single Principal & Distributed

## Parameters (`models.secondary_surf_displ.parameters`)

| Name | Type | Units | Default | Allowed | Required? | Description |
| --- | --- | --- | --- | --- | --- | --- |
| `percentile` | string | – | `"85"` | `"85"`, `"95"` | No | Percentile used for scaling the distributed displacement relative to maximum displacement. Different percentiles use different scaling factors. |

## Notes and cautions

- **Style limitation:** This model is designed specifically for normal faulting. Using it for other fault styles may produce invalid results.
- **Percentile validation:** Only `"85"` and `"95"` are accepted. Invalid values raise `ValueError`.
- **Distance range:** The model is calibrated for distances 0 ≤ r ≤ 15 km. Results beyond this range may be less reliable.
- **Hanging wall/footwall:** The model distinguishes between hanging wall (rx > 0) and footwall (rx ≤ 0) locations, using different scaling factors for each.
- **Numerical integration:** The model uses 100-point discretization for numerical integration over the maximum displacement distribution. This is an internal implementation detail and not user-configurable.
- **Truncation:** The MD distribution is truncated at ±3 sigma to ensure numerical stability in the integration.
- **Gamma shape parameter:** The Gamma distribution uses a constant shape parameter `a = 2.5` for all sites and distances.
- **Vectorization:** The model supports vectorized inputs for `d`, `rx`, and `r`, returning output of shape `(n_sites, n_displacements)`.
- **Wells & Coppersmith scaling:** The model uses Wells & Coppersmith (1994) "Normal faulting" coefficients for maximum displacement scaling.
- **Caching:** The model caches normalization factors for magnitude values (rounded to 2 decimal places) to improve performance.

## References

- Youngs, R. R., et al. (2003). A methodology for probabilistic fault displacement hazard analysis (PFDHA). Earthquake Spectra, 19(1), 191-219.
- Wells, D. L., & Coppersmith, K. J. (1994). New empirical relationships among magnitude, rupture length, rupture width, rupture area, and surface displacement. Bulletin of the Seismological Society of America, 84(4), 974-1002.
- Sarmiento, A., et al. (2025). Comparisons of FDHI fault displacement models for principal and aggregate displacement. Earthquake Spectra. https://doi.org/10.1177/87552930251327894




