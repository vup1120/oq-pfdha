# Configuration guide for Petersen et al. (2011) — Secondary Displacement

Reference: Petersen, M. D., et al. (2011). Fault displacement hazard for strike-slip faults. Bulletin of the Seismological Society of America, 101(2), 805-825. https://doi.org/10.1785/0120100035

The Petersen et al. (2011) secondary (distributed) surface displacement model predicts displacement exceedance probabilities for strike-slip faults at off-fault locations. The model uses a log-normal distribution in ln(cm) space and is limited to distances within 2.5 km of the principal fault.

## Model selection keys

| Path | Required? | Allowed values | Purpose |
| --- | --- | --- | --- |
| `models.secondary_surf_displ.type` | Yes | `Petersen2011SecondaryFD` | Chooses the Petersen (2011) secondary displacement exceedance model. |

## Scope and behavior

- **Supported fault styles:** Strike-slip
- **Magnitude range:** 6–8 (recommended for strike-slip faults)
- **Displacement metric:** Distributed (secondary) displacement
- **Statistical distribution:** Normal distribution in ln(cm) space
- **Slip component:** Lateral displacement
- **Distance range:** 0 ≤ r ≤ 2.5 km (limited to triggered ruptures, no triggered ruptures included)
- **Classification:** Single Principal & Distributed

## Parameters (`models.secondary_surf_displ.parameters`)

| Name | Type | Units | Default | Allowed | Required? | Description |
| --- | --- | --- | --- | --- | --- | --- |
| `pixel_size` | integer | meters | `25` | `25`, `50`, `100`, `150`, `200` | No | Size of the pixel ("cell" in the paper) for distributed rupture modeling. Different pixel sizes have different regression coefficients. The former name `cell_size` is still accepted as a deprecated alias. |

**Note:** While the pixel_size parameter is accepted, the current implementation uses a single regression equation (Page 818, Eqn 18) for all pixel sizes. The pixel size-specific coefficients from Table 4 are defined in the code but not currently used in the regression calculation.

## Notes and cautions

- **Style limitation:** This model is designed specifically for strike-slip faults only. Using it for other fault styles may produce invalid results.
- **Magnitude validation:** Magnitude must be between 6 and 8 for strike-slip faults. Values outside this range raise `ValueError`.
- **Distance range:** The model is calibrated for distances 0 ≤ r ≤ 2.5 km. Results beyond this range may be less reliable. The model does not include triggered ruptures.
- **Pixel size:** The `pixel_size` parameter is accepted (25, 50, 100, 150, or 200 m) but the current implementation uses a single regression equation regardless of pixel size. Invalid pixel sizes may raise errors. (`cell_size` remains a deprecated alias.)
- **Units:** The model internally works in centimeters and log space (ln(cm)). All input displacements are converted from meters to centimeters, and distances are converted from kilometers to meters.
- **Zero distance handling:** Distances of exactly zero are replaced with 0.1 m to avoid numerical issues with log(0).
- **Vectorization:** The model fully supports vectorized inputs for `d` and `r`, with proper broadcasting to produce output of shape `(n_sites, n_displacements)`.
- **Output shape:** For single sites or single displacements, the model may return reduced dimensions (1D arrays) for backward compatibility.
- **Wells & Coppersmith scaling:** The model references Wells & Coppersmith (1994) for strike-slip average displacement, though the regression directly predicts ln(d) rather than normalized displacement.

## References

- Petersen, M. D., et al. (2011). Fault displacement hazard for strike-slip faults. Bulletin of the Seismological Society of America, 101(2), 805-825. https://doi.org/10.1785/0120100035
- Wells, D. L., & Coppersmith, K. J. (1994). New empirical relationships among magnitude, rupture length, rupture width, rupture area, and surface displacement. Bulletin of the Seismological Society of America, 84(4), 974-1002.
- Sarmiento, A., et al. (2025). Comparisons of FDHI fault displacement models for principal and aggregate displacement. Earthquake Spectra. https://doi.org/10.1177/87552930251327894




