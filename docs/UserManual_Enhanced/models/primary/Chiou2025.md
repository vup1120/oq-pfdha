# Configuration guide for Chiou et al. (2025)

Reference: Chiou, B. S. J., et al. (2025). Fault displacement model for surface principal rupture of strike-slip faults. Earthquake Spectra. https://doi.org/10.1177/87552930251337703

The Chiou et al. (2025) primary surface displacement model predicts sum-of-principal displacement exceedance probabilities for strike-slip faults using a negative exponentially modified Gaussian (NEMG) distribution.

## Model selection keys

| Path | Required? | Allowed values | Purpose |
| --- | --- | --- | --- |
| `models.primary_surf_displ.type` | Yes | `Chiou2025PrimaryFD` | Chooses the Chiou (2025) displacement exceedance model. |

## Scope and behavior

- **Supported fault styles:** Strike-slip only
- **Magnitude range:** 6.0–8.3 (recommended). Magnitudes outside this range emit a warning but computation proceeds.
- **Displacement metric:** Sum-of-principal displacement
- **Statistical distribution:** Negative exponentially modified Gaussian (NEMG)
- **Predicted variable:** ln(D) where D is displacement in meters
- **Classification:** Sum of Principal (supersedes Petersen et al. 2011 for strike-slip)

## Parameters (`models.primary_surf_displ.parameters`)

| Name | Type | Units | Default | Allowed | Required? | Description |
| --- | --- | --- | --- | --- | --- | --- |
| `style` | string | - | `"strike-slip"` | `"strike-slip"`, `"strikeslip"`, `"ss"` (case-insensitive) | No | Fault style. Must be strike-slip; otherwise raises ValueError. |
| `version` | string | - | `"model7"` | `"model7"`, `"model8.1"`, `"model8.2"`, `"model8.3"` (case-insensitive) | No | Model formulation version. Coefficients are loaded from `chiou_2025_coefficients.csv`. |

### Runtime inputs

- `d`: Target displacement threshold(s) in meters (sourced from `parameters.target_displacement` at run-time)
- `X_L_ratio`: Normalized position x/L in [0, 1] (computed from rupture geometry at run-time)
- `mag`: Moment magnitude Mw (sourced from each rupture in the input source model at run-time)

## Implementation details

### Position transformation

The model uses a symmetric folding of x/L to [0, 0.5] using modulus reflection:
```python
r = x_L_ratio - np.floor(x_L_ratio)  # map to [0,1)
folded_x = 0.5 - np.abs(r - 0.5)  # symmetric fold to [0,0.5]
```

An elliptical x* coordinate is computed:
```python
x_star = sqrt(max(0, 1 - 4 * (x_L_ratio - 0.5)^2))
```

### Magnitude scaling

The magnitude scaling function f_M(M) uses a smooth transition:
```python
f_M = m2 * (M - m3) + (m2 - m1) / cn * log(0.5 * (1 + exp(-cn * (M - m3))))
```

### Mean and aleatory variability

- Mean: `μ = c0 + f_M + c1 * (x_star - 1.0)`
- Aleatory variability includes magnitude-dependent and position-dependent components:
  - `σ_mag = max(0.4, cv1 * exp(cv2 * max(0, M - 6.1)))`
  - `σ_xl = cv3 * exp(cv4 * max(0, folded_x - ccap))`
  - `σ' = sqrt(σ_mag^2 + σ_xl^2)`

### Distribution

The model uses a negative exponentially modified Gaussian (NEMG) distribution:
- Exponential mixing parameter: `ν = cv5`
- Shape parameter: `K = ν / σ'`
- Exceedance probability computed using `scipy.stats.exponnorm.cdf(-ln(d); loc=-μ, scale=σ', K)`

### Coefficients

Model coefficients are loaded from `openquake/fdha/primary_surf_displ/data/chiou_2025_coefficients.csv`. Each version (model7, model8.1, model8.2, model8.3) has its own coefficient set.

## Notes and cautions

- **Style validation:** Only strike-slip faulting is supported. Providing any other style raises a `ValueError`.
- **Magnitude warnings:** Magnitudes outside the recommended range (6.0, 8.3) emit a warning but computation proceeds.
- **Version validation:** Invalid version strings raise a `ValueError` with a list of available versions.
- **X_L_ratio clipping:** Input X_L_ratio values are automatically clipped to [0, 1] before processing.
- **Displacement units:** All displacement values are in meters. The model internally works in log space (ln(D)).
- **Vectorization:** The model supports vectorized inputs for `d` and `X_L_ratio`, returning exceedance probabilities with appropriate broadcasting.

## References

- Chiou, B. S. J., et al. (2025). Fault displacement model for surface principal rupture of strike-slip faults. Earthquake Spectra. https://doi.org/10.1177/87552930251337703
- Sarmiento, A., et al. (2025). Comparisons of FDHI fault displacement models for principal and aggregate displacement. Earthquake Spectra. https://doi.org/10.1177/87552930251327894




