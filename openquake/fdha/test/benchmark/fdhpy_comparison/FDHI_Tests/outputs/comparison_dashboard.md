# PFDHA vs fdhpy Comparison Dashboard

**Generated:** 2025-11-27 12:16

## Overview

This report presents a comprehensive comparison between the `pfdha` (OpenQuake) 
implementation and the `fdhpy` (FDHI) reference implementation for Probabilistic 
Fault Displacement Hazard Analysis models.

The validation covers multiple magnitudes, normalized positions (x/L), and model variants. 
For each test case, exceedance probabilities are compared across a displacement grid, 
and quantitative metrics are computed to assess numerical agreement.

## Summary Table

| Model | Variant | Tests | Max Rel Diff | Mean Rel Diff | Pass Rate | Status |
|-------|---------|-------|--------------|---------------|-----------|--------|
| Chiou2025 | model7 | 6 | 5.00e-13 | 1.61e-13 | 100% | ✅ PASS |
| Chiou2025 | model8.2 | 6 | 3.44e-13 | 7.93e-14 | 100% | ✅ PASS |
| Kuehn2024 | normal | 9 | 1.11e-14 | 4.67e-15 | 100% | ✅ PASS |
| Kuehn2024 | reverse | 9 | 1.27e-14 | 2.94e-15 | 100% | ✅ PASS |
| Kuehn2024 | strike-slip | 9 | 1.32e-14 | 4.75e-15 | 100% | ✅ PASS |
| Lavrentiadis2023 | aggregate_no_zero | 9 | 2.49e-09 | 5.46e-10 | 100% | ✅ PASS |
| Lavrentiadis2023 | aggregate_with_zero | 9 | 2.49e-09 | 5.46e-10 | 100% | ✅ PASS |
| Moss2024 | d/ad_EQS | 6 | 1.95e-11 | 3.49e-12 | 100% | ✅ PASS |
| Moss2024 | d/ad_GIRS | 6 | 5.23e-11 | 1.10e-11 | 100% | ✅ PASS |
| Moss2024 | d/md_EQS | 6 | 2.51e-10 | 4.37e-11 | 100% | ✅ PASS |
| Moss2024 | d/md_GIRS | 6 | 3.35e-10 | 8.44e-11 | 100% | ✅ PASS |
| Petersen2011 | elliptical | 9 | 0.00e+00 | 0.00e+00 | 100% | ✅ PASS |
| Petersen2011 | quadratic | 9 | 0.00e+00 | 0.00e+00 | 100% | ✅ PASS |
| Youngs2003 | D/AD | 12 | 1.32e-05 | 2.67e-06 | 100% | ✅ PASS |

## Overall Statistics

- **Total test cases:** 111
- **Models validated:** 6
- **Tests passed:** 111 (100.0%)
- **Max relative difference (all):** 1.32e-05
- **Mean relative difference (all):** 2.89e-07

## Youngs2003

- **Variants tested:** D/AD
- **Total test cases:** 12
- **Tests passed:** 12 / 12
- **Max relative difference:** 1.32e-05

### Representative Hazard Curve Comparisons

![Youngs2003_D_AD_M7.0_XL0.25.png](plots/Youngs2003_D_AD_M7.0_XL0.25.png)

![Youngs2003_D_AD_M7.0_XL0.5.png](plots/Youngs2003_D_AD_M7.0_XL0.5.png)


## Petersen2011

- **Variants tested:** elliptical, quadratic
- **Total test cases:** 18
- **Tests passed:** 18 / 18
- **Max relative difference:** 0.00e+00

### Representative Hazard Curve Comparisons

![Petersen2011_elliptical_M7.0_XL0.3.png](plots/Petersen2011_elliptical_M7.0_XL0.3.png)

![Petersen2011_elliptical_M7.0_XL0.5.png](plots/Petersen2011_elliptical_M7.0_XL0.5.png)

![Petersen2011_quadratic_M7.0_XL0.3.png](plots/Petersen2011_quadratic_M7.0_XL0.3.png)

![Petersen2011_quadratic_M7.0_XL0.5.png](plots/Petersen2011_quadratic_M7.0_XL0.5.png)


## Moss2024

- **Variants tested:** d/ad_GIRS, d/ad_EQS, d/md_GIRS, d/md_EQS
- **Total test cases:** 24
- **Tests passed:** 24 / 24
- **Max relative difference:** 3.35e-10

### Representative Hazard Curve Comparisons

![Moss2024_d_ad_GIRS_M7.0_XL0.25.png](plots/Moss2024_d_ad_GIRS_M7.0_XL0.25.png)

![Moss2024_d_ad_GIRS_M7.0_XL0.5.png](plots/Moss2024_d_ad_GIRS_M7.0_XL0.5.png)

![Moss2024_d_ad_EQS_M7.0_XL0.25.png](plots/Moss2024_d_ad_EQS_M7.0_XL0.25.png)

![Moss2024_d_ad_EQS_M7.0_XL0.5.png](plots/Moss2024_d_ad_EQS_M7.0_XL0.5.png)

![Moss2024_d_md_GIRS_M7.0_XL0.25.png](plots/Moss2024_d_md_GIRS_M7.0_XL0.25.png)

![Moss2024_d_md_GIRS_M7.0_XL0.5.png](plots/Moss2024_d_md_GIRS_M7.0_XL0.5.png)

![Moss2024_d_md_EQS_M7.0_XL0.25.png](plots/Moss2024_d_md_EQS_M7.0_XL0.25.png)

![Moss2024_d_md_EQS_M7.0_XL0.5.png](plots/Moss2024_d_md_EQS_M7.0_XL0.5.png)


## Kuehn2024

- **Variants tested:** normal, reverse, strike-slip
- **Total test cases:** 27
- **Tests passed:** 27 / 27
- **Max relative difference:** 1.32e-14

### Representative Hazard Curve Comparisons

![Kuehn2024_normal_M7.0_XL0.3.png](plots/Kuehn2024_normal_M7.0_XL0.3.png)

![Kuehn2024_normal_M7.0_XL0.5.png](plots/Kuehn2024_normal_M7.0_XL0.5.png)

![Kuehn2024_reverse_M7.0_XL0.3.png](plots/Kuehn2024_reverse_M7.0_XL0.3.png)

![Kuehn2024_reverse_M7.0_XL0.5.png](plots/Kuehn2024_reverse_M7.0_XL0.5.png)

![Kuehn2024_strike-slip_M7.0_XL0.3.png](plots/Kuehn2024_strike-slip_M7.0_XL0.3.png)

![Kuehn2024_strike-slip_M7.0_XL0.5.png](plots/Kuehn2024_strike-slip_M7.0_XL0.5.png)


## Lavrentiadis2023

- **Variants tested:** aggregate_with_zero, aggregate_no_zero
- **Total test cases:** 18
- **Tests passed:** 18 / 18
- **Max relative difference:** 2.49e-09

### Representative Hazard Curve Comparisons

![Lavrentiadis2023_aggregate_with_zero_M7.0_XL0.3.png](plots/Lavrentiadis2023_aggregate_with_zero_M7.0_XL0.3.png)

![Lavrentiadis2023_aggregate_with_zero_M7.0_XL0.5.png](plots/Lavrentiadis2023_aggregate_with_zero_M7.0_XL0.5.png)

![Lavrentiadis2023_aggregate_no_zero_M7.0_XL0.3.png](plots/Lavrentiadis2023_aggregate_no_zero_M7.0_XL0.3.png)

![Lavrentiadis2023_aggregate_no_zero_M7.0_XL0.5.png](plots/Lavrentiadis2023_aggregate_no_zero_M7.0_XL0.5.png)


## Chiou2025

- **Variants tested:** model7, model8.2
- **Total test cases:** 12
- **Tests passed:** 12 / 12
- **Max relative difference:** 5.00e-13

### Representative Hazard Curve Comparisons

![Chiou2025_model7_M7.0_XL0.25.png](plots/Chiou2025_model7_M7.0_XL0.25.png)

![Chiou2025_model7_M7.0_XL0.5.png](plots/Chiou2025_model7_M7.0_XL0.5.png)

![Chiou2025_model8.2_M7.0_XL0.25.png](plots/Chiou2025_model8.2_M7.0_XL0.25.png)

![Chiou2025_model8.2_M7.0_XL0.5.png](plots/Chiou2025_model8.2_M7.0_XL0.5.png)


---

## Notes

- Tolerances vary by model: Youngs 2003 uses `rtol=1e-4` due to integration method differences; 
  all other models use `rtol=1e-6`.
- Relative differences are computed as `|pfdha - fdhpy| / |fdhpy|` where `|fdhpy| > 1e-12`.
- All comparisons use the same displacement grid as the automated pytest suite.
