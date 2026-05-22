# PFDHA vs fdhpy Comparison Report

**Date:** November 2025  
**Status:** ✅ All Models Validated

## Executive Summary

This report documents the comparison between `pfdha` (OpenQuake implementation) and `fdhpy` (reference implementation) for Probabilistic Fault Displacement Hazard Analysis models.

**All 6 models have been validated** with numerical differences within acceptable tolerances.

## Model Comparison Results

| Model | Tests | Status | Max Relative Diff |
|-------|-------|--------|-------------------|
| Youngs et al. (2003) | 12 | ✅ PASS | < 0.01% |
| Petersen et al. (2011) | 18 | ✅ PASS | < 1e-6 |
| Moss et al. (2024) | 24 | ✅ PASS | < 1e-6 |
| Kuehn et al. (2024) | 27 | ✅ PASS | < 1e-6 |
| Lavrentiadis & Abrahamson (2023) | 18 | ✅ PASS | < 1e-6 |
| Chiou et al. (2025) | 12 | ✅ PASS | < 1e-6 |

**Total: 111 test cases**

## Model Details

### Youngs et al. (2003)

- **Versions tested:** D/AD (normalized by average displacement)
- **Style:** "all" faulting style (Wells & Coppersmith 1994 coefficients)
- **Parameters:** M = 6.0-7.5, x/L = 0.1-0.5
- **Note:** Slight integration differences at M=6.0 (~0.001%), within tolerance

### Petersen et al. (2011)

- **Versions tested:** elliptical, quadratic
- **Parameters:** M = 6.5-7.5, x/L = 0.1-0.5
- **bilinear version:** Available in pfdha only (not in fdhpy)

### Moss et al. (2024)

- **Versions tested:** D/AD, D/MD
- **Sources:** GIRS, EQS
- **Parameters:** M = 6.5-7.5, x/L = 0.25-0.5

### Kuehn et al. (2024)

- **Styles tested:** normal, reverse, strike-slip
- **Parameters:** M = 7.0-7.6, x/L = 0.3-0.7
- **Folding:** Correctly implements x and (1-x) averaging

### Lavrentiadis & Abrahamson (2023)

- **Metric:** aggregate
- **Version:** full rupture
- **Parameters:** M = 6.5-7.5, x/L = 0.3-0.6
- **Zero slip:** Both include/exclude options tested

### Chiou et al. (2025)

- **Versions tested:** model7, model8.2
- **Parameters:** M = 6.5-7.5, x/L = 0.25-0.5

## Parameter Mapping Reference

| fdhpy Parameter | pfdha Parameter | Notes |
|-----------------|-----------------|-------|
| `magnitude` | `mag` | Earthquake magnitude |
| `xl` | `X_L_ratio` | Normalized position (array) |
| `version` | varies | Model-specific |
| `displ_array` | `d` | Displacement array |
| `use_girs` | `source="GIRS"/"EQS"` | Moss 2024 |
| `complete` | `completeness` | Moss 2024 |
| `include_prob_zero` | `include_zero_slip` | LA23 |
| `folded` | `folded` | Kuehn 2024 |
| `epistemic_uncertainty` | `epistemic_uncertainty` | Kuehn 2024 |

## Running Tests

```bash
# Navigate to the FDHI_Tests folder
cd /path/to/pfdha/openquake/fdha/test/FDHI_Tests

# Quick standalone test
python run_standalone.py

# Full pytest suite
python -m pytest tests/test_simple_reference.py -v
```

## Conclusion

The `pfdha` implementation has been validated against the `fdhpy` reference implementation. All models produce numerically equivalent results within scientific precision limits.
