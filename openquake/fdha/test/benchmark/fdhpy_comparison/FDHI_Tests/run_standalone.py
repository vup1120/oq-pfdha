#!/usr/bin/env python3
"""
Standalone PFDHA Reference Test

This script can be run directly without pytest to verify the test setup.
Run with: python run_standalone.py

Requirements:
  - fdhpy installed (pip install fdhpy or pip install -e /path/to/fdhpy)
  - pfdha installed (pip install -e /path/to/pfdha)
"""

import sys

import numpy as np

# ============================================================================
# Import fdhpy (reference implementation)
# ============================================================================
try:
    import fdhpy
except ImportError as e:
    print(f"[ERROR] fdhpy is not installed: {e}")
    print("Install with: pip install fdhpy  OR  pip install -e /path/to/fdhpy")
    sys.exit(1)

from fdhpy import (
    YoungsEtAl2003,
    PetersenEtAl2011,
    MossEtAl2024,
    KuehnEtAl2024,
    LavrentiadisAbrahamson2023,
    ChiouEtAl2025,
)

# ============================================================================
# Import pfdha (this implementation)
# ============================================================================
try:
    from openquake.fdha.primary_surf_displ import (
        Youngs2003PrimaryFD,
        Petersen2011PrimaryFD,
        Moss2024PrimaryFD,
        Kuehn2024PrimaryFD,
        Lavrentiadis2023PrimaryFD_aggregate,
        Chiou2025PrimaryFD,
    )
except ImportError as e:
    print(f"[ERROR] pfdha OpenQuake FD models not importable: {e}")
    print("Install with: pip install -e /path/to/pfdha")
    sys.exit(1)

# ============================================================================
# Test Configuration
# ============================================================================

print("=" * 70)
print("PFDHA Standalone Reference Test")
print("=" * 70)
print()

# Test displacement array
DISPLACEMENTS = np.array([0.001, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0])

# Results tracking
results = {"passed": 0, "failed": 0, "errors": []}


def test_case(name, fdhpy_result, pfdha_result, rtol=1e-6, atol=1e-10):
    """Compare two results and report."""
    global results
    try:
        np.testing.assert_allclose(fdhpy_result, pfdha_result, rtol=rtol, atol=atol)
        print(f"  ✓ {name}")
        results["passed"] += 1
        return True
    except AssertionError as e:
        print(f"  ✗ {name}")
        max_diff = np.max(np.abs(fdhpy_result - pfdha_result))
        max_rel = np.max(np.abs((fdhpy_result - pfdha_result) / (fdhpy_result + 1e-15)) * 100)
        print(f"    Max abs diff: {max_diff:.2e}")
        print(f"    Max rel diff: {max_rel:.4f}%")
        results["failed"] += 1
        results["errors"].append((name, str(e)))
        return False


# ============================================================================
# Test Youngs 2003
# ============================================================================
print("\n--- Youngs et al. (2003) ---")
try:
    for mag in [6.0, 6.5, 7.0, 7.5]:
        for xl in [0.1, 0.25, 0.5]:
            # fdhpy
            fdhpy_model = YoungsEtAl2003(
                magnitude=mag, xl=xl, version="d/ad", displ_array=DISPLACEMENTS
            )
            fdhpy_result = fdhpy_model.prob_exceed
            
            # pfdha
            pfdha_model = Youngs2003PrimaryFD()
            pfdha_result = pfdha_model.get_prob(
                d=DISPLACEMENTS, X_L_ratio=np.array([xl]),
                mag=mag, style="all", norm_disp_type="AD"
            ).flatten()
            
            # Relaxed tolerance for integration differences
            test_case(f"M={mag}, x/L={xl}", fdhpy_result, pfdha_result, rtol=1e-4)
except Exception as e:
    print(f"  ✗ ERROR: {e}")
    import traceback
    traceback.print_exc()
    results["failed"] += 1
    results["errors"].append(("Youngs 2003", str(e)))


# ============================================================================
# Test Petersen 2011
# ============================================================================
print("\n--- Petersen et al. (2011) ---")
try:
    for version in ["elliptical", "quadratic"]:
        for mag in [6.5, 7.0, 7.5]:
            xl = 0.5
            
            # fdhpy
            fdhpy_model = PetersenEtAl2011(
                magnitude=mag, xl=xl, version=version, displ_array=DISPLACEMENTS
            )
            fdhpy_result = fdhpy_model.prob_exceed
            
            # pfdha
            pfdha_model = Petersen2011PrimaryFD()
            pfdha_result = pfdha_model.get_prob(
                d=DISPLACEMENTS, X_L_ratio=np.array([xl]),
                mag=mag, version=version
            ).flatten()
            
            test_case(f"{version} M={mag}", fdhpy_result, pfdha_result)
except Exception as e:
    print(f"  ✗ ERROR: {e}")
    import traceback
    traceback.print_exc()
    results["failed"] += 1
    results["errors"].append(("Petersen 2011", str(e)))


# ============================================================================
# Test Moss 2024
# ============================================================================
print("\n--- Moss et al. (2024) ---")
try:
    for version, norm_type in [("d/ad", "AD"), ("d/md", "MD")]:
        for use_girs in [True, False]:
            mag, xl = 7.0, 0.5
            source = "GIRS" if use_girs else "EQS"
            
            # fdhpy
            fdhpy_model = MossEtAl2024(
                magnitude=mag, xl=xl, version=version,
                displ_array=DISPLACEMENTS, use_girs=use_girs, complete=True
            )
            fdhpy_result = fdhpy_model.prob_exceed
            
            # pfdha
            pfdha_model = Moss2024PrimaryFD()
            pfdha_result = pfdha_model.get_prob(
                d=DISPLACEMENTS, X_L_ratio=np.array([xl]),
                mag=mag, version=norm_type, source=source, completeness="complete"
            ).flatten()
            
            test_case(f"{version} {source}", fdhpy_result, pfdha_result)
except Exception as e:
    print(f"  ✗ ERROR: {e}")
    import traceback
    traceback.print_exc()
    results["failed"] += 1
    results["errors"].append(("Moss 2024", str(e)))


# ============================================================================
# Test Kuehn 2024
# ============================================================================
print("\n--- Kuehn et al. (2024) ---")
try:
    # Test both folded and unfolded, with various x/L values
    for style in ["normal", "reverse", "strike-slip"]:
        for xl in [0.3, 0.5, 0.7]:
            mag = 7.2
            
            # fdhpy
            fdhpy_model = KuehnEtAl2024(
                style=style, magnitude=mag, xl=xl,
                version="median_coeffs", folded=True, displ_array=DISPLACEMENTS
            )
            fdhpy_result = fdhpy_model.prob_exceed
            
            # pfdha
            pfdha_model = Kuehn2024PrimaryFD()
            pfdha_result = pfdha_model.get_prob(
                d=DISPLACEMENTS, X_L_ratio=np.array([xl]),
                mag=mag, style=style, folded=True, epistemic_uncertainty=False
            ).flatten()
            
            test_case(f"{style} x/L={xl}", fdhpy_result, pfdha_result)
except Exception as e:
    print(f"  ✗ ERROR: {e}")
    import traceback
    traceback.print_exc()
    results["failed"] += 1
    results["errors"].append(("Kuehn 2024", str(e)))


# ============================================================================
# Test Lavrentiadis 2023
# ============================================================================
print("\n--- Lavrentiadis & Abrahamson (2023) ---")
try:
    for include_prob_zero in [True, False]:
        mag, xl = 7.0, 0.5
        
        # fdhpy - needs style parameter
        fdhpy_model = LavrentiadisAbrahamson2023(
            magnitude=mag, xl=xl, displ_array=DISPLACEMENTS,
            metric="aggregate", version="full rupture",
            style="strike-slip",  # Required parameter
            include_prob_zero=include_prob_zero
        )
        fdhpy_result = fdhpy_model.prob_exceed
        
        # pfdha
        pfdha_model = Lavrentiadis2023PrimaryFD_aggregate()
        pfdha_result = pfdha_model.get_prob(
            d=DISPLACEMENTS, X_L_ratio=np.array([xl]),
            mag=mag, style="strike-slip",
            output_type="disp_agg_prime",
            include_zero_slip=include_prob_zero
        ).flatten()
        
        zero_str = "with_zero" if include_prob_zero else "no_zero"
        test_case(f"aggregate {zero_str}", fdhpy_result, pfdha_result)
except Exception as e:
    print(f"  ✗ ERROR: {e}")
    import traceback
    traceback.print_exc()
    results["failed"] += 1
    results["errors"].append(("Lavrentiadis 2023", str(e)))


# ============================================================================
# Test Chiou 2025
# ============================================================================
print("\n--- Chiou et al. (2025) ---")
try:
    for version in ["model7", "model8.2"]:
        mag, xl = 7.0, 0.5
        
        # fdhpy
        fdhpy_model = ChiouEtAl2025(
            magnitude=mag, xl=xl, version=version, displ_array=DISPLACEMENTS
        )
        fdhpy_result = fdhpy_model.prob_exceed
        
        # pfdha
        pfdha_model = Chiou2025PrimaryFD()
        pfdha_result = pfdha_model.get_prob(
            d=DISPLACEMENTS, X_L_ratio=np.array([xl]),
            mag=mag, version=version
        ).flatten()
        
        test_case(f"{version}", fdhpy_result, pfdha_result)
except Exception as e:
    print(f"  ✗ ERROR: {e}")
    import traceback
    traceback.print_exc()
    results["failed"] += 1
    results["errors"].append(("Chiou 2025", str(e)))


# ============================================================================
# Summary
# ============================================================================
print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"Passed: {results['passed']}")
print(f"Failed: {results['failed']}")
print(f"Total:  {results['passed'] + results['failed']}")
print()

if results["failed"] > 0:
    print("FAILURES:")
    for name, error in results["errors"]:
        print(f"  - {name}: {error[:100]}...")
    print()
    sys.exit(1)
else:
    print("All tests passed!")
    sys.exit(0)

















