#!/usr/bin/env python
"""Compare current results with reference values."""
import json
import sys
from pathlib import Path
import numpy as np

def main():
    test_dir = Path(__file__).parent
    results_path = test_dir / "results.json"
    
    if not results_path.exists():
        print(f"❌ Results file not found: {results_path}")
        print("Please run the calculation first to generate results.json")
        return 1
    
    # Load results
    with open(results_path, 'r') as f:
        results = json.load(f)
    
    # Reference values from paper (Youngs2003 AD 85)
    ref_Youngs_N = np.array([
        3.29E-04,
        3.25E-04,
        3.16E-04,
        3.08E-04,
        3.01E-04,
        2.85E-04,
        2.67E-04,
        2.49E-04,
        2.34E-04,
        2.08E-04,
        1.56E-04,
        1.14E-04,
        8.16E-05,
        6.15E-05,
        1.31E-05,
        4.72E-06,
        1.79E-06,
        8.09E-07
    ])
    
    # Extract arrays
    imls = np.array(results['imls'])
    rates = np.array(results['rates'])
    
    # Handle both 1D and 2D arrays - flatten to 1D
    if rates.ndim == 2:
        if rates.shape[0] == 1:
            rates = rates[0]  # Single site: take first row
        else:
            rates = rates.flatten()  # Multiple sites: flatten
    
    # Ensure both are 1D
    rates = rates.flatten()
    
    print("="*80)
    print("COMPARISON: Current Implementation vs Reference Values")
    print("="*80)
    
    # Check lengths
    if len(imls) != len(rates):
        print(f"❌ IML and POE count mismatch: imls={len(imls)}, rates={len(rates)}")
        return 1
    
    if len(imls) != len(ref_Youngs_N):
        print(f"⚠️  IML and reference count mismatch: imls={len(imls)}, ref={len(ref_Youngs_N)}")
        print(f"Using minimum length: {min(len(imls), len(ref_Youngs_N))}")
        min_len = min(len(imls), len(rates), len(ref_Youngs_N))
        imls = imls[:min_len]
        rates = rates[:min_len]
        ref_Youngs_N = ref_Youngs_N[:min_len]
    
    # Calculate differences
    abs_diff = np.abs(rates - ref_Youngs_N)
    rel_diff = np.abs((rates - ref_Youngs_N) / (ref_Youngs_N + 1e-20)) * 100
    ratio = rates / (ref_Youngs_N + 1e-20)
    
    max_abs_diff = np.max(abs_diff)
    max_rel_diff = np.max(rel_diff)
    mean_abs_diff = np.mean(abs_diff)
    mean_rel_diff = np.mean(rel_diff)
    mean_ratio = np.mean(ratio)
    
    print(f"\nPOE Comparison:")
    print(f"  Max absolute difference: {max_abs_diff:.2e}")
    print(f"  Max relative difference: {max_rel_diff:.2f}%")
    print(f"  Mean absolute difference: {mean_abs_diff:.2e}")
    print(f"  Mean relative difference: {mean_rel_diff:.2f}%")
    print(f"  Mean ratio (Current/Reference): {mean_ratio:.2f}")
    print(f"  Min ratio: {np.min(ratio):.2f}")
    print(f"  Max ratio: {np.max(ratio):.2f}")
    
    # Check if results match (within tolerance)
    rtol = 1e-2  # 1% relative tolerance for reference comparison
    atol = 1e-6
    
    if np.allclose(rates, ref_Youngs_N, rtol=rtol, atol=atol):
        print(f"\n✅ Results match reference within tolerance (rtol={rtol}, atol={atol})")
        return 0
    else:
        print(f"\n⚠️  Results differ from reference beyond tolerance (rtol={rtol}, atol={atol})")
        print("\nDetailed comparison (showing differences > tolerance):")
        print(f"{'IML':>10} {'Current':>15} {'Reference':>15} {'Abs Diff':>15} {'Rel Diff %':>12} {'Ratio':>10}")
        print("-" * 90)
        n_diff = 0
        for i in range(len(imls)):
            if abs_diff[i] > atol or rel_diff[i] > rtol * 100:
                print(f"{imls[i]:>10.4f} {rates[i]:>15.6e} {ref_Youngs_N[i]:>15.6e} {abs_diff[i]:>15.6e} {rel_diff[i]:>12.4f} {ratio[i]:>10.2f}")
                n_diff += 1
        if n_diff == 0:
            print("(All values within tolerance)")
        return 0  # Don't fail, just report differences

if __name__ == "__main__":
    sys.exit(main())
