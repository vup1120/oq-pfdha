#!/usr/bin/env python
"""Plot Norcia hazard curve comparison with reference values."""
import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def load_results(json_path):
    """Load results from JSON file."""
    with open(json_path, 'r') as f:
        data = json.load(f)
    return data

def plot_comparison():
    """Plot comparison of current results with reference values."""
    test_dir = Path(__file__).parent
    
    # Load current results
    results = load_results(test_dir / "results.json")
    
    # Extract data
    imls = np.array(results['imls'])
    poes = np.array(results['poes'])
    
    # Handle both 1D and 2D arrays
    if poes.ndim == 2:
        poes = poes[0] if poes.shape[0] == 1 else poes.flatten()
    poes = poes.flatten()
    
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
    
    # Ensure same length
    min_len = min(len(imls), len(poes), len(ref_Youngs_N))
    imls_common = imls[:min_len]
    poes_common = poes[:min_len]
    ref_common = ref_Youngs_N[:min_len]
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    
    # Plot 1: Hazard Curve Comparison
    ax1.loglog(imls_common, poes_common, 'r--s', linewidth=2, markersize=6, 
               label='Current Implementation', alpha=0.8)
    ax1.loglog(imls_common, ref_common, 'g-^', linewidth=2, markersize=6, 
               label='Reference (Paper - Youngs2003 AD 85)', alpha=0.8)
    
    ax1.set_xlabel('Displacement, d (m)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Annual Exceedance Rate, λ(d) (yr⁻¹)', fontsize=12, fontweight='bold')
    ax1.set_title('Hazard Curve Comparison\nNorcia Sensitivity - Youngs2003 AD 85', 
                  fontsize=14, fontweight='bold', pad=15)
    ax1.legend(loc='upper right', fontsize=10, framealpha=0.9)
    ax1.grid(True, which='major', linestyle='-', alpha=0.3, linewidth=0.8)
    ax1.grid(True, which='minor', linestyle='--', alpha=0.2, linewidth=0.5)
    
    # Plot 2: Ratio
    ratio = poes_common / (ref_common + 1e-20)
    
    ax2.semilogx(imls_common, ratio, 'b-o', linewidth=2, markersize=6, 
                 alpha=0.8, label='Current / Reference')
    ax2.axhline(y=1.0, color='k', linestyle='--', linewidth=1, alpha=0.5, label='Ratio = 1.0')
    
    ax2.set_xlabel('Displacement, d (m)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Ratio (Current / Reference)', fontsize=12, fontweight='bold')
    ax2.set_title('Ratio Comparison', fontsize=14, fontweight='bold', pad=15)
    ax2.legend(loc='upper left', fontsize=10, framealpha=0.9)
    ax2.grid(True, which='major', linestyle='-', alpha=0.3, linewidth=0.8)
    ax2.grid(True, which='minor', linestyle='--', alpha=0.2, linewidth=0.5)
    
    # Calculate statistics
    abs_diff = np.abs(poes_common - ref_common)
    rel_diff = np.abs((poes_common - ref_common) / (ref_common + 1e-20)) * 100
    
    stats_text = (
        f"Statistics (Current vs Reference):\n"
        f"Max abs diff: {np.max(abs_diff):.2e}\n"
        f"Max rel diff: {np.max(rel_diff):.2f}%\n"
        f"Mean abs diff: {np.mean(abs_diff):.2e}\n"
        f"Mean rel diff: {np.mean(rel_diff):.2f}%\n"
        f"Mean ratio: {np.mean(ratio):.2f}\n"
        f"Min ratio: {np.min(ratio):.2f}\n"
        f"Max ratio: {np.max(ratio):.2f}"
    )
    
    ax2.text(0.02, 0.98, stats_text,
             transform=ax2.transAxes,
             fontsize=9,
             verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    
    # Save plot
    output_path = test_dir / "hazard_curve_comparison.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✅ Plot saved to: {output_path}")
    
    plt.close()

if __name__ == "__main__":
    plot_comparison()
