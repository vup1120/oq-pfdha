"""
Test Case 1 after segment length bounds fix.
Compare results with reference values.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../../..'))

import subprocess
import json

# Configuration file
config_file = Path(__file__).parent / "config_fig13_case1.toml"
output_dir = Path(__file__).parent
output_json = output_dir / "case1_results_after_fix.json"

# Reference CSV file (place in reference_data subdirectory or update path)
ref_csv = str(Path(__file__).parent / "reference_data" / "visini2025_case1.csv")

print("=" * 80)
print("Running Case 1 Test with Segment Length Bounds Fix")
print("=" * 80)
print(f"Config file: {config_file}")
print(f"Output directory: {output_dir}")

# Run calculation using CLI
try:
    cmd = [
        "python", "-m", "openquake.fdha.main",
        str(config_file),
        "--output", str(output_json)
    ]
    print(f"\nRunning: {' '.join(cmd)}")
    # Use repository root (6 levels up from this script's directory)
    repo_root = Path(__file__).parent.parent.parent.parent.parent.parent
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(repo_root))
    
    if result.returncode != 0:
        print(f"Error running calculation:")
        print(result.stderr)
        raise RuntimeError(f"Calculation failed with return code {result.returncode}")
    
    print("Calculation completed successfully!")
    
    # Load results from JSON
    with open(output_json, 'r') as f:
        result_data = json.load(f)
    
    # Extract results - adjust based on actual JSON structure
    if 'imls' in result_data and 'poes' in result_data:
        # Standard format: imls (displacements) and poes (probabilities)
        displacements = np.array(result_data['imls'])
        # poes is a list of lists, take the first curve
        if isinstance(result_data['poes'][0], list):
            probabilities = np.array(result_data['poes'][0])
        else:
            probabilities = np.array(result_data['poes'])
    elif 'hazard_curves' in result_data:
        # OpenQuake format
        hc = result_data['hazard_curves'][0]
        displacements = np.array([p['poe'] for p in hc['poes']])
        probabilities = np.array([p['poe'] for p in hc['poes']])
    elif 'displacements' in result_data:
        # Custom format
        displacements = np.array(result_data['displacements'])
        probabilities = np.array(result_data['probabilities'])
    else:
        # Try to find the data
        displacements = np.array(result_data.get('displacements', []))
        probabilities = np.array(result_data.get('probabilities', []))
    
    if len(displacements) == 0 or len(probabilities) == 0:
        print(f"\n⚠️  Warning: Could not extract results from JSON. Structure:")
        print(json.dumps(result_data, indent=2)[:500])
        raise ValueError("Could not extract displacement/probability data from results")
    
    print(f"\nResults:")
    print(f"  Number of displacement points: {len(displacements)}")
    print(f"  Probability range: {np.min(probabilities):.6f} to {np.max(probabilities):.6f}")
    
    # Load reference
    if os.path.exists(ref_csv):
        ref_data = np.loadtxt(ref_csv, delimiter=',')
        ref_d = ref_data[:, 0]
        ref_p = ref_data[:, 1]
        
        print(f"\nReference data loaded:")
        print(f"  Number of points: {len(ref_d)}")
        print(f"  Displacement range: {np.min(ref_d):.4f} to {np.max(ref_d):.4f} m")
        print(f"  Probability range: {np.min(ref_p):.6f} to {np.max(ref_p):.6f}")
        
        # Create comparison plot
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Plot reference
        ax.loglog(ref_d, ref_p, 'k--', linewidth=3, label='Reference (Case 1)', alpha=0.8)
        
        # Plot implementation
        ax.loglog(displacements, probabilities, 'r-', linewidth=2, marker='o', 
                 markersize=4, label='Implementation (after fix)', alpha=0.8)
        
        ax.set_xlabel('Displacement (m)', fontsize=12)
        ax.set_ylabel('Probability of Exceedance', fontsize=12)
        ax.set_title('Case 1: Comparison After Segment Length Bounds Fix', fontsize=13)
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_xlim([0.0001, 11])
        
        plt.tight_layout()
        output_plot = output_dir / "case1_after_segment_bounds_fix.png"
        plt.savefig(output_plot, dpi=150, bbox_inches='tight')
        print(f"\nPlot saved to: {output_plot}")
        
        # Compare at key points
        print("\n" + "=" * 80)
        print("Comparison at Key Displacement Values")
        print("=" * 80)
        print(f"{'Displacement (m)':>18} | {'Reference':>12} | {'Implementation':>15} | {'Ratio':>8}")
        print("-" * 60)
        
        key_displacements = [0.01, 0.05, 0.1, 0.5, 1.0, 3.0, 5.0, 10.0]
        for d in key_displacements:
            # Interpolate reference
            ref_idx = np.argmin(np.abs(ref_d - d))
            ref_val = ref_p[ref_idx] if ref_p[ref_idx] > 0 else np.nan
            
            # Interpolate implementation
            impl_idx = np.argmin(np.abs(displacements - d))
            impl_val = probabilities[impl_idx] if probabilities[impl_idx] > 0 else np.nan
            
            ratio = impl_val / ref_val if not np.isnan(ref_val) and not np.isnan(impl_val) and ref_val > 0 else np.nan
            
            ratio_str = f"{ratio:.2f}" if not np.isnan(ratio) else "N/A"
            ref_str = f"{ref_val:.6f}" if not np.isnan(ref_val) else "N/A"
            impl_str = f"{impl_val:.6f}" if not np.isnan(impl_val) else "N/A"
            
            print(f"{d:>18.2f} | {ref_str:>12} | {impl_str:>15} | {ratio_str:>8}")
        
        plt.close()
    else:
        print(f"\n⚠️  Reference file not found: {ref_csv}")
        print("  Creating plot without reference comparison...")
        
        # Plot without reference
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.loglog(displacements, probabilities, 'r-', linewidth=2, marker='o', 
                 markersize=4, label='Implementation (after fix)')
        ax.set_xlabel('Displacement (m)', fontsize=12)
        ax.set_ylabel('Probability of Exceedance', fontsize=12)
        ax.set_title('Case 1: Results After Segment Length Bounds Fix', fontsize=13)
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_xlim([0.0001, 11])
        plt.tight_layout()
        output_plot = output_dir / "case1_after_segment_bounds_fix.png"
        plt.savefig(output_plot, dpi=150, bbox_inches='tight')
        print(f"Plot saved to: {output_plot}")
        plt.close()
    
    # Save results to JSON
    import json
    results_dict = {
        'displacements': displacements.tolist(),
        'probabilities': probabilities.tolist()
    }
    results_file = output_dir / "case1_results_after_segment_bounds_fix.json"
    with open(results_file, 'w') as f:
        json.dump(results_dict, f, indent=2)
    print(f"Results saved to: {results_file}")
    
except Exception as e:
    print(f"\n❌ Error running calculation: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\n" + "=" * 80)
print("Test completed!")
print("=" * 80)
