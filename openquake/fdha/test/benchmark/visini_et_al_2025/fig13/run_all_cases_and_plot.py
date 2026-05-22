"""
Run all three cases (1, 2, 3) and create comparison plot with reference values.
Similar to all_cases_total_ref_vs_impl_colorcases_axes_v8.png
"""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import subprocess
import json
import os
import sys

# Configuration files (using INI format)
base_dir = Path(__file__).parent
config_files = {
    'case1': base_dir / "job_case1.ini",
    'case2': base_dir / "job_case2.ini",
    'case3': base_dir / "job_case3.ini"
}

# Reference CSV files (relative to this script's directory)
# Place reference data in a 'reference_data' subdirectory or update paths as needed
ref_data_dir = base_dir / "reference_data"
ref_files = {
    'case1': ref_data_dir / "visini2025_case1.csv",
    'case2': ref_data_dir / "visini2025_case2.csv",  # Updated to match actual filename
    'case3': ref_data_dir / "visini2025_case3.csv"
}

# Output files
output_jsons = {
    'case1': base_dir / "case1_results_after_fix.json",
    'case2': base_dir / "case2_results_after_fix.json",
    'case3': base_dir / "case3_results_after_fix.json"
}

print("=" * 80)
print("Running All Cases (1, 2, 3) using INI configuration files")
print("=" * 80)

results = {}

# Run each case
for case_name, config_file in config_files.items():
    print(f"\n{'='*80}")
    print(f"Running {case_name.upper()}...")
    print(f"{'='*80}")
    print(f"Config: {config_file}")
    
    output_json = output_jsons[case_name]
    
    cmd = [
        sys.executable, "-m", "openquake.fdha.main",
        str(config_file.name),
        "--output", str(output_json.name),
    ]
    
    try:
        # cwd = fig13 so rank1p5_traces.xml, logic trees, and source resolve like CLI users
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=str(base_dir),
        )
        
        if result.returncode != 0:
            print(f"❌ Error running {case_name}:")
            print(result.stderr)
            continue
        
        with open(output_json, 'r') as f:
            result_data = json.load(f)

        displacements = None
        probabilities = None
        if 'imls' in result_data and 'poes' in result_data:
            displacements = np.array(result_data['imls'])
            if isinstance(result_data['poes'][0], list):
                probabilities = np.array(result_data['poes'][0])
            else:
                probabilities = np.array(result_data['poes'])
        elif result_data.get('logic_tree') and result_data.get('aggregate_hazard_csv'):
            ag_path = Path(result_data['aggregate_hazard_csv'])
            with open(ag_path, newline='') as f:
                reader = csv.DictReader(f)
                d0_list, mean_list = [], []
                for row in reader:
                    d0_list.append(float(row['D0']))
                    mean_list.append(float(row['mean']))
            displacements = np.array(d0_list)
            probabilities = np.array(mean_list)
        else:
            print(f"⚠️  Unexpected format for {case_name}")
            continue
        
        results[case_name] = {
            'displacements': displacements,
            'probabilities': probabilities
        }
        
        print(f"✅ {case_name.upper()} completed: {len(displacements)} points")
        print(f"   Probability range: {np.min(probabilities):.6f} to {np.max(probabilities):.6f}")
        
    except Exception as e:
        print(f"❌ Error processing {case_name}: {e}")
        import traceback
        traceback.print_exc()
        continue

# Load reference data
print(f"\n{'='*80}")
print("Loading Reference Data")
print(f"{'='*80}")

ref_data = {}
for case_name, ref_file in ref_files.items():
    if os.path.exists(ref_file):
        try:
            # Try loading with skiprows=1 in case there's a header
            try:
                ref_arr = np.loadtxt(ref_file, delimiter=',', skiprows=1)
            except (ValueError, IndexError):
                # If that fails, try without skiprows
                ref_arr = np.loadtxt(ref_file, delimiter=',')
            
            ref_data[case_name] = {
                'displacements': ref_arr[:, 0],
                'probabilities': ref_arr[:, 1]
            }
            print(f"✅ Loaded {case_name} reference: {len(ref_arr)} points")
        except Exception as e:
            print(f"⚠️  Error loading {case_name} reference: {e}")
    else:
        print(f"⚠️  Reference file not found: {ref_file}")

# Create comparison plot
print(f"\n{'='*80}")
print("Creating Comparison Plot")
print(f"{'='*80}")

fig, ax = plt.subplots(figsize=(9, 6))

# Color scheme: Case 1 (reds), Case 2 (blues), Case 3 (greens)
colors = {
    'case1': {'ref': '#8B0000', 'impl': '#DC143C'},  # Dark red, Crimson
    'case2': {'ref': '#00008B', 'impl': '#4169E1'},  # Dark blue, Royal blue
    'case3': {'ref': '#006400', 'impl': '#32CD32'}  # Dark green, Lime green
}

# Plot each case
for case_name in ['case1', 'case2', 'case3']:
    if case_name not in results:
        continue
    
    impl_d = results[case_name]['displacements']
    impl_p = results[case_name]['probabilities']
    
    # Plot implementation (solid line with triangle marker)
    ax.loglog(impl_d, impl_p, 
             color=colors[case_name]['impl'],
             linestyle='-',
             linewidth=2,
             marker='^',
             markersize=5,
             label=f"Case {case_name[-1]} (our work)",
             alpha=0.9)
    
    # Plot reference if available (dashed line, thicker)
    if case_name in ref_data:
        ref_d = ref_data[case_name]['displacements']
        ref_p = ref_data[case_name]['probabilities']
        
        ax.loglog(ref_d, ref_p,
                 color=colors[case_name]['ref'],
                 linestyle='--',
                 linewidth=3,
                 label=f"Case {case_name[-1]} (reference)",
                 alpha=0.8)

# Formatting
ax.set_xlabel('Displacement (m)', fontsize=12)
ax.set_ylabel('Probability of Exceedance', fontsize=12)
ax.set_xlim([0.01, 11])
ax.legend(fontsize=12, loc='best')
ax.grid(True, alpha=0.3, which='both')

# Remove title (as per previous requirements)
plt.tight_layout()
plt.title("Visini et al., 2025 model validation")

# Save plot
output_plot = base_dir / "all_cases_total_ref_vs_impl_colorcases_axes_v9.png"
plt.savefig(output_plot, dpi=150, bbox_inches='tight')
print(f"\n✅ Plot saved to: {output_plot}")

plt.close()

# Print summary comparison
print(f"\n{'='*80}")
print("Summary Comparison at Key Displacements")
print(f"{'='*80}")

key_displacements = [0.01, 0.05, 0.1, 0.5, 1.0, 3.0, 5.0, 10.0]

for case_name in ['case1', 'case2', 'case3']:
    if case_name not in results:
        continue
    
    print(f"\n{case_name.upper()}:")
    print(f"{'Displacement (m)':>18} | {'Reference':>12} | {'Implementation':>15} | {'Ratio':>8}")
    print("-" * 60)
    
    impl_d = results[case_name]['displacements']
    impl_p = results[case_name]['probabilities']
    
    for d in key_displacements:
        # Interpolate implementation
        impl_idx = np.argmin(np.abs(impl_d - d))
        impl_val = impl_p[impl_idx] if impl_p[impl_idx] > 0 else np.nan
        
        # Interpolate reference if available
        if case_name in ref_data:
            ref_d = ref_data[case_name]['displacements']
            ref_p = ref_data[case_name]['probabilities']
            ref_idx = np.argmin(np.abs(ref_d - d))
            ref_val = ref_p[ref_idx] if ref_p[ref_idx] > 0 else np.nan
            
            ratio = impl_val / ref_val if not np.isnan(ref_val) and not np.isnan(impl_val) and ref_val > 0 else np.nan
            ratio_str = f"{ratio:.2f}" if not np.isnan(ratio) else "N/A"
            ref_str = f"{ref_val:.6f}" if not np.isnan(ref_val) else "N/A"
        else:
            ref_str = "N/A"
            ratio_str = "N/A"
        
        impl_str = f"{impl_val:.6f}" if not np.isnan(impl_val) else "N/A"
        
        print(f"{d:>18.2f} | {ref_str:>12} | {impl_str:>15} | {ratio_str:>8}")

print("\n" + "=" * 80)
print("All cases completed!")
print("=" * 80)
