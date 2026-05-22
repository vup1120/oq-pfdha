#!/usr/bin/env python3
"""
Generate comprehensive comparison plot for Mammarella et al. (2024) model.

Compares pfdha implementation against CPSR.m (MATLAB/Octave) reference
across all 3 MSRs and 3 fault styles (9 combinations).
"""

import csv
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path

from openquake.fdha.primary_surf_rup.mammarella2024 import MammarellaEtAl2024PrimarySR


# Paths
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / "outputs"
# Reference data (place in reference_data subdirectory or update path)
REFERENCE_CSV = SCRIPT_DIR / "reference_data" / "cpsr_reference_long.csv"


def load_reference_data():
    """Load CPSR.m reference data from CSV."""
    data = {}
    
    with open(REFERENCE_CSV, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            msr = int(row["MSR"])
            sof_name = row["SoF_name"].strip()
            # Convert SoF names to match pfdha style names
            style = sof_name.lower().replace('_', '-')
            if style == 'strike-slip':
                style = 'strike-slip'
            
            key = (msr, style)
            if key not in data:
                data[key] = {'mags': [], 'probs': []}
            
            data[key]['mags'].append(float(row["Mw"]))
            data[key]['probs'].append(float(row["CPSR"]))
    
    # Convert to numpy arrays
    for key in data:
        data[key]['mags'] = np.array(data[key]['mags'])
        data[key]['probs'] = np.array(data[key]['probs'])
    
    return data


def generate_plot():
    """Generate the comprehensive comparison plot."""
    
    print("=" * 70)
    print("Mammarella et al. (2024) - Comprehensive Validation Plot")
    print("=" * 70)
    
    # Load reference data
    ref_data = load_reference_data()
    
    # Set publication-quality matplotlib parameters
    rcParams['font.family'] = 'serif'
    rcParams['font.serif'] = ['Times New Roman', 'Times', 'DejaVu Serif']
    rcParams['font.size'] = 12
    rcParams['axes.labelsize'] = 14
    rcParams['axes.titlesize'] = 14
    rcParams['xtick.labelsize'] = 11
    rcParams['ytick.labelsize'] = 11
    rcParams['legend.fontsize'] = 10
    rcParams['figure.dpi'] = 300
    rcParams['savefig.dpi'] = 600
    rcParams['savefig.bbox'] = 'tight'
    
    # Model and reference parameters (from the MATLAB run)
    model = MammarellaEtAl2024PrimarySR()
    
    # Parameters matching the reference data
    HDD = 'AGG_R'  # μ=0.66, σ=0.24
    dip_mu = 40.0
    dip_sigma = 2.0
    t_d = 2.5
    Zs_mu = 14.0
    Zs_sigma = 2.0
    t_z = 1.0
    
    # Magnitude range matching reference
    magnitudes = np.arange(5.0, 8.01, 0.1)
    
    # MSR configurations
    msr_config = {
        0: {'name': 'WC94 (Wells & Coppersmith)', 'short': 'WC94'},
        1: {'name': 'LSRF (Leonard 2010)', 'short': 'LSRF'},
        2: {'name': 'TMG (Thingbaijam 2017)', 'short': 'TMG'},
    }
    
    # Fault styles configuration
    styles_config = {
        'normal': {'color': '#2ca02c', 'label': 'Normal'},
        'reverse': {'color': '#d62728', 'label': 'Reverse'},
        'strike-slip': {'color': '#1f77b4', 'label': 'Strike-Slip'},
    }
    
    # Create figure with 3x3 subplots
    fig, axes = plt.subplots(3, 3, figsize=(14, 12))
    
    print(f"\nReference params: HDD={HDD}, dip_mu={dip_mu}°, dip_sigma={dip_sigma}°")
    print(f"                  Zs_mu={Zs_mu} km, Zs_sigma={Zs_sigma} km")
    print(f"                  t_d={t_d}, t_z={t_z}\n")
    
    # Store max differences for summary table
    max_diffs = {}
    
    for row_idx, (msr, msr_info) in enumerate(msr_config.items()):
        print(f"\nMSR={msr} ({msr_info['short']}):")
        
        for col_idx, (style, config) in enumerate(styles_config.items()):
            ax = axes[row_idx, col_idx]
            
            # Get reference data
            ref_key = (msr, style)
            ref_mags = ref_data[ref_key]['mags']
            ref_probs = ref_data[ref_key]['probs']
            
            # Run pfdha model
            pfdha_probs = []
            for mag in magnitudes:
                p = model.get_prob(
                    mag=mag,
                    MSR=msr,
                    HDD_str=HDD,
                    dip_mu=dip_mu,
                    dip_sigma=dip_sigma,
                    t_d=t_d,
                    Zs_sigma=Zs_sigma,
                    t_z=t_z,
                    style=style,
                    seismothickness=Zs_mu,
                )
                pfdha_probs.append(p)
            pfdha_probs = np.array(pfdha_probs)
            
            # Interpolate reference to match pfdha magnitudes for comparison
            ref_probs_interp = np.interp(magnitudes, ref_mags, ref_probs)
            
            # Calculate max relative difference
            mask = ref_probs_interp > 1e-6
            if np.any(mask):
                rel_diff = np.max(np.abs(pfdha_probs[mask] - ref_probs_interp[mask]) / ref_probs_interp[mask])
            else:
                rel_diff = 0.0
            max_diffs[(msr, style)] = rel_diff
            
            print(f"  {config['label']:12s}: max rel diff = {rel_diff:.2e}")
            
            # Plot reference (CPSR.m) - thick dashed line
            ax.plot(ref_mags, ref_probs,
                    color=config['color'],
                    linestyle='--',
                    linewidth=2.5,
                    label='CPSR.m (reference)')
            
            # Plot pfdha - solid line with markers
            marker_idx = np.linspace(0, len(magnitudes)-1, 8, dtype=int)
            ax.plot(magnitudes, pfdha_probs,
                    color=config['color'],
                    linestyle='-',
                    linewidth=1.5,
                    marker='^',
                    markersize=5,
                    markevery=list(marker_idx),
                    label='pfdha (our work)')
            
            # Format subplot
            if row_idx == 2:
                ax.set_xlabel('Magnitude, Mw', fontweight='bold')
            if col_idx == 0:
                ax.set_ylabel('P(Surface Rupture)', fontweight='bold')
            
            # Title: MSR name for top row, style for first column
            if row_idx == 0:
                ax.set_title(f'{config["label"]}', fontweight='bold', fontsize=13)
            
            # Add MSR label on left
            if col_idx == 0:
                ax.text(-0.25, 0.5, msr_info['short'],
                       transform=ax.transAxes,
                       fontsize=12, fontweight='bold',
                       rotation=90, va='center', ha='center')
            
            ax.set_xlim([5.0, 8.0])
            ax.set_ylim([0, 1.05])
            ax.grid(True, alpha=0.3)
            
            # Add max diff annotation
            ax.text(0.97, 0.03, f'Max diff:\n{rel_diff:.2e}',
                    transform=ax.transAxes, fontsize=9,
                    verticalalignment='bottom',
                    horizontalalignment='right',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                             edgecolor='gray', alpha=0.9))
            
            # Legend only for first subplot
            if row_idx == 0 and col_idx == 0:
                ax.legend(loc='upper left', framealpha=0.95, fontsize=9)
    
    # Add overall title
    fig.suptitle(
        'Mammarella et al. (2024) - pfdha vs CPSR.m Reference Validation\n'
        f'Parameters: HDD={HDD}, $\\delta_{{\\mu}}$={dip_mu}°, $Z_{{s,\\mu}}$={Zs_mu} km',
        fontsize=14, fontweight='bold', y=0.995
    )
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.92, left=0.08)
    
    # Save figure
    OUTPUT_DIR.mkdir(exist_ok=True)
    
    output_path = OUTPUT_DIR / "mammarella2024_comprehensive_validation.png"
    plt.savefig(output_path, dpi=600, bbox_inches='tight', facecolor='white')
    print(f"\nPlot saved to: {output_path}")
    
    pdf_path = OUTPUT_DIR / "mammarella2024_comprehensive_validation.pdf"
    plt.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    print(f"PDF saved to: {pdf_path}")
    
    plt.close(fig)
    
    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY - Max Relative Differences:")
    print("=" * 70)
    for (msr, style), diff in max_diffs.items():
        print(f"  MSR={msr} ({msr_config[msr]['short']:4s}), {style:12s}: {diff:.2e}")
    
    print("\nDone!")
    return 0


if __name__ == "__main__":
    sys.exit(generate_plot())




