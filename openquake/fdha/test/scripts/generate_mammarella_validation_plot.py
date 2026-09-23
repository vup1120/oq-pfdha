#!/usr/bin/env python3
"""
Generate validation plot for Mammarella et al. (2024) PrimarySR model.

Compares pfdha implementation against CPSR.m golden reference values
across three fault styles: reverse, normal, and strike-slip.
"""

import csv
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path

from openquake.pfd.primary_surf_rup.mammarella2024 import MammarellaEtAl2024PrimarySR


# Paths
SCRIPT_DIR = Path(__file__).parent
GOLDEN_PATH = SCRIPT_DIR / "expected" / "mammarella2024_cpsr.csv"
OUTPUT_DIR = SCRIPT_DIR / "outputs"


def load_golden_data():
    """Load golden reference data from CSV."""
    data = {'reverse': [], 'normal': [], 'strike-slip': []}
    
    with open(GOLDEN_PATH, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            style = row["style"].strip()
            entry = {
                'Mw': float(row["Mw"]),
                'MSR': int(row["MSR"]),
                'HDD': row["HDD"].strip(),
                'Zs_mu': float(row["Zs_mu"]),
                'dip_mu': float(row["dip_mu"]),
                'dip_sigma': float(row["dip_sigma"]),
                't_d': float(row["t_d"]),
                'Zs_sigma': float(row["Zs_sigma"]),
                't_z': float(row["t_z"]),
                'expected': float(row["expected"]),
            }
            data[style].append(entry)
    
    return data


def run_pfdha_model(entries):
    """Run pfdha model for given entries and return probabilities."""
    model = MammarellaEtAl2024PrimarySR()
    probs = []
    
    for entry in entries:
        p = model.get_prob(
            mag=entry['Mw'],
            MSR=entry['MSR'],
            HDD_str=entry['HDD'],
            dip_mu=entry['dip_mu'],
            dip_sigma=entry['dip_sigma'],
            t_d=entry['t_d'],
            Zs_sigma=entry['Zs_sigma'],
            t_z=entry['t_z'],
            style=entries[0]['Mw'],  # This should be style, fix below
            seismothickness=entry['Zs_mu'],
        )
        probs.append(p)
    
    return np.array(probs)


def generate_plot():
    """Generate the validation comparison plot."""
    
    print("=" * 60)
    print("Mammarella et al. (2024) Validation Plot Generator")
    print("=" * 60)
    
    # Load data
    data = load_golden_data()
    
    # Set publication-quality matplotlib parameters
    rcParams['font.family'] = 'serif'
    rcParams['font.serif'] = ['Times New Roman', 'Times', 'DejaVu Serif']
    rcParams['font.size'] = 12
    rcParams['axes.labelsize'] = 14
    rcParams['axes.titlesize'] = 14
    rcParams['xtick.labelsize'] = 12
    rcParams['ytick.labelsize'] = 12
    rcParams['legend.fontsize'] = 11
    rcParams['figure.dpi'] = 300
    rcParams['savefig.dpi'] = 600
    rcParams['savefig.bbox'] = 'tight'
    
    # Create figure with 1x3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    
    styles = ['reverse', 'normal', 'strike-slip']
    style_titles = ['Reverse Faulting', 'Normal Faulting', 'Strike-Slip Faulting']
    colors = ['#d62728', '#2ca02c', '#1f77b4']
    
    model = MammarellaEtAl2024PrimarySR()
    max_diffs = []
    
    for idx, (style, title, color) in enumerate(zip(styles, style_titles, colors)):
        ax = axes[idx]
        entries = data[style]
        
        # Sort by magnitude
        entries = sorted(entries, key=lambda x: x['Mw'])
        
        # Get magnitudes and reference probabilities
        mags = np.array([e['Mw'] for e in entries])
        ref_probs = np.array([e['expected'] for e in entries])
        
        # Run pfdha model
        pfdha_probs = []
        for entry in entries:
            p = model.get_prob(
                mag=entry['Mw'],
                MSR=entry['MSR'],
                HDD_str=entry['HDD'],
                dip_mu=entry['dip_mu'],
                dip_sigma=entry['dip_sigma'],
                t_d=entry['t_d'],
                Zs_sigma=entry['Zs_sigma'],
                t_z=entry['t_z'],
                style=style,
                seismothickness=entry['Zs_mu'],
            )
            pfdha_probs.append(p)
        pfdha_probs = np.array(pfdha_probs)
        
        # Calculate max relative difference
        mask = ref_probs > 1e-10
        if np.any(mask):
            rel_diff = np.max(np.abs(pfdha_probs[mask] - ref_probs[mask]) / ref_probs[mask])
        else:
            rel_diff = 0.0
        max_diffs.append(rel_diff)
        
        print(f"\n{title}:")
        print(f"  Mw: {mags}")
        print(f"  Reference: {ref_probs}")
        print(f"  pfdha:     {pfdha_probs}")
        print(f"  Max rel diff: {rel_diff:.2e}")
        
        # Plot reference (CPSR.m) - thick dashed line with circle markers
        ax.plot(mags, ref_probs, 
                color=color, linestyle='--', linewidth=2.5,
                marker='o', markersize=10,
                label='CPSR.m (reference)')
        
        # Plot pfdha - solid line with triangle markers
        ax.plot(mags, pfdha_probs,
                color=color, linestyle='-', linewidth=1.5,
                marker='^', markersize=8,
                label='pfdha (our work)')
        
        # Format subplot
        ax.set_xlabel('Magnitude, Mw', fontweight='bold')
        if idx == 0:
            ax.set_ylabel('Surface Rupture Probability', fontweight='bold')
        ax.set_title(title, fontweight='bold')
        ax.set_xlim([6.3, 7.7])
        ax.set_ylim([0, 1.05])
        ax.set_xticks([6.5, 7.0, 7.5])
        ax.grid(True, alpha=0.3)
        ax.legend(loc='lower right')
        
        # Add max diff annotation
        ax.text(0.05, 0.95, f'Max diff: {rel_diff:.2e}',
                transform=ax.transAxes, fontsize=10,
                verticalalignment='top',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                         edgecolor='gray', alpha=0.9))
    
    # Add overall title with common parameters
    fig.suptitle(
        'Mammarella et al. (2024) Primary Surface Rupture Model Validation\n'
        r'Common parameters: $\delta_{mu}$=45°, $\delta_{\sigma}$=10°, $Z_{s,mu}$=12 km, $Z_{s,\sigma}$=2 km',
        fontsize=13, fontweight='bold', y=1.02
    )
    
    plt.tight_layout()
    
    # Save figure
    OUTPUT_DIR.mkdir(exist_ok=True)
    
    output_path = OUTPUT_DIR / "mammarella2024_validation.png"
    plt.savefig(output_path, dpi=600, bbox_inches='tight', facecolor='white')
    print(f"\nPlot saved to: {output_path}")
    
    pdf_path = OUTPUT_DIR / "mammarella2024_validation.pdf"
    plt.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    print(f"PDF saved to: {pdf_path}")
    
    plt.close(fig)
    
    print("\nDone!")
    return 0


if __name__ == "__main__":
    sys.exit(generate_plot())
