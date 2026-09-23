#!/usr/bin/env python3
"""
Mammarella et al. (2024) validation - Grouped by MSR.
Each subplot shows all 3 fault styles for one MSR.
"""

import csv
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path

from openquake.pfd.primary_surf_rup.mammarella2024 import MammarellaEtAl2024PrimarySR

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / "outputs"
# Reference data (place in reference_data subdirectory or update path)
REFERENCE_CSV = SCRIPT_DIR / "reference_data" / "cpsr_reference_long.csv"


def load_reference_data():
    data = {}
    with open(REFERENCE_CSV, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            msr = int(row["MSR"])
            style = row["SoF_name"].strip().lower().replace('_', '-')
            key = (msr, style)
            if key not in data:
                data[key] = {'mags': [], 'probs': []}
            data[key]['mags'].append(float(row["Mw"]))
            data[key]['probs'].append(float(row["CPSR"]))
    for key in data:
        data[key]['mags'] = np.array(data[key]['mags'])
        data[key]['probs'] = np.array(data[key]['probs'])
    return data


def generate_plot():
    print("=" * 70)
    print("Mammarella et al. (2024) - Grouped by MSR")
    print("=" * 70)
    
    ref_data = load_reference_data()
    model = MammarellaEtAl2024PrimarySR()
    
    rcParams['font.family'] = 'serif'
    rcParams['font.size'] = 12
    rcParams['axes.labelsize'] = 14
    rcParams['axes.titlesize'] = 14
    rcParams['legend.fontsize'] = 11
    
    # Parameters matching reference
    HDD, dip_mu, dip_sigma, t_d = 'AGG_R', 40.0, 2.0, 2.5
    Zs_mu, Zs_sigma, t_z = 14.0, 2.0, 1.0
    magnitudes = np.arange(5.0, 8.01, 0.1)
    
    msr_config = {
        0: 'WC94 (Wells & Coppersmith 1994)',
        1: 'LSRF (Leonard 2010)',
        2: 'TMG (Thingbaijam 2017)',
    }
    
    styles_config = {
        'normal': {'color': '#2ca02c', 'label': 'Normal'},
        'reverse': {'color': '#d62728', 'label': 'Reverse'},
        'strike-slip': {'color': '#1f77b4', 'label': 'Strike-Slip'},
    }
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    for idx, (msr, msr_name) in enumerate(msr_config.items()):
        ax = axes[idx]
        print(f"\nMSR={msr} ({msr_name}):")
        
        for style, config in styles_config.items():
            ref_key = (msr, style)
            ref_mags = ref_data[ref_key]['mags']
            ref_probs = ref_data[ref_key]['probs']
            
            pfdha_probs = []
            for mag in magnitudes:
                p = model.get_prob(mag=mag, MSR=msr, HDD_str=HDD,
                                   dip_mu=dip_mu, dip_sigma=dip_sigma, t_d=t_d,
                                   Zs_sigma=Zs_sigma, t_z=t_z, style=style,
                                   seismothickness=Zs_mu)
                pfdha_probs.append(p)
            pfdha_probs = np.array(pfdha_probs)
            
            # Reference - dashed
            ax.plot(ref_mags, ref_probs, color=config['color'], linestyle='--',
                    linewidth=2.5, label=f'{config["label"]} (ref)')
            
            # pfdha - solid with markers
            marker_idx = np.linspace(0, len(magnitudes)-1, 8, dtype=int)
            ax.plot(magnitudes, pfdha_probs, color=config['color'], linestyle='-',
                    linewidth=1.5, marker='^', markersize=5, markevery=list(marker_idx),
                    label=f'{config["label"]} (pfdha)')
            
            # Max diff
            ref_interp = np.interp(magnitudes, ref_mags, ref_probs)
            mask = ref_interp > 1e-6
            rel_diff = np.max(np.abs(pfdha_probs[mask] - ref_interp[mask]) / ref_interp[mask]) if np.any(mask) else 0
            print(f"  {config['label']:12s}: max rel diff = {rel_diff:.2e}")
        
        ax.set_xlabel('Magnitude, Mw', fontweight='bold')
        if idx == 0:
            ax.set_ylabel('P(Surface Rupture)', fontweight='bold')
        ax.set_title(msr_name, fontweight='bold')
        ax.set_xlim([5.0, 8.0])
        ax.set_ylim([0, 1.05])
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper left', fontsize=9, ncol=2)
    
    fig.suptitle('Mammarella et al. (2024) Validation - Grouped by MSR\n'
                 'Dashed: CPSR.m reference | Solid+△: pfdha',
                 fontsize=13, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    OUTPUT_DIR.mkdir(exist_ok=True)
    
    output_path = OUTPUT_DIR / "mammarella2024_by_msr.png"
    plt.savefig(output_path, dpi=600, bbox_inches='tight', facecolor='white')
    print(f"\nSaved: {output_path}")
    
    plt.savefig(OUTPUT_DIR / "mammarella2024_by_msr.pdf", bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(generate_plot())




