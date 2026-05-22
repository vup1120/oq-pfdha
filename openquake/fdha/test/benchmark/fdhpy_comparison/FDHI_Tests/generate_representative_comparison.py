#!/usr/bin/env python3
"""
Generate Representative FDHI Comparison Plot

This script creates a publication-quality figure comparing hazard curves from
6 primary fault displacement models, showing both pfdha (OpenQuake implementation)
and fdhpy (FDHI reference implementation) results to demonstrate good agreement.

Models included:
1. Youngs et al. (2003) - D/AD version
2. Petersen et al. (2011) - elliptical version
3. Moss et al. (2024) - D/AD, GIRS version
4. Kuehn et al. (2024) - normal style
5. Lavrentiadis & Abrahamson (2023) - aggregate metric
6. Chiou et al. (2025) - model7 version

Usage:
    cd openquake/fdha/test/FDHI_Tests
    python generate_representative_comparison.py
"""

import sys
from pathlib import Path

import numpy as np

# ============================================================================
# Verify dependencies
# ============================================================================

try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
except ImportError:
    print("[ERROR] matplotlib is required: pip install matplotlib")
    sys.exit(1)

try:
    from fdhpy import (
        YoungsEtAl2003,
        PetersenEtAl2011,
        MossEtAl2024,
        KuehnEtAl2024,
        LavrentiadisAbrahamson2023,
        ChiouEtAl2025,
    )
except ImportError as e:
    print(f"[ERROR] fdhpy is not installed: {e}")
    print("Install with: pip install fdhpy")
    sys.exit(1)

try:
    from openquake.fdha.primary_surf_displ import (
        Youngs2003PrimaryFD,
        Petersen2011PrimaryFD,
        Moss2024PrimaryFD,
        Kuehn2024PrimaryFD,
        Lavrentiadis2023PrimaryFD,
        Chiou2025PrimaryFD,
    )
except ImportError as e:
    print(f"[ERROR] pfdha models not importable: {e}")
    print("Install with: pip install -e /path/to/pfdha")
    sys.exit(1)


# ============================================================================
# Configuration
# ============================================================================

# Use consistent parameters for all models
MAG = 7.0
XL = 0.5
DISPLACEMENTS = np.logspace(-3, 1, 50)  # 0.001 to 10 m

# Define colors for each model (colorblind-friendly palette)
COLORS = {
    'Youngs2003': '#1f77b4',       # blue
    'Petersen2011': '#ff7f0e',     # orange  
    'Moss2024': '#2ca02c',         # green
    'Kuehn2024': '#d62728',        # red
    'Lavrentiadis2023': '#9467bd', # purple
    'Chiou2025': '#8c564b',        # brown
}

# Model display names for legend
MODEL_NAMES = {
    'Youngs2003': 'Youngs et al. (2003)',
    'Petersen2011': 'Petersen et al. (2011)',
    'Moss2024': 'Moss et al. (2024)',
    'Kuehn2024': 'Kuehn et al. (2024)',
    'Lavrentiadis2023': 'Lavrentiadis & Abrahamson (2023)',
    'Chiou2025': 'Chiou et al. (2025)',
}


# ============================================================================
# Model Runners
# ============================================================================

def run_youngs2003(displacements: np.ndarray, mag: float, xl: float):
    """Run Youngs 2003 for both implementations."""
    # fdhpy
    fdhpy_model = YoungsEtAl2003(
        magnitude=mag,
        xl=xl,
        version="d/ad",
        displ_array=displacements,
    )
    fdhpy_probs = fdhpy_model.prob_exceed
    
    # pfdha
    pfdha_model = Youngs2003PrimaryFD()
    pfdha_probs = pfdha_model.get_prob(
        d=displacements,
        X_L_ratio=np.array([xl]),
        mag=mag,
        style="all",
        norm_disp_type="AD"
    ).flatten()
    
    return fdhpy_probs, pfdha_probs


def run_petersen2011(displacements: np.ndarray, mag: float, xl: float):
    """Run Petersen 2011 for both implementations."""
    # fdhpy
    fdhpy_model = PetersenEtAl2011(
        magnitude=mag,
        xl=xl,
        version="elliptical",
        displ_array=displacements,
    )
    fdhpy_probs = fdhpy_model.prob_exceed
    
    # pfdha
    pfdha_model = Petersen2011PrimaryFD()
    pfdha_probs = pfdha_model.get_prob(
        d=displacements,
        X_L_ratio=np.array([xl]),
        mag=mag,
        version="elliptical"
    ).flatten()
    
    return fdhpy_probs, pfdha_probs


def run_moss2024(displacements: np.ndarray, mag: float, xl: float):
    """Run Moss 2024 for both implementations."""
    # fdhpy
    fdhpy_model = MossEtAl2024(
        magnitude=mag,
        xl=xl,
        version="d/ad",
        displ_array=displacements,
        use_girs=True,
        complete=True,
    )
    fdhpy_probs = fdhpy_model.prob_exceed
    
    # pfdha
    pfdha_model = Moss2024PrimaryFD()
    pfdha_probs = pfdha_model.get_prob(
        d=displacements,
        X_L_ratio=np.array([xl]),
        mag=mag,
        version="AD",
        source="GIRS",
        completeness="complete"
    ).flatten()
    
    return fdhpy_probs, pfdha_probs


def run_kuehn2024(displacements: np.ndarray, mag: float, xl: float):
    """Run Kuehn 2024 for both implementations."""
    # fdhpy
    fdhpy_model = KuehnEtAl2024(
        style="normal",
        magnitude=mag,
        xl=xl,
        version="median_coeffs",
        folded=True,
        displ_array=displacements,
    )
    fdhpy_probs = fdhpy_model.prob_exceed
    
    # pfdha
    pfdha_model = Kuehn2024PrimaryFD()
    pfdha_probs = pfdha_model.get_prob(
        d=displacements,
        X_L_ratio=np.array([xl]),
        mag=mag,
        style="normal",
        folded=True,
        epistemic_uncertainty=False
    ).flatten()
    
    return fdhpy_probs, pfdha_probs


def run_lavrentiadis2023(displacements: np.ndarray, mag: float, xl: float):
    """Run Lavrentiadis 2023 for both implementations."""
    # fdhpy
    fdhpy_model = LavrentiadisAbrahamson2023(
        magnitude=mag,
        xl=xl,
        displ_array=displacements,
        metric="aggregate",
        version="full rupture",
        style="strike-slip",
        include_prob_zero=True,
    )
    fdhpy_probs = fdhpy_model.prob_exceed
    
    # pfdha
    pfdha_model = Lavrentiadis2023PrimaryFD()
    pfdha_probs = pfdha_model.get_prob(
        d=displacements,
        X_L_ratio=np.array([xl]),
        mag=mag,
        style="strike-slip",
        output_type="disp_agg_prime",
        include_zero_slip=True,
    ).flatten()
    
    return fdhpy_probs, pfdha_probs


def run_chiou2025(displacements: np.ndarray, mag: float, xl: float):
    """Run Chiou 2025 for both implementations."""
    # fdhpy
    fdhpy_model = ChiouEtAl2025(
        magnitude=mag,
        xl=xl,
        version="model7",
        displ_array=displacements,
    )
    fdhpy_probs = fdhpy_model.prob_exceed
    
    # pfdha
    pfdha_model = Chiou2025PrimaryFD()
    pfdha_probs = pfdha_model.get_prob(
        d=displacements,
        X_L_ratio=np.array([xl]),
        mag=mag,
        version="model7"
    ).flatten()
    
    return fdhpy_probs, pfdha_probs


# ============================================================================
# Main Plot Generation
# ============================================================================

def generate_comparison_plot():
    """Generate the representative comparison plot with summary table."""
    
    print("=" * 70)
    print("FDHI Reference Comparison Plot Generator")
    print("=" * 70)
    print(f"\nParameters: M = {MAG}, x/L = {XL}")
    print(f"Displacement range: {DISPLACEMENTS.min():.3f} to {DISPLACEMENTS.max():.1f} m")
    
    # Set publication-quality matplotlib parameters
    rcParams['font.family'] = 'serif'
    rcParams['font.serif'] = ['Times New Roman', 'Times', 'DejaVu Serif']
    rcParams['font.size'] = 14
    rcParams['axes.labelsize'] = 16
    rcParams['axes.titlesize'] = 18
    rcParams['xtick.labelsize'] = 13
    rcParams['ytick.labelsize'] = 13
    rcParams['legend.fontsize'] = 14
    rcParams['figure.dpi'] = 600
    rcParams['savefig.dpi'] = 600
    rcParams['savefig.bbox'] = 'tight'
    
    # Create figure with single plot (table inside)
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Define model runners
    model_runners = {
        'Youngs2003': run_youngs2003,
        'Petersen2011': run_petersen2011,
        'Moss2024': run_moss2024,
        'Kuehn2024': run_kuehn2024,
        'Lavrentiadis2023': run_lavrentiadis2023,
        'Chiou2025': run_chiou2025,
    }
    
    print("\nRunning models...")
    
    # Store results for table
    table_data = []
    
    # Run each model and plot
    for model_key, runner in model_runners.items():
        print(f"  {MODEL_NAMES[model_key]}...", end=" ", flush=True)
        
        try:
            fdhpy_probs, pfdha_probs = runner(DISPLACEMENTS, MAG, XL)
            
            color = COLORS[model_key]
            display_name = MODEL_NAMES[model_key]
            
            # Plot reference (fdhpy) - thick dashed line
            ax.loglog(
                DISPLACEMENTS, fdhpy_probs,
                color=color,
                linestyle='--',
                linewidth=2.5,
                label=f'{display_name} (reference)',
                alpha=0.9,
            )
            
            # Plot implementation (pfdha) - solid line with triangle markers
            # Use fewer markers for clarity
            marker_indices = np.linspace(0, len(DISPLACEMENTS)-1, 8, dtype=int)
            ax.loglog(
                DISPLACEMENTS, pfdha_probs,
                color=color,
                linestyle='-',
                linewidth=1.5,
                marker='^',
                markersize=6,
                markevery=list(marker_indices),
                label=f'{display_name} (our work)',
                alpha=0.9,
            )
            
            # Calculate max relative difference
            mask = fdhpy_probs > 1e-10
            if np.any(mask):
                rel_diff = np.max(np.abs(pfdha_probs[mask] - fdhpy_probs[mask]) / fdhpy_probs[mask])
                print(f"max rel diff = {rel_diff:.2e}")
            else:
                rel_diff = 0.0
                print("OK")
            
            # Store for table
            table_data.append({
                'model': display_name,
                'color': color,
                'rel_diff': rel_diff,
            })
                
        except Exception as e:
            print(f"ERROR: {e}")
            continue
    
    # Format plot
    ax.set_xlabel('Displacement, d (m)', fontweight='bold')
    ax.set_ylabel('Exceedance Probability, P(D > d)', fontweight='bold')
    ax.set_title(
        f'FDHI Reference Suite Comparison\n'
        f'M = {MAG}, x/L = {XL}',
        fontweight='bold'
    )
    
    # Set axis limits
    ax.set_xlim([0.05, 11])
    ax.set_ylim([1e-4, 1.1])
    
    # Grid
    ax.grid(True, which='major', alpha=0.4, linestyle='-')
    ax.grid(True, which='minor', alpha=0.2, linestyle=':')
    
    # Legend - organize in two columns
    handles, labels = ax.get_legend_handles_labels()
    
    # Create custom legend with reference and implementation together
    ax.legend(
        handles, labels,
        loc='lower left',
        ncol=2,
        framealpha=0.95,
        edgecolor='gray',
        fancybox=True,
        fontsize=9,
    )
    
    
    # Create summary table inside the plot (upper right)
    col_labels = ['Model', 'Max Rel. Diff.']
    
    # Table data
    cell_text = []
    for item in table_data:
        # Format relative difference
        if item['rel_diff'] < 1e-10:
            diff_str = '< 1e-10'
        else:
            diff_str = f"{item['rel_diff']:.2e}"
        cell_text.append([item['model'], diff_str])
    
    # Create table inside plot (left middle)
    table = ax.table(
        cellText=cell_text,
        colLabels=col_labels,
        cellLoc='left',
        bbox=[0.02, 0.2, 0.46, 0.45],  # [left, bottom, width, height] - lower left
        colWidths=[0.65, 0.35],
    )
    
    # Style the table
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.0, 1.6)
    
    # Style header
    for j, label in enumerate(col_labels):
        cell = table[(0, j)]
        cell.set_text_props(fontweight='bold')
        cell.set_facecolor('#e0e0e0')
    
    # Style data rows with model colors
    for i, item in enumerate(table_data):
        # Model name cell - use light model color
        table[(i + 1, 0)].set_facecolor(item['color'] + '20')
        # Difference cell - white background
        table[(i + 1, 1)].set_facecolor('white')
    
    # Save figure
    script_dir = Path(__file__).parent
    output_dir = script_dir / "outputs"
    output_dir.mkdir(exist_ok=True)
    
    output_path = output_dir / "representative_fdhi_comparison.png"
    plt.savefig(output_path, dpi=1200, bbox_inches='tight', facecolor='white')
    print(f"\nPlot saved to: {output_path}")
    
    # Also save PDF version
    pdf_path = output_dir / "representative_fdhi_comparison.pdf"
    plt.savefig(pdf_path, bbox_inches='tight', facecolor='white')
    print(f"PDF saved to: {pdf_path}")
    
    plt.close(fig)
    
    print("\nDone!")
    return 0


if __name__ == "__main__":
    sys.exit(generate_comparison_plot())








