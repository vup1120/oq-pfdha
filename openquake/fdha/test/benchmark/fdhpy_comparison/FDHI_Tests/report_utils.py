"""
PFDHA Comparison Report Utilities

This module provides reusable utilities for generating comparison reports,
plots, and summary tables for validating pfdha against fdhpy.

Usage:
    These utilities are used by generate_comparison_dashboard.py to create
    visual and tabular comparison outputs.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional, Tuple
from datetime import datetime
import json
import warnings


# ============================================================================
# Configuration
# ============================================================================

# Standard displacement grid (same as tests)
DISPLACEMENTS = np.array([0.001, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0])

# Fine displacement grid for plotting
FINE_DISPLACEMENTS = np.logspace(-3, 1.5, 100)

# Default tolerances (same as tests)
DEFAULT_TOLERANCES = {
    "youngs2003": {"rtol": 1e-4, "atol": 1e-8},
    "petersen2011": {"rtol": 1e-6, "atol": 1e-10},
    "moss2024": {"rtol": 1e-6, "atol": 1e-10},
    "kuehn2024": {"rtol": 1e-6, "atol": 1e-10},
    "lavrentiadis2023": {"rtol": 1e-6, "atol": 1e-10},
    "chiou2025": {"rtol": 1e-6, "atol": 1e-10},
}


# ============================================================================
# Parameter Grids (shared with tests)
# ============================================================================

def get_parameter_grids() -> Dict[str, Dict[str, Any]]:
    """
    Return parameter grids for all models.
    These match the parameters used in test_simple_reference.py.
    """
    return {
        "youngs2003": {
            "magnitudes": [6.0, 6.5, 7.0, 7.5],
            "x_l_ratios": [0.1, 0.25, 0.5],
            "version": "d/ad",
            "style": "all",
        },
        "petersen2011": {
            "magnitudes": [6.5, 7.0, 7.5],
            "x_l_ratios": [0.1, 0.3, 0.5],
            "versions": ["elliptical", "quadratic"],
        },
        "moss2024": {
            "magnitudes": [6.5, 7.0, 7.5],
            "x_l_ratios": [0.25, 0.5],
            "versions": [("d/ad", "AD"), ("d/md", "MD")],
            "use_girs": [True, False],
        },
        "kuehn2024": {
            "magnitudes": [7.0, 7.2, 7.6],
            "x_l_ratios": [0.3, 0.5, 0.7],
            "styles": ["normal", "reverse", "strike-slip"],
        },
        "lavrentiadis2023": {
            "magnitudes": [6.5, 7.0, 7.5],
            "x_l_ratios": [0.3, 0.5, 0.6],
            "include_prob_zero": [True, False],
        },
        "chiou2025": {
            "magnitudes": [6.5, 7.0, 7.5],
            "x_l_ratios": [0.25, 0.5],
            "versions": ["model7", "model8.2"],
        },
    }


# ============================================================================
# Data Classes
# ============================================================================

@dataclass
class ComparisonMetrics:
    """Metrics from comparing two probability arrays."""
    max_abs_diff: float
    mean_abs_diff: float
    max_rel_diff: float
    mean_rel_diff: float
    rmse: float
    n_values: int
    within_tolerance: bool = True
    
    @classmethod
    def compute(
        cls,
        expected: np.ndarray,
        actual: np.ndarray,
        rtol: float = 1e-6,
        atol: float = 1e-10,
    ) -> "ComparisonMetrics":
        """Compute comparison metrics between two arrays."""
        diff = np.abs(actual - expected)
        max_abs_diff = float(np.max(diff))
        mean_abs_diff = float(np.mean(diff))
        
        # Relative difference (handle zeros safely)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            rel_diff = np.abs(diff / (expected + 1e-15))
            # Mask out where expected is very small
            mask = np.abs(expected) > 1e-12
            if np.any(mask):
                max_rel_diff = float(np.max(rel_diff[mask]))
                mean_rel_diff = float(np.mean(rel_diff[mask]))
            else:
                max_rel_diff = 0.0
                mean_rel_diff = 0.0
        
        rmse = float(np.sqrt(np.mean(diff ** 2)))
        
        # Check tolerance
        within_tolerance = np.allclose(expected, actual, rtol=rtol, atol=atol)
        
        return cls(
            max_abs_diff=max_abs_diff,
            mean_abs_diff=mean_abs_diff,
            max_rel_diff=max_rel_diff,
            mean_rel_diff=mean_rel_diff,
            rmse=rmse,
            n_values=len(expected),
            within_tolerance=within_tolerance,
        )


@dataclass
class TestCaseResult:
    """Result from a single test case comparison."""
    model_name: str
    variant: str
    magnitude: float
    x_l_ratio: float
    extra_params: Dict[str, Any]
    metrics: ComparisonMetrics
    fdhpy_probs: np.ndarray
    pfdha_probs: np.ndarray
    displacements: np.ndarray
    tolerance_rtol: float
    tolerance_atol: float


# ============================================================================
# Plotting Utilities
# ============================================================================

def create_hazard_curve_plot(
    result: TestCaseResult,
    output_path: Path,
    show_residuals: bool = True,
    figsize: Tuple[float, float] = (10, 8),
) -> None:
    """
    Create a hazard curve comparison plot.
    
    Args:
        result: TestCaseResult with comparison data
        output_path: Path to save the PNG file
        show_residuals: If True, add a residual subplot
        figsize: Figure size in inches
    """
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    
    if show_residuals:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, height_ratios=[3, 1], sharex=True)
    else:
        fig, ax1 = plt.subplots(1, 1, figsize=figsize)
        ax2 = None
    
    d = result.displacements
    fdhpy = result.fdhpy_probs
    pfdha = result.pfdha_probs
    
    # Main hazard curve plot
    ax1.loglog(d, fdhpy, 'b-', linewidth=2, label='fdhpy (reference)', marker='o', markersize=4)
    ax1.loglog(d, pfdha, 'r--', linewidth=2, label='pfdha (OpenQuake)', marker='s', markersize=4)
    
    ax1.set_ylabel('P[D > d]', fontsize=12)
    ax1.set_title(
        f'{result.model_name} ({result.variant})\n'
        f'M = {result.magnitude}, x/L = {result.x_l_ratio}',
        fontsize=14
    )
    ax1.legend(loc='upper right', fontsize=10)
    ax1.grid(True, which='both', alpha=0.3)
    ax1.set_ylim([1e-6, 1.1])
    
    # Add metrics annotation
    metrics_text = (
        f"Max abs diff: {result.metrics.max_abs_diff:.2e}\n"
        f"Max rel diff: {result.metrics.max_rel_diff:.2e}\n"
        f"Status: {'PASS ✓' if result.metrics.within_tolerance else 'FAIL ✗'}"
    )
    ax1.text(
        0.02, 0.02, metrics_text,
        transform=ax1.transAxes,
        fontsize=9,
        verticalalignment='bottom',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
    )
    
    # Residual plot
    if ax2 is not None:
        residuals = pfdha - fdhpy
        ax2.semilogx(d, residuals, 'g-', linewidth=1.5, marker='d', markersize=3)
        ax2.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
        ax2.set_xlabel('Displacement d (m)', fontsize=12)
        ax2.set_ylabel('Residual\n(pfdha - fdhpy)', fontsize=10)
        ax2.grid(True, which='both', alpha=0.3)
        
        # Add tolerance bands
        max_expected = np.max(fdhpy)
        tol_band = result.tolerance_rtol * max_expected + result.tolerance_atol
        ax2.axhline(y=tol_band, color='orange', linestyle='--', alpha=0.7, label=f'±tolerance')
        ax2.axhline(y=-tol_band, color='orange', linestyle='--', alpha=0.7)
    else:
        ax1.set_xlabel('Displacement d (m)', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)


def create_summary_plot(
    df: pd.DataFrame,
    output_path: Path,
    figsize: Tuple[float, float] = (12, 8),
) -> None:
    """
    Create a summary bar chart of max relative differences by model.
    """
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')
    
    # Aggregate by model
    summary = df.groupby('model_name').agg({
        'max_rel_diff': ['max', 'mean'],
        'within_tolerance': 'mean',
    }).reset_index()
    summary.columns = ['model_name', 'max_rel_diff_max', 'max_rel_diff_mean', 'pass_rate']
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Bar chart of max relative differences
    x = range(len(summary))
    width = 0.35
    ax1.bar([i - width/2 for i in x], summary['max_rel_diff_max'], width, label='Max', color='coral')
    ax1.bar([i + width/2 for i in x], summary['max_rel_diff_mean'], width, label='Mean', color='steelblue')
    ax1.set_ylabel('Relative Difference')
    ax1.set_title('Max Relative Difference by Model')
    ax1.set_xticks(x)
    ax1.set_xticklabels(summary['model_name'], rotation=45, ha='right')
    ax1.legend()
    ax1.set_yscale('log')
    ax1.grid(True, alpha=0.3)
    
    # Pass rate by model
    colors = ['green' if r == 1.0 else 'orange' if r > 0.9 else 'red' for r in summary['pass_rate']]
    ax2.bar(x, summary['pass_rate'] * 100, color=colors)
    ax2.set_ylabel('Pass Rate (%)')
    ax2.set_title('Test Pass Rate by Model')
    ax2.set_xticks(x)
    ax2.set_xticklabels(summary['model_name'], rotation=45, ha='right')
    ax2.set_ylim([0, 105])
    ax2.axhline(y=100, color='green', linestyle='--', alpha=0.5)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)


# ============================================================================
# Summary Generation
# ============================================================================

def results_to_dataframe(results: List[TestCaseResult]) -> pd.DataFrame:
    """Convert list of TestCaseResult to a pandas DataFrame."""
    records = []
    for r in results:
        record = {
            'model_name': r.model_name,
            'variant': r.variant,
            'magnitude': r.magnitude,
            'x_l_ratio': r.x_l_ratio,
            'max_abs_diff': r.metrics.max_abs_diff,
            'mean_abs_diff': r.metrics.mean_abs_diff,
            'max_rel_diff': r.metrics.max_rel_diff,
            'mean_rel_diff': r.metrics.mean_rel_diff,
            'rmse': r.metrics.rmse,
            'within_tolerance': r.metrics.within_tolerance,
            'tolerance_rtol': r.tolerance_rtol,
            'tolerance_atol': r.tolerance_atol,
            'n_values': r.metrics.n_values,
        }
        # Add extra params
        for k, v in r.extra_params.items():
            record[k] = v
        records.append(record)
    
    return pd.DataFrame(records)


def generate_model_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Generate per-model summary statistics."""
    summary = df.groupby(['model_name', 'variant']).agg({
        'max_abs_diff': ['max', 'mean'],
        'max_rel_diff': ['max', 'mean'],
        'within_tolerance': ['sum', 'count'],
    }).reset_index()
    
    # Flatten column names
    summary.columns = [
        'model_name', 'variant',
        'max_abs_diff_max', 'max_abs_diff_mean',
        'max_rel_diff_max', 'max_rel_diff_mean',
        'tests_passed', 'total_tests',
    ]
    
    summary['pass_rate'] = summary['tests_passed'] / summary['total_tests']
    summary['status'] = summary['pass_rate'].apply(
        lambda x: '✅ PASS' if x == 1.0 else ('⚠️ PARTIAL' if x > 0.5 else '❌ FAIL')
    )
    
    return summary


def summary_to_markdown(summary: pd.DataFrame) -> str:
    """Convert summary DataFrame to Markdown table."""
    lines = [
        "| Model | Variant | Tests | Max Rel Diff | Mean Rel Diff | Pass Rate | Status |",
        "|-------|---------|-------|--------------|---------------|-----------|--------|",
    ]
    
    for _, row in summary.iterrows():
        lines.append(
            f"| {row['model_name']} | {row['variant']} | "
            f"{int(row['total_tests'])} | "
            f"{row['max_rel_diff_max']:.2e} | "
            f"{row['max_rel_diff_mean']:.2e} | "
            f"{row['pass_rate']*100:.0f}% | "
            f"{row['status']} |"
        )
    
    return "\n".join(lines)


def generate_markdown_report(
    summary: pd.DataFrame,
    df: pd.DataFrame,
    plots_dir: Path,
    representative_plots: Dict[str, List[str]],
) -> str:
    """Generate a complete Markdown comparison report."""
    
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    
    lines = [
        "# PFDHA vs fdhpy Comparison Dashboard",
        "",
        f"**Generated:** {now}",
        "",
        "## Overview",
        "",
        "This report presents a comprehensive comparison between the `pfdha` (OpenQuake) ",
        "implementation and the `fdhpy` (FDHI) reference implementation for Probabilistic ",
        "Fault Displacement Hazard Analysis models.",
        "",
        "The validation covers multiple magnitudes, normalized positions (x/L), and model variants. ",
        "For each test case, exceedance probabilities are compared across a displacement grid, ",
        "and quantitative metrics are computed to assess numerical agreement.",
        "",
        "## Summary Table",
        "",
        summary_to_markdown(summary),
        "",
        "## Overall Statistics",
        "",
        f"- **Total test cases:** {len(df)}",
        f"- **Models validated:** {df['model_name'].nunique()}",
        f"- **Tests passed:** {df['within_tolerance'].sum()} ({df['within_tolerance'].mean()*100:.1f}%)",
        f"- **Max relative difference (all):** {df['max_rel_diff'].max():.2e}",
        f"- **Mean relative difference (all):** {df['max_rel_diff'].mean():.2e}",
        "",
    ]
    
    # Per-model sections
    for model_name in df['model_name'].unique():
        model_df = df[df['model_name'] == model_name]
        model_summary = summary[summary['model_name'] == model_name]
        
        lines.extend([
            f"## {model_name}",
            "",
        ])
        
        # Model statistics
        lines.extend([
            f"- **Variants tested:** {', '.join(model_df['variant'].unique())}",
            f"- **Total test cases:** {len(model_df)}",
            f"- **Tests passed:** {model_df['within_tolerance'].sum()} / {len(model_df)}",
            f"- **Max relative difference:** {model_df['max_rel_diff'].max():.2e}",
            "",
        ])
        
        # Representative plots
        if model_name in representative_plots:
            lines.append("### Representative Hazard Curve Comparisons")
            lines.append("")
            for plot_file in representative_plots[model_name]:
                rel_path = f"plots/{plot_file}"
                lines.append(f"![{plot_file}]({rel_path})")
                lines.append("")
        
        lines.append("")
    
    # Footer
    lines.extend([
        "---",
        "",
        "## Notes",
        "",
        "- Tolerances vary by model: Youngs 2003 uses `rtol=1e-4` due to integration method differences; ",
        "  all other models use `rtol=1e-6`.",
        "- Relative differences are computed as `|pfdha - fdhpy| / |fdhpy|` where `|fdhpy| > 1e-12`.",
        "- All comparisons use the same displacement grid as the automated pytest suite.",
        "",
    ])
    
    return "\n".join(lines)


# ============================================================================
# CSV/JSON Export
# ============================================================================

def export_results(
    df: pd.DataFrame,
    summary: pd.DataFrame,
    output_dir: Path,
) -> None:
    """Export results to CSV and JSON files."""
    
    # Detailed results
    df.to_csv(output_dir / "comparison_results.csv", index=False)
    
    # JSON version
    results_dict = {
        "generated": datetime.now().isoformat(),
        "total_tests": len(df),
        "tests_passed": int(df['within_tolerance'].sum()),
        "models": df['model_name'].unique().tolist(),
        "results": df.to_dict(orient='records'),
    }
    with open(output_dir / "comparison_results.json", 'w') as f:
        json.dump(results_dict, f, indent=2, default=str)
    
    # Summary CSV
    summary.to_csv(output_dir / "comparison_summary.csv", index=False)
    
    # Summary JSON
    summary_dict = {
        "generated": datetime.now().isoformat(),
        "summary": summary.to_dict(orient='records'),
    }
    with open(output_dir / "comparison_summary.json", 'w') as f:
        json.dump(summary_dict, f, indent=2, default=str)

















