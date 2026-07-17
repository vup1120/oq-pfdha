#!/usr/bin/env python3
"""
PFDHA Comparison Dashboard Generator

This script generates a comprehensive comparison dashboard between pfdha
(OpenQuake implementation) and fdhpy (FDHI reference implementation).

Outputs (in outputs/ folder):
- plots/              : PNG hazard curve comparison plots
- comparison_results.csv/json : Detailed per-test-case metrics
- comparison_summary.csv/json : Aggregated per-model summary
- comparison_dashboard.md     : Full Markdown report with embedded figures

Usage:
    cd openquake/fdha/test/FDHI_Tests
    python generate_comparison_dashboard.py

Requirements:
    - fdhpy installed (pip install fdhpy)
    - pfdha installed (pip install -e /path/to/pfdha)
    - matplotlib installed
"""

import sys
from pathlib import Path
from typing import List, Dict, Any

import numpy as np

# ============================================================================
# Verify dependencies
# ============================================================================

try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
except ImportError:
    print("[ERROR] matplotlib is required: pip install matplotlib")
    sys.exit(1)

try:
    import pandas as pd
except ImportError:
    print("[ERROR] pandas is required: pip install pandas")
    sys.exit(1)

try:
    import fdhpy
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
        Lavrentiadis2023PrimaryFD_aggregate,
        Chiou2025PrimaryFD,
    )
except ImportError as e:
    print(f"[ERROR] pfdha models not importable: {e}")
    print("Install with: pip install -e /path/to/pfdha")
    sys.exit(1)

# Import local utilities
from report_utils import (
    DISPLACEMENTS,
    DEFAULT_TOLERANCES,
    get_parameter_grids,
    ComparisonMetrics,
    TestCaseResult,
    create_hazard_curve_plot,
    create_summary_plot,
    results_to_dataframe,
    generate_model_summary,
    generate_markdown_report,
    export_results,
)


# ============================================================================
# Model Runners
# ============================================================================

def run_youngs2003(
    magnitude: float,
    xl: float,
    displacements: np.ndarray,
) -> TestCaseResult:
    """Run Youngs 2003 comparison."""
    # fdhpy
    fdhpy_model = YoungsEtAl2003(
        magnitude=magnitude,
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
        mag=magnitude,
        style="all",
        norm_disp_type="AD"
    ).flatten()
    
    tol = DEFAULT_TOLERANCES["youngs2003"]
    metrics = ComparisonMetrics.compute(fdhpy_probs, pfdha_probs, **tol)
    
    return TestCaseResult(
        model_name="Youngs2003",
        variant="D/AD",
        magnitude=magnitude,
        x_l_ratio=xl,
        extra_params={},
        metrics=metrics,
        fdhpy_probs=fdhpy_probs,
        pfdha_probs=pfdha_probs,
        displacements=displacements,
        tolerance_rtol=tol["rtol"],
        tolerance_atol=tol["atol"],
    )


def run_petersen2011(
    magnitude: float,
    xl: float,
    version: str,
    displacements: np.ndarray,
) -> TestCaseResult:
    """Run Petersen 2011 comparison."""
    # fdhpy
    fdhpy_model = PetersenEtAl2011(
        magnitude=magnitude,
        xl=xl,
        version=version,
        displ_array=displacements,
    )
    fdhpy_probs = fdhpy_model.prob_exceed
    
    # pfdha
    pfdha_model = Petersen2011PrimaryFD()
    pfdha_probs = pfdha_model.get_prob(
        d=displacements,
        X_L_ratio=np.array([xl]),
        mag=magnitude,
        version=version
    ).flatten()
    
    tol = DEFAULT_TOLERANCES["petersen2011"]
    metrics = ComparisonMetrics.compute(fdhpy_probs, pfdha_probs, **tol)
    
    return TestCaseResult(
        model_name="Petersen2011",
        variant=version,
        magnitude=magnitude,
        x_l_ratio=xl,
        extra_params={"version": version},
        metrics=metrics,
        fdhpy_probs=fdhpy_probs,
        pfdha_probs=pfdha_probs,
        displacements=displacements,
        tolerance_rtol=tol["rtol"],
        tolerance_atol=tol["atol"],
    )


def run_moss2024(
    magnitude: float,
    xl: float,
    fdhpy_version: str,
    norm_type: str,
    use_girs: bool,
    displacements: np.ndarray,
) -> TestCaseResult:
    """Run Moss 2024 comparison."""
    source = "GIRS" if use_girs else "EQS"
    
    # fdhpy
    fdhpy_model = MossEtAl2024(
        magnitude=magnitude,
        xl=xl,
        version=fdhpy_version,
        displ_array=displacements,
        use_girs=use_girs,
        complete=True,
    )
    fdhpy_probs = fdhpy_model.prob_exceed
    
    # pfdha
    pfdha_model = Moss2024PrimaryFD()
    pfdha_probs = pfdha_model.get_prob(
        d=displacements,
        X_L_ratio=np.array([xl]),
        mag=magnitude,
        version=norm_type,
        source=source,
        completeness="complete"
    ).flatten()
    
    tol = DEFAULT_TOLERANCES["moss2024"]
    metrics = ComparisonMetrics.compute(fdhpy_probs, pfdha_probs, **tol)
    
    variant = f"{fdhpy_version}_{source}"
    
    return TestCaseResult(
        model_name="Moss2024",
        variant=variant,
        magnitude=magnitude,
        x_l_ratio=xl,
        extra_params={"fdhpy_version": fdhpy_version, "norm_type": norm_type, "source": source},
        metrics=metrics,
        fdhpy_probs=fdhpy_probs,
        pfdha_probs=pfdha_probs,
        displacements=displacements,
        tolerance_rtol=tol["rtol"],
        tolerance_atol=tol["atol"],
    )


def run_kuehn2024(
    magnitude: float,
    xl: float,
    style: str,
    displacements: np.ndarray,
) -> TestCaseResult:
    """Run Kuehn 2024 comparison."""
    # fdhpy
    fdhpy_model = KuehnEtAl2024(
        style=style,
        magnitude=magnitude,
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
        mag=magnitude,
        style=style,
        folded=True,
        epistemic_uncertainty=False
    ).flatten()
    
    tol = DEFAULT_TOLERANCES["kuehn2024"]
    metrics = ComparisonMetrics.compute(fdhpy_probs, pfdha_probs, **tol)
    
    return TestCaseResult(
        model_name="Kuehn2024",
        variant=style,
        magnitude=magnitude,
        x_l_ratio=xl,
        extra_params={"style": style},
        metrics=metrics,
        fdhpy_probs=fdhpy_probs,
        pfdha_probs=pfdha_probs,
        displacements=displacements,
        tolerance_rtol=tol["rtol"],
        tolerance_atol=tol["atol"],
    )


def run_lavrentiadis2023(
    magnitude: float,
    xl: float,
    include_prob_zero: bool,
    displacements: np.ndarray,
) -> TestCaseResult:
    """Run Lavrentiadis 2023 comparison."""
    # fdhpy
    fdhpy_model = LavrentiadisAbrahamson2023(
        magnitude=magnitude,
        xl=xl,
        displ_array=displacements,
        metric="aggregate",
        version="full rupture",
        style="strike-slip",
        include_prob_zero=include_prob_zero,
    )
    fdhpy_probs = fdhpy_model.prob_exceed
    
    # pfdha
    pfdha_model = Lavrentiadis2023PrimaryFD_aggregate()
    pfdha_probs = pfdha_model.get_prob(
        d=displacements,
        X_L_ratio=np.array([xl]),
        mag=magnitude,
        style="strike-slip",
        output_type="disp_agg_prime",
        include_zero_slip=include_prob_zero,
    ).flatten()
    
    tol = DEFAULT_TOLERANCES["lavrentiadis2023"]
    metrics = ComparisonMetrics.compute(fdhpy_probs, pfdha_probs, **tol)
    
    variant = f"aggregate_{'with_zero' if include_prob_zero else 'no_zero'}"
    
    return TestCaseResult(
        model_name="Lavrentiadis2023",
        variant=variant,
        magnitude=magnitude,
        x_l_ratio=xl,
        extra_params={"include_prob_zero": include_prob_zero},
        metrics=metrics,
        fdhpy_probs=fdhpy_probs,
        pfdha_probs=pfdha_probs,
        displacements=displacements,
        tolerance_rtol=tol["rtol"],
        tolerance_atol=tol["atol"],
    )


def run_chiou2025(
    magnitude: float,
    xl: float,
    version: str,
    displacements: np.ndarray,
) -> TestCaseResult:
    """Run Chiou 2025 comparison."""
    # fdhpy
    fdhpy_model = ChiouEtAl2025(
        magnitude=magnitude,
        xl=xl,
        version=version,
        displ_array=displacements,
    )
    fdhpy_probs = fdhpy_model.prob_exceed
    
    # pfdha
    pfdha_model = Chiou2025PrimaryFD()
    pfdha_probs = pfdha_model.get_prob(
        d=displacements,
        X_L_ratio=np.array([xl]),
        mag=magnitude,
        version=version
    ).flatten()
    
    tol = DEFAULT_TOLERANCES["chiou2025"]
    metrics = ComparisonMetrics.compute(fdhpy_probs, pfdha_probs, **tol)
    
    return TestCaseResult(
        model_name="Chiou2025",
        variant=version,
        magnitude=magnitude,
        x_l_ratio=xl,
        extra_params={"version": version},
        metrics=metrics,
        fdhpy_probs=fdhpy_probs,
        pfdha_probs=pfdha_probs,
        displacements=displacements,
        tolerance_rtol=tol["rtol"],
        tolerance_atol=tol["atol"],
    )


# ============================================================================
# Main Runner
# ============================================================================

def run_all_comparisons(displacements: np.ndarray) -> List[TestCaseResult]:
    """Run all model comparisons and return results."""
    results = []
    grids = get_parameter_grids()
    
    print("\nRunning comparisons...")
    
    # Youngs 2003
    print("  Youngs 2003...", end=" ", flush=True)
    count = 0
    for mag in grids["youngs2003"]["magnitudes"]:
        for xl in grids["youngs2003"]["x_l_ratios"]:
            results.append(run_youngs2003(mag, xl, displacements))
            count += 1
    print(f"{count} cases")
    
    # Petersen 2011
    print("  Petersen 2011...", end=" ", flush=True)
    count = 0
    for version in grids["petersen2011"]["versions"]:
        for mag in grids["petersen2011"]["magnitudes"]:
            for xl in grids["petersen2011"]["x_l_ratios"]:
                results.append(run_petersen2011(mag, xl, version, displacements))
                count += 1
    print(f"{count} cases")
    
    # Moss 2024
    print("  Moss 2024...", end=" ", flush=True)
    count = 0
    for fdhpy_ver, norm_type in grids["moss2024"]["versions"]:
        for use_girs in grids["moss2024"]["use_girs"]:
            for mag in grids["moss2024"]["magnitudes"]:
                for xl in grids["moss2024"]["x_l_ratios"]:
                    results.append(run_moss2024(mag, xl, fdhpy_ver, norm_type, use_girs, displacements))
                    count += 1
    print(f"{count} cases")
    
    # Kuehn 2024
    print("  Kuehn 2024...", end=" ", flush=True)
    count = 0
    for style in grids["kuehn2024"]["styles"]:
        for mag in grids["kuehn2024"]["magnitudes"]:
            for xl in grids["kuehn2024"]["x_l_ratios"]:
                results.append(run_kuehn2024(mag, xl, style, displacements))
                count += 1
    print(f"{count} cases")
    
    # Lavrentiadis 2023
    print("  Lavrentiadis 2023...", end=" ", flush=True)
    count = 0
    for include_zero in grids["lavrentiadis2023"]["include_prob_zero"]:
        for mag in grids["lavrentiadis2023"]["magnitudes"]:
            for xl in grids["lavrentiadis2023"]["x_l_ratios"]:
                results.append(run_lavrentiadis2023(mag, xl, include_zero, displacements))
                count += 1
    print(f"{count} cases")
    
    # Chiou 2025
    print("  Chiou 2025...", end=" ", flush=True)
    count = 0
    for version in grids["chiou2025"]["versions"]:
        for mag in grids["chiou2025"]["magnitudes"]:
            for xl in grids["chiou2025"]["x_l_ratios"]:
                results.append(run_chiou2025(mag, xl, version, displacements))
                count += 1
    print(f"{count} cases")
    
    return results


def select_representative_cases(
    results: List[TestCaseResult],
) -> Dict[str, List[TestCaseResult]]:
    """
    Select representative test cases for plotting.
    Choose mid-range M and x/L for each model/variant.
    """
    selected = {}
    
    for result in results:
        key = (result.model_name, result.variant)
        if key not in selected:
            selected[key] = []
        
        # Select mid-range cases
        is_mid_mag = result.magnitude in [7.0, 7.2]
        is_mid_xl = result.x_l_ratio in [0.5, 0.3, 0.25]
        
        if is_mid_mag and is_mid_xl:
            selected[key].append(result)
    
    # Flatten to model-level
    by_model = {}
    for (model_name, variant), cases in selected.items():
        if model_name not in by_model:
            by_model[model_name] = []
        # Take up to 2 cases per variant
        by_model[model_name].extend(cases[:2])
    
    return by_model


def generate_plots(
    results: List[TestCaseResult],
    plots_dir: Path,
) -> Dict[str, List[str]]:
    """Generate hazard curve plots for representative cases."""
    
    representative = select_representative_cases(results)
    plot_files = {}
    
    print("\nGenerating plots...")
    
    for model_name, cases in representative.items():
        print(f"  {model_name}...", end=" ", flush=True)
        plot_files[model_name] = []
        
        for case in cases:
            # Create filename
            safe_variant = case.variant.replace("/", "_").replace(" ", "_")
            filename = f"{model_name}_{safe_variant}_M{case.magnitude}_XL{case.x_l_ratio}.png"
            filepath = plots_dir / filename
            
            create_hazard_curve_plot(case, filepath)
            plot_files[model_name].append(filename)
        
        print(f"{len(cases)} plots")
    
    return plot_files


def main():
    """Main entry point."""
    print("=" * 70)
    print("PFDHA Comparison Dashboard Generator")
    print("=" * 70)
    
    # Setup output directories
    script_dir = Path(__file__).parent
    output_dir = script_dir / "outputs"
    plots_dir = output_dir / "plots"
    
    output_dir.mkdir(exist_ok=True)
    plots_dir.mkdir(exist_ok=True)
    
    print(f"\nOutput directory: {output_dir}")
    
    # Run all comparisons
    results = run_all_comparisons(DISPLACEMENTS)
    
    print(f"\nTotal test cases: {len(results)}")
    
    # Convert to DataFrame
    df = results_to_dataframe(results)
    
    # Generate summary
    summary = generate_model_summary(df)
    
    # Export CSV/JSON
    print("\nExporting results...")
    export_results(df, summary, output_dir)
    print(f"  - comparison_results.csv")
    print(f"  - comparison_results.json")
    print(f"  - comparison_summary.csv")
    print(f"  - comparison_summary.json")
    
    # Generate plots
    plot_files = generate_plots(results, plots_dir)
    
    # Generate summary plot
    create_summary_plot(df, output_dir / "comparison_overview.png")
    print(f"  - comparison_overview.png")
    
    # Generate Markdown report
    print("\nGenerating Markdown report...")
    markdown = generate_markdown_report(summary, df, plots_dir, plot_files)
    
    md_path = output_dir / "comparison_dashboard.md"
    with open(md_path, 'w') as f:
        f.write(markdown)
    print(f"  - comparison_dashboard.md")
    
    # Print summary to console
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    print(f"\nModels validated: {df['model_name'].nunique()}")
    print(f"Total test cases: {len(df)}")
    print(f"Tests passed: {df['within_tolerance'].sum()} / {len(df)} ({df['within_tolerance'].mean()*100:.1f}%)")
    
    print("\nPer-model results:")
    for _, row in summary.iterrows():
        status = "✓" if row['pass_rate'] == 1.0 else "✗"
        print(f"  {status} {row['model_name']:20} {row['variant']:20} "
              f"max_rel={row['max_rel_diff_max']:.2e} "
              f"({int(row['tests_passed'])}/{int(row['total_tests'])} passed)")
    
    print(f"\nOutputs written to: {output_dir}")
    print("\nDone!")
    
    # Return success/failure
    if df['within_tolerance'].all():
        return 0
    else:
        print("\n[WARNING] Some tests exceeded tolerance!")
        return 1


if __name__ == "__main__":
    sys.exit(main())

















