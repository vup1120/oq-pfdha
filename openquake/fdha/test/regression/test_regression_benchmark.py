# -*- coding: utf-8 -*-
"""
Regression Test Suite for PFDHA Calculators - JSON Output Comparison (v5)

This module compares COMPLETE JSON output files, not just single arrays.
This provides more comprehensive validation and human-readable golden files.

What the golden files ARE (and are not):
- They are SELF-REFERENTIAL snapshots of THIS code's own output, frozen by
  running with GENERATE_GOLDEN=1. Their job is to catch ACCIDENTAL drift:
  "does the code still produce the same numbers it did last time?"
- They are NOT reference values from the model authors. Author/paper
  benchmarks live in the benchmark/ suites and are compared within
  published-tolerance bands - do not confuse the two.
- When you INTENTIONALLY change the numerics (e.g. a more accurate distance
  frame), the code's output legitimately shifts, so you re-freeze these
  goldens with GENERATE_GOLDEN=1 and review the diff. Re-freezing after a
  deliberate change is the designed workflow, not a violation - the guard
  here is against UNINTENDED changes.

Advantages of JSON comparison:
- Complete output validation (all keys, not just poes)
- Human-readable golden files (can be inspected in any editor)
- Git diff friendly (can see exactly what changed)
- Includes metadata (site coordinates, fault contributions, etc.)

Usage:
------
1. Generate Golden Truths:
   GENERATE_GOLDEN=1 pytest test_regression_benchmark.py -v

2. Run Regression Check:
   pytest test_regression_benchmark.py -v

3. View golden files:
   cat openquake/fdha/test/regression/data/golden_minimal_curve.json
"""

import os
import json
import tempfile
import pytest
import numpy as np
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, List, Dict, Any, Tuple

# =============================================================================
# Constants
# =============================================================================

ENV_GENERATE_GOLDEN = "GENERATE_GOLDEN"
REGRESSION_RTOL = 1e-7  # Strict tolerance
REGRESSION_ATOL = 0.0
MONTE_CARLO_RTOL = 1e-3  # Relaxed for stochastic models
DATA_DIR = Path(__file__).parent / "data"

# =============================================================================
# Benchmark Definitions
# =============================================================================

@dataclass
class BenchmarkConfig:
    name: str
    config_path: str
    calculator_type: str  # 'curve', 'curve_rectangular', 'map'
    description: str
    fixed_sites: Optional[List[Tuple[float, float]]] = None
    source_model_override: Optional[str] = None
    rtol_override: Optional[float] = None
    # Keys to compare (None = compare all numeric keys)
    compare_keys: Optional[List[str]] = None


BENCHMARKS = [
    # --- 1. Basic Functionality ---
    BenchmarkConfig(
        name="minimal_curve",
        config_path="examples/hazard_curve_minimal.ini",
        calculator_type="curve",
        description="Minimal hazard curve example (Youngs 2003).",
    ),

    # --- 2. Norcia (Youngs 2003) ---
    BenchmarkConfig(
        name="norcia_youngs2003",
        config_path="openquake/fdha/test/benchmark/norcia_sensitivity_youngs2003/job_norcia_youngs2003_AD_85.ini",
        calculator_type="curve",
        description="Norcia sensitivity (Youngs 2003 AD 85th).",
    ),

    # --- 3. Kumamoto (Chiou 2025) ---
    BenchmarkConfig(
        name="kumamoto_chiou2025",
        config_path="openquake/fdha/test/benchmark/valentini_et_al_2025/data/configuration/job_kumamoto_case2_Chiou2025.ini",
        calculator_type="curve",
        description="Kumamoto Case 2 (Chiou 2025 Primary FD).",
    ),

    # --- 4. FDHI Model Demos ---
    BenchmarkConfig(
        name="demo_kuehn2024",
        config_path="openquake/fdha/demo/hazard_curve/job_kuehn2024.ini",
        calculator_type="curve",
        description="Demo: Kuehn 2024 Primary FD.",
    ),
    BenchmarkConfig(
        name="demo_moss2024",
        config_path="openquake/fdha/demo/hazard_curve/job_moss2024.ini",
        calculator_type="curve",
        description="Demo: Moss 2024 Primary FD.",
    ),
]

# =============================================================================
# Utility Functions
# =============================================================================

def get_golden_path(benchmark_name: str) -> Path:
    """Get path to golden JSON file."""
    return DATA_DIR / f"golden_{benchmark_name}.json"


def should_generate_golden() -> bool:
    return os.environ.get(ENV_GENERATE_GOLDEN, "").lower() in ("1", "true", "yes")


def ensure_data_dir():
    DATA_DIR.mkdir(parents=True, exist_ok=True)


def resolve_config_path(config_path: str) -> Path:
    """Find config file by searching up directory tree."""
    potential_roots = [
        Path(__file__).parent.parent.parent.parent.parent,
        Path(__file__).parent.parent.parent.parent,
        Path(__file__).parent.parent.parent,
        Path.cwd(),
    ]
    for root in potential_roots:
        full_path = root / config_path
        if full_path.exists():
            return full_path.resolve()
    raise FileNotFoundError(f"Config not found: {config_path}")


def find_source_model(config_path: Path, config: Dict[str, Any]) -> List[str]:
    """Find source model files from config."""
    config_dir = config_path.parent
    calc_cfg = config.get('calculation', config.get('parameters', {}))
    source_file = calc_cfg.get('source_model_file', calc_cfg.get('source_model', ''))
    
    if source_file:
        files = [source_file] if isinstance(source_file, str) else list(source_file)
        resolved = []
        for f in files:
            p = Path(f)
            if p.is_absolute() and p.exists():
                resolved.append(str(p))
            elif (config_dir / f).exists():
                resolved.append(str((config_dir / f).resolve()))
            else:
                raise ValueError(f"Source model not found: {f}")
        return resolved
    raise ValueError(f"No source model in config: {config_path}")


# =============================================================================
# JSON Serialization Helpers
# =============================================================================

def numpy_to_json_serializable(obj: Any) -> Any:
    """Convert numpy types to JSON-serializable Python types."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64, np.float32)):
        return float(obj)
    elif isinstance(obj, dict):
        return {k: numpy_to_json_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [numpy_to_json_serializable(item) for item in obj]
    return obj


def save_result_json(result: Dict[str, Any], path: Path):
    """Save result dictionary to JSON file."""
    serializable = numpy_to_json_serializable(result)
    with open(path, 'w') as f:
        json.dump(serializable, f, indent=2)


def load_result_json(path: Path) -> Dict[str, Any]:
    """Load result dictionary from JSON file."""
    with open(path, 'r') as f:
        return json.load(f)


# =============================================================================
# JSON Comparison Functions
# =============================================================================

def compare_values(current: Any, golden: Any, path: str, rtol: float, atol: float) -> List[str]:
    """
    Recursively compare two values and return list of differences.
    
    Returns:
        List of difference descriptions (empty if values match)
    """
    differences = []
    
    if isinstance(golden, dict) and isinstance(current, dict):
        # Compare dictionaries
        all_keys = set(golden.keys()) | set(current.keys())
        for key in all_keys:
            key_path = f"{path}.{key}" if path else key
            if key not in golden:
                differences.append(f"NEW KEY: {key_path} (not in golden)")
            elif key not in current:
                differences.append(f"MISSING KEY: {key_path} (not in current)")
            else:
                differences.extend(compare_values(current[key], golden[key], key_path, rtol, atol))
    
    elif isinstance(golden, list) and isinstance(current, list):
        # Compare lists
        if len(golden) != len(current):
            differences.append(f"LENGTH MISMATCH at {path}: golden={len(golden)}, current={len(current)}")
        else:
            for i, (c, g) in enumerate(zip(current, golden)):
                differences.extend(compare_values(c, g, f"{path}[{i}]", rtol, atol))
    
    elif isinstance(golden, (int, float)) and isinstance(current, (int, float)):
        # Compare numeric values
        if golden == 0 and current == 0:
            pass  # Both zero, OK
        elif golden == 0:
            if abs(current) > atol:
                differences.append(f"VALUE MISMATCH at {path}: golden=0, current={current:.6e}")
        else:
            rel_diff = abs(current - golden) / abs(golden)
            abs_diff = abs(current - golden)
            if rel_diff > rtol and abs_diff > atol:
                differences.append(
                    f"VALUE MISMATCH at {path}: golden={golden:.6e}, current={current:.6e}, "
                    f"rel_diff={rel_diff:.2e}, abs_diff={abs_diff:.2e}"
                )
    
    elif isinstance(golden, str) and isinstance(current, str):
        # Compare strings exactly
        if golden != current:
            differences.append(f"STRING MISMATCH at {path}: golden='{golden}', current='{current}'")
    
    elif type(golden) != type(current):
        differences.append(f"TYPE MISMATCH at {path}: golden={type(golden).__name__}, current={type(current).__name__}")
    
    return differences


def compare_results(current: Dict, golden: Dict, rtol: float, atol: float, 
                    compare_keys: Optional[List[str]] = None) -> List[str]:
    """
    Compare current result with golden result.
    
    Args:
        current: Current calculation result
        golden: Golden truth result
        rtol: Relative tolerance
        atol: Absolute tolerance
        compare_keys: If specified, only compare these keys (None = all keys)
    
    Returns:
        List of differences (empty if results match)
    """
    if compare_keys:
        # Only compare specified keys
        current_filtered = {k: current[k] for k in compare_keys if k in current}
        golden_filtered = {k: golden[k] for k in compare_keys if k in golden}
        return compare_values(current_filtered, golden_filtered, "", rtol, atol)
    else:
        # Compare all keys
        return compare_values(current, golden, "", rtol, atol)


# =============================================================================
# Calculator Runners
# =============================================================================

def run_curve_calculator(bm: BenchmarkConfig) -> Dict[str, Any]:
    """Run curve calculator and return full result dict."""
    from openquake.fdha.calc.calculators import FaultRuptureProbabilityCalculator
    from openquake.fdha.calc.config_loader import load_config
    
    cfg_path = resolve_config_path(bm.config_path)
    config = load_config(str(cfg_path))
    sources = [bm.source_model_override] if bm.source_model_override else find_source_model(cfg_path, config)
    
    calc = FaultRuptureProbabilityCalculator(str(cfg_path), sources)
    return calc.run()


def run_curve_rectangular_calculator(bm: BenchmarkConfig) -> Dict[str, Any]:
    """Run rectangular site calculator and return full result dict.
    
    NOTE: Rectangular site support has been removed. This function raises an error.
    """
    raise NotImplementedError(
        "Rectangular site (corner_points) support has been removed. "
        "Please use point sites (latitude/longitude) instead."
    )


def run_map_calculator(bm: BenchmarkConfig) -> Dict[str, Any]:
    """Run map calculator with fixed sites and return full result dict."""
    from openquake.fdha.calc.calculators import BaseFaultRuptureCalculator
    from openquake.fdha.calc.hazard import calculate_fdha_hazard
    from openquake.fdha.calc.config_loader import load_config
    from openquake.hazardlib.site import Site, SiteCollection
    from openquake.hazardlib.geo import Point
    
    if not bm.fixed_sites:
        raise ValueError(f"Map benchmark '{bm.name}' requires fixed_sites.")
    
    cfg_path = resolve_config_path(bm.config_path)
    config = load_config(str(cfg_path))
    sources = [bm.source_model_override] if bm.source_model_override else find_source_model(cfg_path, config)
    
    calc = BaseFaultRuptureCalculator(str(cfg_path), sources)
    
    # Use fixed sites for deterministic regression
    sites = [Site(Point(lon, lat)) for lon, lat in bm.fixed_sites]
    sitecol = SiteCollection(sites)
    
    return calculate_fdha_hazard(calc, sitecol)


def run_benchmark(bm: BenchmarkConfig) -> Dict[str, Any]:
    """Dispatch to appropriate calculator and return full result dict."""
    cfg_path = resolve_config_path(bm.config_path)

    config = None
    try:
        from openquake.fdha.calc.config_loader import (
            calculation_requests_fdha_logic_tree,
            load_config,
        )
        config = load_config(str(cfg_path))
    except Exception:
        config = None

    calc_cfg = config.get("calculation", {}) if isinstance(config, dict) else {}
    if calculation_requests_fdha_logic_tree(calc_cfg):
        return run_logic_tree_calculator(cfg_path)

    if bm.calculator_type == "curve":
        return run_curve_calculator(bm)
    elif bm.calculator_type == "curve_rectangular":
        return run_curve_rectangular_calculator(bm)
    elif bm.calculator_type == "map":
        return run_map_calculator(bm)
    raise ValueError(f"Unknown calculator type: {bm.calculator_type}")


def run_logic_tree_calculator(config_path: Path) -> Dict[str, Any]:
    """Run a canonical logic-tree INI and normalize output to regression JSON."""
    from openquake.fdha.logic_tree.driver import FdhaLogicTree

    with tempfile.TemporaryDirectory(prefix="pfdha-regression-") as tmp:
        result = FdhaLogicTree.from_ini(config_path).run(outdir=tmp)

        if result.mode == "hazard_curve":
            return {
                "imls": list(result.d0),
                "rates": result.mean_rates,
            }

        if result.mode == "hazard_map":
            if result.displ_mean is None:
                raise ValueError(f"Map result missing displacement mean: {config_path}")
            return {
                "displacements": list(result.d0),
                "annual_rate_total": result.displ_mean,
            }

        raise ValueError(f"Unsupported logic-tree result mode: {result.mode}")


# =============================================================================
# Core Test Logic
# =============================================================================

def regression_test_benchmark(benchmark: BenchmarkConfig):
    """Core regression test with JSON comparison."""
    ensure_data_dir()
    golden_path = get_golden_path(benchmark.name)
    
    try:
        if should_generate_golden() or not golden_path.exists():
            # === GENERATE MODE ===
            print(f"\n[GENERATE] {benchmark.name}")
            result = run_benchmark(benchmark)
            save_result_json(result, golden_path)
            print(f"[GENERATE] Saved: {golden_path}")
            
            # Print summary
            print(f"[GENERATE] Keys: {list(result.keys())}")
            if 'rates' in result:
                rates = np.array(result['rates'])
                print(f"[GENERATE] rates shape: {rates.shape}, range: [{rates.min():.6e}, {rates.max():.6e}]")
            
            if should_generate_golden():
                pytest.skip(f"Generated golden file for {benchmark.name}")
        
        else:
            # === COMPARE MODE ===
            print(f"\n[COMPARE] {benchmark.name}")
            golden = load_result_json(golden_path)
            result = run_benchmark(benchmark)
            
            rtol = benchmark.rtol_override or REGRESSION_RTOL
            
            differences = compare_results(result, golden, rtol, REGRESSION_ATOL, benchmark.compare_keys)
            
            if differences:
                print(f"\n{'='*60}")
                print(f"REGRESSION FAILURE: {benchmark.name}")
                print(f"{'='*60}")
                print(f"Found {len(differences)} difference(s):")
                for i, diff in enumerate(differences[:20], 1):  # Show first 20
                    print(f"  {i}. {diff}")
                if len(differences) > 20:
                    print(f"  ... and {len(differences) - 20} more differences")
                print(f"{'='*60}")
                print(f"Golden file: {golden_path}")
                print(f"To regenerate: GENERATE_GOLDEN=1 pytest ... -k '{benchmark.name}'")
                print(f"{'='*60}\n")
                
                pytest.fail(f"Regression failure: {len(differences)} differences found")
            
            print(f"[PASS] {benchmark.name} - all values match within tolerance")

    except FileNotFoundError as e:
        pytest.fail(f"Missing regression input for {benchmark.name}: {e}")


# =============================================================================
# Explicit Test Functions
# =============================================================================

def test_minimal_curve():
    regression_test_benchmark(next(b for b in BENCHMARKS if b.name == "minimal_curve"))

def test_norcia_youngs2003():
    regression_test_benchmark(next(b for b in BENCHMARKS if b.name == "norcia_youngs2003"))

def test_kumamoto_chiou2025():
    regression_test_benchmark(next(b for b in BENCHMARKS if b.name == "kumamoto_chiou2025"))

def test_demo_kuehn2024():
    regression_test_benchmark(next(b for b in BENCHMARKS if b.name == "demo_kuehn2024"))

def test_demo_moss2024():
    regression_test_benchmark(next(b for b in BENCHMARKS if b.name == "demo_moss2024"))

# =============================================================================
# Main Entry Point
# =============================================================================

if __name__ == "__main__":
    import sys
    if "--generate" in sys.argv:
        os.environ[ENV_GENERATE_GOLDEN] = "1"
        sys.argv.remove("--generate")
    sys.exit(pytest.main([__file__, "-v"]))
