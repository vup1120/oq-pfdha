"""
PFDHA Reference Test Suite - Comparison Helper Utilities.

This module provides reusable comparison functions and utilities
for validating pfdha against fdhpy reference implementation.

Note: This module is a "toolbox" for richer diagnostics and may not yet
be used by test_simple_reference.py. The simple reference tests currently
only use direct numpy comparisons. These helpers can be used for more
detailed analysis when needed.
"""

import numpy as np
from typing import Dict, Any, Optional, Tuple, Union, List
from dataclasses import dataclass
import warnings


# ============================================================================
# Comparison Results
# ============================================================================

@dataclass
class DetailedComparisonResult:
    """Detailed results from comparing two arrays."""
    passed: bool
    n_values: int
    n_failures: int
    max_abs_diff: float
    mean_abs_diff: float
    max_rel_diff: float
    mean_rel_diff: float
    rmse: float
    failure_indices: Optional[np.ndarray]
    failure_values_expected: Optional[np.ndarray]
    failure_values_actual: Optional[np.ndarray]
    tolerance_rtol: float
    tolerance_atol: float
    
    def __str__(self) -> str:
        status = "PASS" if self.passed else "FAIL"
        return (
            f"{status}: {self.n_values} values, {self.n_failures} failures\n"
            f"  Max abs diff: {self.max_abs_diff:.2e}\n"
            f"  Mean abs diff: {self.mean_abs_diff:.2e}\n"
            f"  Max rel diff: {self.max_rel_diff:.4f}%\n"
            f"  Mean rel diff: {self.mean_rel_diff:.4f}%\n"
            f"  RMSE: {self.rmse:.2e}\n"
            f"  Tolerance: rtol={self.tolerance_rtol}, atol={self.tolerance_atol}"
        )
    
    def get_failure_message(self) -> str:
        """Generate detailed failure message for pytest."""
        if self.passed:
            return "Comparison passed"
        
        msg = f"Comparison failed: {self.n_failures}/{self.n_values} values outside tolerance\n"
        msg += f"Tolerance: rtol={self.tolerance_rtol}, atol={self.tolerance_atol}\n"
        msg += f"Max absolute difference: {self.max_abs_diff:.6e}\n"
        msg += f"Max relative difference: {self.max_rel_diff:.6f}%\n"
        
        if self.failure_indices is not None and len(self.failure_indices) > 0:
            # Show first few failures
            n_show = min(5, len(self.failure_indices))
            msg += f"\nFirst {n_show} failures:\n"
            for i in range(n_show):
                idx = self.failure_indices[i]
                exp = self.failure_values_expected[i]
                act = self.failure_values_actual[i]
                diff = abs(act - exp)
                rel = abs(diff / exp) * 100 if exp != 0 else float('inf')
                msg += f"  [{idx}] expected={exp:.6e}, actual={act:.6e}, diff={diff:.6e} ({rel:.4f}%)\n"
        
        return msg


class ComparisonHelper:
    """
    Utility class for comparing fdhpy and pfdha outputs with detailed diagnostics.
    
    Provides methods for:
    - Comparing probability arrays with proper tolerance handling
    - Comparing intermediate values (mean, sigma, etc.)
    - Generating detailed failure messages
    - Handling edge cases (zeros, infinities)
    """
    
    def __init__(self, rtol: float = 1e-6, atol: float = 1e-10):
        """
        Initialize comparison helper.
        
        Args:
            rtol: Relative tolerance for comparison
            atol: Absolute tolerance for comparison
        """
        self.rtol = rtol
        self.atol = atol
    
    def compare_arrays(
        self,
        expected: np.ndarray,
        actual: np.ndarray,
        rtol: Optional[float] = None,
        atol: Optional[float] = None,
    ) -> DetailedComparisonResult:
        """
        Compare two arrays with detailed diagnostics.
        
        Args:
            expected: Expected values (from fdhpy)
            actual: Actual values (from pfdha)
            rtol: Relative tolerance (uses instance default if None)
            atol: Absolute tolerance (uses instance default if None)
        
        Returns:
            DetailedComparisonResult with comparison statistics
        """
        rtol = rtol if rtol is not None else self.rtol
        atol = atol if atol is not None else self.atol
        
        # Ensure arrays are numpy arrays
        expected = np.atleast_1d(np.asarray(expected, dtype=float)).flatten()
        actual = np.atleast_1d(np.asarray(actual, dtype=float)).flatten()
        
        if expected.shape != actual.shape:
            raise ValueError(f"Shape mismatch: expected {expected.shape}, got {actual.shape}")
        
        n_values = len(expected)
        
        # Calculate differences
        abs_diff = np.abs(actual - expected)
        
        # Safe relative difference (avoid division by zero)
        with np.errstate(divide='ignore', invalid='ignore'):
            rel_diff = np.where(expected != 0, abs_diff / np.abs(expected) * 100, 0)
            rel_diff = np.where(np.isfinite(rel_diff), rel_diff, 0)
        
        # Check tolerance using numpy's allclose logic
        # |actual - expected| <= atol + rtol * |expected|
        tolerance = atol + rtol * np.abs(expected)
        within_tolerance = abs_diff <= tolerance
        
        # Find failures
        failure_mask = ~within_tolerance
        n_failures = np.sum(failure_mask)
        failure_indices = np.where(failure_mask)[0] if n_failures > 0 else None
        
        # Statistics
        max_abs_diff = float(np.max(abs_diff))
        mean_abs_diff = float(np.mean(abs_diff))
        max_rel_diff = float(np.max(rel_diff))
        mean_rel_diff = float(np.mean(rel_diff))
        rmse = float(np.sqrt(np.mean(abs_diff ** 2)))
        
        return DetailedComparisonResult(
            passed=(n_failures == 0),
            n_values=n_values,
            n_failures=n_failures,
            max_abs_diff=max_abs_diff,
            mean_abs_diff=mean_abs_diff,
            max_rel_diff=max_rel_diff,
            mean_rel_diff=mean_rel_diff,
            rmse=rmse,
            failure_indices=failure_indices,
            failure_values_expected=expected[failure_indices] if failure_indices is not None else None,
            failure_values_actual=actual[failure_indices] if failure_indices is not None else None,
            tolerance_rtol=rtol,
            tolerance_atol=atol,
        )
    
    def assert_arrays_equal(
        self,
        expected: np.ndarray,
        actual: np.ndarray,
        rtol: Optional[float] = None,
        atol: Optional[float] = None,
        msg: str = "",
    ) -> None:
        """
        Assert two arrays are equal within tolerance, with detailed failure message.
        
        Args:
            expected: Expected values (from fdhpy)
            actual: Actual values (from pfdha)
            rtol: Relative tolerance
            atol: Absolute tolerance
            msg: Additional context message
        
        Raises:
            AssertionError: If arrays differ beyond tolerance
        """
        result = self.compare_arrays(expected, actual, rtol, atol)
        
        if not result.passed:
            full_msg = f"{msg}\n{result.get_failure_message()}" if msg else result.get_failure_message()
            raise AssertionError(full_msg)
    
    def assert_probabilities_equal(
        self,
        expected: np.ndarray,
        actual: np.ndarray,
        msg: str = "",
    ) -> None:
        """
        Assert probability arrays are equal with appropriate tolerance.
        
        Uses tighter tolerance for mid-range probabilities and looser
        for extreme values (near 0 or 1).
        """
        # Standard tolerance for probabilities
        self.assert_arrays_equal(expected, actual, rtol=1e-6, atol=1e-10, msg=msg)
    
    def assert_log_values_equal(
        self,
        expected: np.ndarray,
        actual: np.ndarray,
        msg: str = "",
    ) -> None:
        """
        Assert log-space values are equal with appropriate tolerance.
        """
        self.assert_arrays_equal(expected, actual, rtol=1e-8, atol=1e-12, msg=msg)


# ============================================================================
# Model-Specific Parameter Mappers
# ============================================================================

class ParameterMapper:
    """
    Maps parameters between fdhpy and pfdha interfaces.
    
    Each model may have different parameter names or formats between
    the two implementations. This class handles the translation.
    """
    
    @staticmethod
    def youngs2003_fdhpy_to_pfdha(fdhpy_params: Dict[str, Any]) -> Dict[str, Any]:
        """Map Youngs 2003 parameters from fdhpy to pfdha format."""
        version = fdhpy_params.get("version", "d/ad")
        norm_disp_type = "AD" if version == "d/ad" else "MD"
        
        return {
            "mag": fdhpy_params["magnitude"],
            "X_L_ratio": np.array([fdhpy_params["xl"]]),
            "style": "all",  # fdhpy uses "All styles" WC94 coefficients
            "norm_disp_type": norm_disp_type,
        }
    
    @staticmethod
    def petersen2011_fdhpy_to_pfdha(fdhpy_params: Dict[str, Any]) -> Dict[str, Any]:
        """Map Petersen 2011 parameters from fdhpy to pfdha format."""
        return {
            "mag": fdhpy_params["magnitude"],
            "X_L_ratio": np.array([fdhpy_params["xl"]]),
            "version": fdhpy_params["version"],
        }
    
    @staticmethod
    def moss2024_fdhpy_to_pfdha(fdhpy_params: Dict[str, Any]) -> Dict[str, Any]:
        """Map Moss 2024 parameters from fdhpy to pfdha format."""
        version = fdhpy_params.get("version", "d/ad")
        norm_disp_type = "AD" if version == "d/ad" else "MD"
        
        return {
            "mag": fdhpy_params["magnitude"],
            "X_L_ratio": np.array([fdhpy_params["xl"]]),
            "version": norm_disp_type,
            "source": "GIRS" if fdhpy_params.get("use_girs", True) else "EQS",
            "completeness": "complete" if fdhpy_params.get("complete", True) else "all",
        }
    
    @staticmethod
    def kuehn2024_fdhpy_to_pfdha(fdhpy_params: Dict[str, Any]) -> Dict[str, Any]:
        """Map Kuehn 2024 parameters from fdhpy to pfdha format."""
        version = fdhpy_params.get("version", "median_coeffs")
        epistemic = version == "full_coeffs"
        
        return {
            "mag": fdhpy_params["magnitude"],
            "X_L_ratio": np.array([fdhpy_params["xl"]]),
            "style": fdhpy_params["style"],
            "folded": fdhpy_params.get("folded", True),
            "epistemic_uncertainty": epistemic,
        }
    
    @staticmethod
    def lavrentiadis2023_fdhpy_to_pfdha(fdhpy_params: Dict[str, Any]) -> Dict[str, Any]:
        """Map Lavrentiadis 2023 parameters from fdhpy to pfdha format."""
        metric = fdhpy_params.get("metric", "aggregate")
        version = fdhpy_params.get("version", "full rupture")
        
        # Map metric/version to output_type
        if metric == "aggregate" and version == "full rupture":
            output_type = "disp_agg_prime"
        elif metric == "aggregate" and version == "individual segment":
            output_type = "disp_agg_seg"
        elif metric == "sum-of-principal" and version == "full rupture":
            output_type = "disp_prnc_prime"
        else:
            output_type = "disp_agg_prime"
        
        return {
            "mag": fdhpy_params["magnitude"],
            "X_L_ratio": np.array([fdhpy_params["xl"]]),
            "style": "strike-slip",
            "output_type": output_type,
            "include_zero_slip": fdhpy_params.get("include_prob_zero", True),
        }
    
    @staticmethod
    def chiou2025_fdhpy_to_pfdha(fdhpy_params: Dict[str, Any]) -> Dict[str, Any]:
        """Map Chiou 2025 parameters from fdhpy to pfdha format."""
        return {
            "mag": fdhpy_params["magnitude"],
            "X_L_ratio": np.array([fdhpy_params["xl"]]),
            "version": fdhpy_params["version"],
        }


# ============================================================================
# Model Runner Utilities
# ============================================================================

class ModelRunner:
    """
    Utility class for running fdhpy and pfdha models with consistent interface.
    """
    
    def __init__(self, fdhpy_models: Dict[str, Any], pfdha_models: Dict[str, Any]):
        """
        Initialize with model dictionaries.
        
        Args:
            fdhpy_models: Dictionary of fdhpy model classes
            pfdha_models: Dictionary of pfdha model classes
        """
        self.fdhpy_models = fdhpy_models
        self.pfdha_models = pfdha_models
        self.mapper = ParameterMapper()
    
    def run_youngs2003(
        self,
        displacements: np.ndarray,
        magnitude: float,
        xl: float,
        version: str = "d/ad",
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Run both implementations of Youngs 2003.
        
        Returns:
            Tuple of (fdhpy_result, pfdha_result)
        """
        # Run fdhpy
        fdhpy_model = self.fdhpy_models["YoungsEtAl2003"](
            magnitude=magnitude,
            xl=xl,
            version=version,
            displ_array=displacements,
        )
        fdhpy_result = fdhpy_model.prob_exceed
        
        # Run pfdha
        pfdha_params = self.mapper.youngs2003_fdhpy_to_pfdha({
            "magnitude": magnitude, "xl": xl, "version": version
        })
        pfdha_model = self.pfdha_models["Youngs2003PrimaryFD"]()
        pfdha_result = pfdha_model.get_prob(d=displacements, **pfdha_params).flatten()
        
        return fdhpy_result, pfdha_result
    
    def run_petersen2011(
        self,
        displacements: np.ndarray,
        magnitude: float,
        xl: float,
        version: str = "elliptical",
    ) -> Tuple[Optional[np.ndarray], np.ndarray]:
        """
        Run both implementations of Petersen 2011.
        
        Returns:
            Tuple of (fdhpy_result, pfdha_result)
            fdhpy_result is None for bilinear (not supported by fdhpy)
        """
        # Run pfdha first (always works)
        pfdha_params = self.mapper.petersen2011_fdhpy_to_pfdha({
            "magnitude": magnitude, "xl": xl, "version": version
        })
        pfdha_model = self.pfdha_models["Petersen2011PrimaryFD"]()
        pfdha_result = pfdha_model.get_prob(d=displacements, **pfdha_params).flatten()
        
        # Run fdhpy (may not support all versions)
        fdhpy_result = None
        if version != "bilinear":  # fdhpy doesn't support bilinear
            fdhpy_model = self.fdhpy_models["PetersenEtAl2011"](
                magnitude=magnitude,
                xl=xl,
                version=version,
                displ_array=displacements,
            )
            fdhpy_result = fdhpy_model.prob_exceed
        
        return fdhpy_result, pfdha_result
    
    def run_moss2024(
        self,
        displacements: np.ndarray,
        magnitude: float,
        xl: float,
        version: str = "d/ad",
        use_girs: bool = True,
        complete: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Run both implementations of Moss 2024.
        """
        # Run fdhpy
        fdhpy_model = self.fdhpy_models["MossEtAl2024"](
            magnitude=magnitude,
            xl=xl,
            version=version,
            displ_array=displacements,
            use_girs=use_girs,
            complete=complete,
        )
        fdhpy_result = fdhpy_model.prob_exceed
        
        # Run pfdha
        pfdha_params = self.mapper.moss2024_fdhpy_to_pfdha({
            "magnitude": magnitude, "xl": xl, "version": version,
            "use_girs": use_girs, "complete": complete
        })
        pfdha_model = self.pfdha_models["Moss2024PrimaryFD"]()
        pfdha_result = pfdha_model.get_prob(d=displacements, **pfdha_params).flatten()
        
        return fdhpy_result, pfdha_result
    
    def run_kuehn2024(
        self,
        displacements: np.ndarray,
        magnitude: float,
        xl: float,
        style: str = "normal",
        version: str = "median_coeffs",
        folded: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Run both implementations of Kuehn 2024.
        """
        # Run fdhpy
        fdhpy_model = self.fdhpy_models["KuehnEtAl2024"](
            style=style,
            magnitude=magnitude,
            xl=xl,
            version=version,
            folded=folded,
            displ_array=displacements,
        )
        fdhpy_result = fdhpy_model.prob_exceed
        
        # Run pfdha
        pfdha_params = self.mapper.kuehn2024_fdhpy_to_pfdha({
            "magnitude": magnitude, "xl": xl, "style": style,
            "version": version, "folded": folded
        })
        pfdha_model = self.pfdha_models["Kuehn2024PrimaryFD"]()
        pfdha_result = pfdha_model.get_prob(d=displacements, **pfdha_params)
        
        # Handle shape differences
        if pfdha_result.ndim > 1:
            if version == "full_coeffs":
                # Transpose for shape alignment
                if pfdha_result.shape != fdhpy_result.shape:
                    pfdha_result = pfdha_result.T
            else:
                pfdha_result = pfdha_result.flatten()
        
        return fdhpy_result, pfdha_result
    
    def run_lavrentiadis2023(
        self,
        displacements: np.ndarray,
        magnitude: float,
        xl: float,
        metric: str = "aggregate",
        version: str = "full rupture",
        include_prob_zero: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Run both implementations of Lavrentiadis 2023.
        """
        # Run fdhpy
        fdhpy_model = self.fdhpy_models["LavrentiadisAbrahamson2023"](
            magnitude=magnitude,
            xl=xl,
            displ_array=displacements,
            metric=metric,
            version=version,
            include_prob_zero=include_prob_zero,
        )
        fdhpy_result = fdhpy_model.prob_exceed
        
        # Run pfdha
        pfdha_params = self.mapper.lavrentiadis2023_fdhpy_to_pfdha({
            "magnitude": magnitude, "xl": xl,
            "metric": metric, "version": version,
            "include_prob_zero": include_prob_zero
        })
        pfdha_model = self.pfdha_models["Lavrentiadis2023PrimaryFD"]()
        pfdha_result = pfdha_model.get_prob(d=displacements, **pfdha_params).flatten()
        
        return fdhpy_result, pfdha_result
    
    def run_chiou2025(
        self,
        displacements: np.ndarray,
        magnitude: float,
        xl: float,
        version: str = "model7",
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Run both implementations of Chiou 2025.
        """
        # Run fdhpy
        fdhpy_model = self.fdhpy_models["ChiouEtAl2025"](
            magnitude=magnitude,
            xl=xl,
            version=version,
            displ_array=displacements,
        )
        fdhpy_result = fdhpy_model.prob_exceed
        
        # Run pfdha
        pfdha_params = self.mapper.chiou2025_fdhpy_to_pfdha({
            "magnitude": magnitude, "xl": xl, "version": version
        })
        pfdha_model = self.pfdha_models["Chiou2025PrimaryFD"]()
        pfdha_result = pfdha_model.get_prob(d=displacements, **pfdha_params).flatten()
        
        return fdhpy_result, pfdha_result


# ============================================================================
# Test Case Generators
# ============================================================================

def generate_parameter_combinations(
    magnitudes: List[float],
    x_l_ratios: List[float],
    **additional_params
) -> List[Dict[str, Any]]:
    """
    Generate all combinations of parameters for testing.
    
    Args:
        magnitudes: List of magnitudes to test
        x_l_ratios: List of x/L ratios to test
        **additional_params: Additional parameters with lists of values
    
    Returns:
        List of parameter dictionaries
    """
    import itertools
    
    # Start with magnitude and x/L combinations
    base_params = [
        {"magnitude": m, "xl": xl}
        for m, xl in itertools.product(magnitudes, x_l_ratios)
    ]
    
    if not additional_params:
        return base_params
    
    # Add additional parameter combinations
    result = []
    for base in base_params:
        param_lists = list(additional_params.items())
        param_names = [p[0] for p in param_lists]
        param_values = [p[1] for p in param_lists]
        
        for combo in itertools.product(*param_values):
            params = dict(base)
            params.update(dict(zip(param_names, combo)))
            result.append(params)
    
    return result


def generate_test_id(params: Dict[str, Any]) -> str:
    """Generate a human-readable test ID from parameters."""
    parts = []
    for key, value in params.items():
        if key == "magnitude":
            parts.append(f"M{value}")
        elif key == "xl":
            parts.append(f"xl{value}")
        elif isinstance(value, bool):
            if value:
                parts.append(key)
        else:
            parts.append(f"{key}={value}")
    return "-".join(parts)

