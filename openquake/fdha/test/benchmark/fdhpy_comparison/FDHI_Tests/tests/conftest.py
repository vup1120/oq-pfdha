"""
PFDHA Reference Test Suite - pytest configuration and fixtures.

This module provides comprehensive fixtures for comparing fdhpy (reference)
and pfdha (this implementation) outputs with proper CI/CD structure.

Note: Markers are defined in pytest.ini; do not duplicate them here.
"""

import os
import sys
from pathlib import Path
from typing import Dict, Any, Optional, Tuple, Callable
from dataclasses import dataclass
import warnings

import numpy as np
import pytest

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False
    yaml = None

# ============================================================================
# Path Configuration
# ============================================================================

TESTS_DIR = Path(__file__).parent
FDHI_TESTS_DIR = TESTS_DIR.parent
CONFIG_DIR = TESTS_DIR

# Expected data directories (derived from installed packages when available)
# These are placeholders; tests should handle missing directories gracefully.
try:
    import fdhpy
    FDHPY_EXPECTED_DIR = Path(fdhpy.__file__).resolve().parent.parent / "tests" / "expected"
except ImportError:
    FDHPY_EXPECTED_DIR = Path("expected")  # Placeholder if fdhpy not installed

try:
    import openquake.fdha
    PFDHA_EXPECTED_DIR = Path(openquake.fdha.__file__).resolve().parent / "test" / "expected"
except ImportError:
    PFDHA_EXPECTED_DIR = Path("expected")  # Placeholder if pfdha not installed


# ============================================================================
# Tolerance Configuration
# ============================================================================

@dataclass
class ToleranceConfig:
    """Tolerance configuration for a specific comparison type."""
    rtol: float  # Relative tolerance
    atol: float  # Absolute tolerance
    description: str = ""


# Default tolerances by value type
TOLERANCES = {
    # Final probability values
    "prob_exceed": ToleranceConfig(rtol=1e-6, atol=1e-10, description="Exceedance probability"),
    "prob_zero": ToleranceConfig(rtol=1e-6, atol=1e-10, description="Zero-slip probability"),
    
    # Intermediate values
    "mean_ln": ToleranceConfig(rtol=1e-8, atol=1e-12, description="Mean in ln space"),
    "sigma_ln": ToleranceConfig(rtol=1e-8, atol=1e-12, description="Sigma in ln space"),
    "mean_log10": ToleranceConfig(rtol=1e-8, atol=1e-12, description="Mean in log10 space"),
    "sigma_log10": ToleranceConfig(rtol=1e-8, atol=1e-12, description="Sigma in log10 space"),
    
    # Distribution parameters
    "alpha": ToleranceConfig(rtol=1e-8, atol=1e-12, description="Gamma/Beta alpha parameter"),
    "beta": ToleranceConfig(rtol=1e-8, atol=1e-12, description="Gamma/Beta beta parameter"),
    
    # Epistemic uncertainty (looser due to sampling)
    "epistemic_mean": ToleranceConfig(rtol=1e-3, atol=1e-6, description="Mean across epistemic branches"),
    "epistemic_std": ToleranceConfig(rtol=5e-2, atol=1e-4, description="Std across epistemic branches"),
    
    # Default fallback
    "default": ToleranceConfig(rtol=1e-6, atol=1e-10, description="Default tolerance"),
}


def get_tolerance(value_type: str = "default") -> ToleranceConfig:
    """Get tolerance configuration for a value type."""
    return TOLERANCES.get(value_type, TOLERANCES["default"])


# ============================================================================
# Model Import Fixtures
# ============================================================================

@pytest.fixture(scope="session")
def fdhpy_models() -> Dict[str, Any]:
    """Import and return all fdhpy model classes."""
    models = {}
    try:
        from fdhpy import (
            YoungsEtAl2003,
            PetersenEtAl2011,
            MossEtAl2024,
            KuehnEtAl2024,
            LavrentiadisAbrahamson2023,
            ChiouEtAl2025,
        )
        models = {
            "YoungsEtAl2003": YoungsEtAl2003,
            "PetersenEtAl2011": PetersenEtAl2011,
            "MossEtAl2024": MossEtAl2024,
            "KuehnEtAl2024": KuehnEtAl2024,
            "LavrentiadisAbrahamson2023": LavrentiadisAbrahamson2023,
            "ChiouEtAl2025": ChiouEtAl2025,
        }
    except ImportError as e:
        pytest.skip(f"fdhpy not available: {e}")
    return models


@pytest.fixture(scope="session")
def pfdha_models() -> Dict[str, Any]:
    """Import and return all pfdha model classes."""
    models = {}
    try:
        from openquake.fdha.primary_surf_displ import (
            Youngs2003PrimaryFD,
            Petersen2011PrimaryFD,
            Moss2024PrimaryFD,
            Kuehn2024PrimaryFD,
            Lavrentiadis2023PrimaryFD,
            Chiou2025PrimaryFD,
        )
        models = {
            "Youngs2003PrimaryFD": Youngs2003PrimaryFD,
            "Petersen2011PrimaryFD": Petersen2011PrimaryFD,
            "Moss2024PrimaryFD": Moss2024PrimaryFD,
            "Kuehn2024PrimaryFD": Kuehn2024PrimaryFD,
            "Lavrentiadis2023PrimaryFD": Lavrentiadis2023PrimaryFD,
            "Chiou2025PrimaryFD": Chiou2025PrimaryFD,
        }
    except ImportError as e:
        pytest.fail(f"pfdha models import failed: {e}")
    return models


# ============================================================================
# Parameter Grid Configuration
# ============================================================================

@pytest.fixture(scope="session")
def parameter_grids() -> Dict[str, Any]:
    """Load parameter grid configuration from YAML."""
    config_path = CONFIG_DIR / "parameter_grids.yaml"
    if HAS_YAML and config_path.exists():
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    else:
        # Default configuration if YAML not found or not installed
        return get_default_parameter_grids()


def get_default_parameter_grids() -> Dict[str, Any]:
    """Return default parameter grids for testing."""
    return {
        "displacements": {
            "standard": [0.001, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0],
            "fine": list(np.logspace(-3, 1.5, 50)),
            "coarse": [0.01, 0.1, 1.0, 10.0],
        },
        "models": {
            "youngs2003": {
                "magnitudes": [6.0, 6.5, 7.0, 7.5, 8.0],
                "x_l_ratios": [0.1, 0.25, 0.5],
                "versions": ["d/ad"],
                "styles": ["all"],
            },
            "petersen2011": {
                "magnitudes": [6.5, 7.0, 7.5],
                "x_l_ratios": [0.1, 0.3, 0.5],
                "versions": ["elliptical", "quadratic"],
            },
            "moss2024": {
                "magnitudes": [6.5, 7.0, 7.5],
                "x_l_ratios": [0.25, 0.5],
                "versions": ["d/ad", "d/md"],
                "use_girs": [True, False],
                "complete": [True],
            },
            "kuehn2024": {
                "magnitudes": [7.0, 7.2, 7.6],
                "x_l_ratios": [0.3, 0.5, 0.7],
                "styles": ["normal", "reverse", "strike-slip"],
                "folded": [True],
                "epistemic": [False, True],
            },
            "lavrentiadis2023": {
                "magnitudes": [6.5, 7.0, 7.5],
                "x_l_ratios": [0.3, 0.5, 0.6],
                "metrics": ["aggregate"],
                "versions": ["full rupture"],
                "include_prob_zero": [True, False],
            },
            "chiou2025": {
                "magnitudes": [6.5, 7.0, 7.5],
                "x_l_ratios": [0.25, 0.5],
                "versions": ["model7", "model8.2"],
            },
        },
    }


# ============================================================================
# Standard Displacement Arrays
# ============================================================================

@pytest.fixture(scope="session")
def standard_displacements() -> np.ndarray:
    """Standard displacement array for testing."""
    return np.array([0.001, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0])


@pytest.fixture(scope="session")
def fine_displacements() -> np.ndarray:
    """Fine displacement array for detailed comparison."""
    return np.logspace(-3, 1.5, 50)


@pytest.fixture(scope="session")
def coarse_displacements() -> np.ndarray:
    """Coarse displacement array for quick checks."""
    return np.array([0.01, 0.1, 1.0, 10.0])


# ============================================================================
# Data Loading Fixtures
# ============================================================================

@pytest.fixture
def load_fdhpy_expected():
    """Factory fixture to load expected values from fdhpy tests."""
    def _load(filename: str) -> np.ndarray:
        filepath = FDHPY_EXPECTED_DIR / filename
        if not filepath.exists():
            pytest.skip(f"Expected file not found: {filepath}")
        return np.genfromtxt(
            filepath, delimiter=",", names=True, encoding="UTF-8-sig", dtype=None
        )
    return _load


@pytest.fixture
def load_pfdha_expected():
    """Factory fixture to load expected values from pfdha tests."""
    import pandas as pd
    def _load(filename: str):
        filepath = PFDHA_EXPECTED_DIR / filename
        if not filepath.exists():
            pytest.skip(f"Expected file not found: {filepath}")
        return pd.read_csv(filepath)
    return _load


# ============================================================================
# Comparison Result Tracking
# ============================================================================

@dataclass
class ComparisonResult:
    """Store detailed comparison results for reporting."""
    model: str
    test_case: str
    parameter_set: Dict[str, Any]
    value_type: str
    n_values: int
    max_abs_diff: float
    mean_abs_diff: float
    max_rel_diff: float
    mean_rel_diff: float
    passed: bool
    tolerance_used: ToleranceConfig
    failure_indices: Optional[np.ndarray] = None
    failure_message: Optional[str] = None


class ComparisonResultCollector:
    """Collect comparison results across tests for reporting."""
    
    def __init__(self):
        self.results = []
    
    def add(self, result: ComparisonResult):
        self.results.append(result)
    
    def get_summary(self) -> Dict[str, Any]:
        """Get summary statistics."""
        passed = sum(1 for r in self.results if r.passed)
        failed = sum(1 for r in self.results if not r.passed)
        return {
            "total": len(self.results),
            "passed": passed,
            "failed": failed,
            "pass_rate": passed / len(self.results) if self.results else 0,
            "by_model": self._group_by_model(),
        }
    
    def _group_by_model(self) -> Dict[str, Dict[str, int]]:
        """Group results by model."""
        groups = {}
        for r in self.results:
            if r.model not in groups:
                groups[r.model] = {"passed": 0, "failed": 0}
            if r.passed:
                groups[r.model]["passed"] += 1
            else:
                groups[r.model]["failed"] += 1
        return groups


@pytest.fixture(scope="session")
def result_collector():
    """Session-scoped result collector for generating summary reports."""
    return ComparisonResultCollector()


# ============================================================================
# pytest Hooks for Test Collection
# ============================================================================

def pytest_collection_modifyitems(config, items):
    """Modify test collection for ordering and filtering."""
    # Add slow marker to epistemic tests
    for item in items:
        if "epistemic" in item.keywords:
            item.add_marker(pytest.mark.slow)
