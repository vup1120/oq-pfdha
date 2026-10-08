"""
Simplified PFDHA Reference Tests

Simple pytest tests that directly compare fdhpy and pfdha without
complex helper classes. Easier to debug and understand.

These tests require both `fdhpy` and `pfdha` to be installed:
  pip install fdhpy
  pip install -e /path/to/pfdha
"""

import numpy as np
import pytest

# Require fdhpy to be installed; skip entire module if not available
fdhpy = pytest.importorskip("fdhpy", reason="fdhpy not installed - FDHI reference tests skipped")

from fdhpy import (
    YoungsEtAl2003,
    PetersenEtAl2011,
    MossEtAl2024,
    KuehnEtAl2024,
    LavrentiadisAbrahamson2023,
    ChiouEtAl2025,
)

# Import pfdha models; fail the module if pfdha is not importable
try:
    from openquake.fdha.primary_surf_displ import (
        Youngs2003PrimaryFD,
        Petersen2011PrimaryFD,
        Moss2024PrimaryFD,
        Kuehn2024PrimaryFD,
        Lavrentiadis2023PrimaryFD_aggregate,
        Lavrentiadis2023PrimaryFD_principal,
        Chiou2025PrimaryFD,
    )
except ImportError as e:
    pytest.skip(f"pfdha OpenQuake FD models not importable: {e}", allow_module_level=True)

# Standard test parameters
DISPLACEMENTS = np.array([0.001, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0])


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture(scope="module")
def fdhpy_models():
    """Import fdhpy models."""
    return {
        "YoungsEtAl2003": YoungsEtAl2003,
        "PetersenEtAl2011": PetersenEtAl2011,
        "MossEtAl2024": MossEtAl2024,
        "KuehnEtAl2024": KuehnEtAl2024,
        "LavrentiadisAbrahamson2023": LavrentiadisAbrahamson2023,
        "ChiouEtAl2025": ChiouEtAl2025,
    }


@pytest.fixture(scope="module")
def pfdha_models():
    """Import pfdha models."""
    return {
        "Youngs2003PrimaryFD": Youngs2003PrimaryFD,
        "Petersen2011PrimaryFD": Petersen2011PrimaryFD,
        "Moss2024PrimaryFD": Moss2024PrimaryFD,
        "Kuehn2024PrimaryFD": Kuehn2024PrimaryFD,
        "Lavrentiadis2023PrimaryFD_aggregate": Lavrentiadis2023PrimaryFD_aggregate,
        "Lavrentiadis2023PrimaryFD_principal": Lavrentiadis2023PrimaryFD_principal,
        "Chiou2025PrimaryFD": Chiou2025PrimaryFD,
    }


# ============================================================================
# Youngs 2003 Tests
# ============================================================================

class TestYoungs2003:
    """Youngs et al. (2003) reference tests."""
    
    @pytest.mark.parametrize("magnitude", [6.0, 6.5, 7.0, 7.5])
    @pytest.mark.parametrize("xl", [0.1, 0.25, 0.5])
    def test_prob_exceed(self, magnitude, xl, fdhpy_models, pfdha_models):
        """Test exceedance probability matches fdhpy."""
        # fdhpy
        fdhpy_model = fdhpy_models["YoungsEtAl2003"](
            magnitude=magnitude,
            xl=xl,
            version="d/ad",
            displ_array=DISPLACEMENTS,
        )
        fdhpy_result = fdhpy_model.prob_exceed
        
        # pfdha
        pfdha_model = pfdha_models["Youngs2003PrimaryFD"]()
        pfdha_result = pfdha_model.get_prob(
            d=DISPLACEMENTS,
            X_L_ratio=np.array([xl]),
            mag=magnitude,
            style="all",
            norm_disp_type="AD"
        ).flatten()
        
        # Compare (relaxed tolerance for integration differences)
        np.testing.assert_allclose(
            fdhpy_result, pfdha_result,
            rtol=1e-4, atol=1e-8,
            err_msg=f"Youngs 2003: M={magnitude}, x/L={xl}"
        )


# ============================================================================
# Petersen 2011 Tests
# ============================================================================

class TestPetersen2011:
    """Petersen et al. (2011) reference tests."""
    
    @pytest.mark.parametrize("magnitude", [6.5, 7.0, 7.5])
    @pytest.mark.parametrize("xl", [0.1, 0.3, 0.5])
    @pytest.mark.parametrize("version", ["elliptical", "quadratic"])
    def test_prob_exceed(self, magnitude, xl, version, fdhpy_models, pfdha_models):
        """Test exceedance probability matches fdhpy."""
        # fdhpy
        fdhpy_model = fdhpy_models["PetersenEtAl2011"](
            magnitude=magnitude,
            xl=xl,
            version=version,
            displ_array=DISPLACEMENTS,
        )
        fdhpy_result = fdhpy_model.prob_exceed
        
        # pfdha
        pfdha_model = pfdha_models["Petersen2011PrimaryFD"]()
        pfdha_result = pfdha_model.get_prob(
            d=DISPLACEMENTS,
            X_L_ratio=np.array([xl]),
            mag=magnitude,
            version=version
        ).flatten()
        
        # Compare
        np.testing.assert_allclose(
            fdhpy_result, pfdha_result,
            rtol=1e-6, atol=1e-10,
            err_msg=f"Petersen 2011 {version}: M={magnitude}, x/L={xl}"
        )


# ============================================================================
# Moss 2024 Tests
# ============================================================================

class TestMoss2024:
    """Moss et al. (2024) reference tests."""
    
    @pytest.mark.parametrize("magnitude", [6.5, 7.0, 7.5])
    @pytest.mark.parametrize("xl", [0.25, 0.5])
    @pytest.mark.parametrize("version,norm_type", [("d/ad", "AD"), ("d/md", "MD")])
    @pytest.mark.parametrize("use_girs", [True, False])
    def test_prob_exceed(self, magnitude, xl, version, norm_type, use_girs, fdhpy_models, pfdha_models):
        """Test exceedance probability matches fdhpy."""
        # fdhpy
        fdhpy_model = fdhpy_models["MossEtAl2024"](
            magnitude=magnitude,
            xl=xl,
            version=version,
            displ_array=DISPLACEMENTS,
            use_girs=use_girs,
            complete=True,
        )
        fdhpy_result = fdhpy_model.prob_exceed
        
        # pfdha
        source = "GIRS" if use_girs else "EQS"
        pfdha_model = pfdha_models["Moss2024PrimaryFD"]()
        pfdha_result = pfdha_model.get_prob(
            d=DISPLACEMENTS,
            X_L_ratio=np.array([xl]),
            mag=magnitude,
            version=norm_type,
            source=source,
            completeness="complete"
        ).flatten()
        
        # Compare
        np.testing.assert_allclose(
            fdhpy_result, pfdha_result,
            rtol=1e-6, atol=1e-10,
            err_msg=f"Moss 2024 {version} {source}: M={magnitude}, x/L={xl}"
        )


# ============================================================================
# Kuehn 2024 Tests
# ============================================================================

class TestKuehn2024:
    """Kuehn et al. (2024) reference tests."""
    
    @pytest.mark.parametrize("magnitude", [7.0, 7.2, 7.6])
    @pytest.mark.parametrize("xl", [0.3, 0.5, 0.7])
    @pytest.mark.parametrize("style", ["normal", "reverse", "strike-slip"])
    def test_prob_exceed_median(self, magnitude, xl, style, fdhpy_models, pfdha_models):
        """Test exceedance probability with median coefficients."""
        # fdhpy
        fdhpy_model = fdhpy_models["KuehnEtAl2024"](
            style=style,
            magnitude=magnitude,
            xl=xl,
            version="median_coeffs",
            folded=True,
            displ_array=DISPLACEMENTS,
        )
        fdhpy_result = fdhpy_model.prob_exceed
        
        # pfdha
        pfdha_model = pfdha_models["Kuehn2024PrimaryFD"]()
        pfdha_result = pfdha_model.get_prob(
            d=DISPLACEMENTS,
            X_L_ratio=np.array([xl]),
            mag=magnitude,
            style=style,
            folded=True,
            epistemic_uncertainty=False
        ).flatten()
        
        # Compare
        np.testing.assert_allclose(
            fdhpy_result, pfdha_result,
            rtol=1e-6, atol=1e-10,
            err_msg=f"Kuehn 2024 {style}: M={magnitude}, x/L={xl}"
        )


# ============================================================================
# Lavrentiadis 2023 Tests
# ============================================================================

class TestLavrentiadis2023:
    """Lavrentiadis & Abrahamson (2023) reference tests."""
    
    @pytest.mark.parametrize("magnitude", [6.5, 7.0, 7.5])
    @pytest.mark.parametrize("xl", [0.3, 0.5, 0.6])
    @pytest.mark.parametrize("include_prob_zero", [True, False])
    @pytest.mark.parametrize("style", ["strike-slip", "normal", "reverse"])
    def test_prob_exceed_aggregate(self, magnitude, xl, include_prob_zero, style, fdhpy_models, pfdha_models):
        """Test exceedance probability for aggregate metric."""
        # fdhpy - requires style parameter
        fdhpy_model = fdhpy_models["LavrentiadisAbrahamson2023"](
            magnitude=magnitude,
            xl=xl,
            displ_array=DISPLACEMENTS,
            metric="aggregate",
            version="full rupture",
            style=style,  # Required parameter
            include_prob_zero=include_prob_zero,
        )
        fdhpy_result = fdhpy_model.prob_exceed
        
        # pfdha
        pfdha_model = pfdha_models["Lavrentiadis2023PrimaryFD_aggregate"]()
        pfdha_result = pfdha_model.get_prob(
            d=DISPLACEMENTS,
            X_L_ratio=np.array([xl]),
            mag=magnitude,
            style=style,
            output_type="disp_agg_prime",
            include_zero_slip=include_prob_zero,
        ).flatten()
        
        # Compare
        np.testing.assert_allclose(
            fdhpy_result, pfdha_result,
            rtol=1e-6, atol=1e-10,
            err_msg=f"LA23 aggregate: {style}, M={magnitude}, x/L={xl}, zero={include_prob_zero}"
        )

    @pytest.mark.parametrize("magnitude", [6.5, 7.0, 7.5])
    @pytest.mark.parametrize("xl", [0.3, 0.5, 0.6])
    @pytest.mark.parametrize("include_prob_zero", [True, False])
    @pytest.mark.parametrize("style", ["strike-slip", "normal", "reverse"])
    def test_prob_exceed_sum_of_principal(self, magnitude, xl, include_prob_zero, style, fdhpy_models, pfdha_models):
        """Sum-of-principal class (disp_prnc_prime) against the fdhpy
        'sum-of-principal' metric, full-rupture version; with the zero terms
        both models scale by (1 - P_gap)(1 - P_zero_slip) (paper Eq. 31)."""
        fdhpy_model = fdhpy_models["LavrentiadisAbrahamson2023"](
            magnitude=magnitude,
            xl=xl,
            displ_array=DISPLACEMENTS,
            metric="sum-of-principal",
            version="full rupture",
            style=style,
            include_prob_zero=include_prob_zero,
        )
        pfdha_result = pfdha_models["Lavrentiadis2023PrimaryFD_principal"]().get_prob(
            d=DISPLACEMENTS,
            X_L_ratio=np.array([xl]),
            mag=magnitude,
            style=style,
            include_zero_slip=include_prob_zero,
        ).flatten()
        np.testing.assert_allclose(
            fdhpy_model.prob_exceed, pfdha_result,
            rtol=1e-6, atol=1e-10,
            err_msg=f"LA23 sum-of-principal: {style}, M={magnitude}, x/L={xl}, zero={include_prob_zero}"
        )

    @pytest.mark.parametrize("magnitude", [6.5, 7.0, 7.5])
    @pytest.mark.parametrize("xl", [0.3, 0.5, 0.6])
    @pytest.mark.parametrize("style", ["strike-slip", "normal", "reverse"])
    def test_prob_exceed_aggregate_individual_segment(self, magnitude, xl, style, fdhpy_models, pfdha_models):
        """Aggregate class, single-segment variant (disp_agg_seg) against the
        fdhpy 'individual segment' version, without zero terms."""
        fdhpy_model = fdhpy_models["LavrentiadisAbrahamson2023"](
            magnitude=magnitude,
            xl=xl,
            displ_array=DISPLACEMENTS,
            metric="aggregate",
            version="individual segment",
            style=style,
            include_prob_zero=False,
        )
        pfdha_result = pfdha_models["Lavrentiadis2023PrimaryFD_aggregate"]().get_prob(
            d=DISPLACEMENTS,
            X_L_ratio=np.array([xl]),
            mag=magnitude,
            style=style,
            output_type="disp_agg_seg",
            include_zero_slip=False,
        ).flatten()
        np.testing.assert_allclose(
            fdhpy_model.prob_exceed, pfdha_result,
            rtol=1e-6, atol=1e-10,
            err_msg=f"LA23 aggregate segment: {style}, M={magnitude}, x/L={xl}"
        )

    @pytest.mark.parametrize("magnitude", [6.5, 7.0, 7.5])
    @pytest.mark.parametrize("xl", [0.3, 0.5, 0.6])
    @pytest.mark.parametrize("style", ["strike-slip", "normal", "reverse"])
    def test_prob_exceed_sum_of_principal_individual_segment(self, magnitude, xl, style, fdhpy_models, pfdha_models):
        """Sum-of-principal class, single-segment variant (disp_prnc_seg)
        against the fdhpy 'individual segment' version, without zero terms
        (with them pfdha applies 1 - P(D_P = 0), decision B2)."""
        fdhpy_model = fdhpy_models["LavrentiadisAbrahamson2023"](
            magnitude=magnitude,
            xl=xl,
            displ_array=DISPLACEMENTS,
            metric="sum-of-principal",
            version="individual segment",
            style=style,
            include_prob_zero=False,
        )
        pfdha_result = pfdha_models["Lavrentiadis2023PrimaryFD_principal"]().get_prob(
            d=DISPLACEMENTS,
            X_L_ratio=np.array([xl]),
            mag=magnitude,
            style=style,
            output_type="disp_prnc_seg",
            include_zero_slip=False,
        ).flatten()
        np.testing.assert_allclose(
            fdhpy_model.prob_exceed, pfdha_result,
            rtol=1e-6, atol=1e-10,
            err_msg=f"LA23 sum-of-principal segment: {style}, M={magnitude}, x/L={xl}"
        )

    def test_prob_exceed_aggregate_individual_segment_with_zero_terms(self, fdhpy_models, pfdha_models):
        """A single segment has no gap (p. 21) and P(D_P = 0) is the
        zero-PRINCIPAL probability (Eqs. 31-32), so neither fdhpy nor pfdha
        applies a zero term to the single-segment aggregate; the flag is
        ignored with a warning."""
        fdhpy_model = fdhpy_models["LavrentiadisAbrahamson2023"](
            magnitude=6.5,
            xl=0.3,
            displ_array=DISPLACEMENTS,
            metric="aggregate",
            version="individual segment",
            style="normal",
            include_prob_zero=True,
        )
        pfdha_result = pfdha_models["Lavrentiadis2023PrimaryFD_aggregate"]().get_prob(
            d=DISPLACEMENTS,
            X_L_ratio=np.array([0.3]),
            mag=6.5,
            style="normal",
            output_type="disp_agg_seg",
            include_zero_slip=True,
        ).flatten()
        np.testing.assert_allclose(
            fdhpy_model.prob_exceed, pfdha_result, rtol=1e-6, atol=1e-10)


# ============================================================================
# Chiou 2025 Tests
# ============================================================================

class TestChiou2025:
    """Chiou et al. (2025) reference tests."""
    
    @pytest.mark.parametrize("magnitude", [6.5, 7.0, 7.5])
    @pytest.mark.parametrize("xl", [0.25, 0.5])
    @pytest.mark.parametrize("version", ["model7", "model8.2"])
    def test_prob_exceed(self, magnitude, xl, version, fdhpy_models, pfdha_models):
        """Test exceedance probability matches fdhpy."""
        # fdhpy
        fdhpy_model = fdhpy_models["ChiouEtAl2025"](
            magnitude=magnitude,
            xl=xl,
            version=version,
            displ_array=DISPLACEMENTS,
        )
        fdhpy_result = fdhpy_model.prob_exceed
        
        # pfdha
        pfdha_model = pfdha_models["Chiou2025PrimaryFD"]()
        pfdha_result = pfdha_model.get_prob(
            d=DISPLACEMENTS,
            X_L_ratio=np.array([xl]),
            mag=magnitude,
            version=version
        ).flatten()
        
        # Compare
        np.testing.assert_allclose(
            fdhpy_result, pfdha_result,
            rtol=1e-6, atol=1e-10,
            err_msg=f"Chiou 2025 {version}: M={magnitude}, x/L={xl}"
        )


# ============================================================================
# Run standalone
# ============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
