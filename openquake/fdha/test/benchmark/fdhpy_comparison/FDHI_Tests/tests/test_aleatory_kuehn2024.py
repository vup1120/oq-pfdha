"""
Aleatory uncertainty tests for Kuehn et al. (2024).

STATUS: TESTABLE (using internal API) - pfdha exposes _calc_params() returning mu and sigma.

fdhpy exposes:
- sd_med property (std dev of predicted median)
- sd_sigma_T property (std dev of predicted sigma)
- _calc_sigma_mag() (magnitude aleatory variability)
- _calc_sigma_xl() -> (sigma_xl_u1, sigma_xl_u2)
- _calc_mean_mu(u_star) (mean in transformed space)

pfdha Kuehn2024PrimaryFD:
- _calc_params(coeffs, mag, X_L_ratio, style) -> (model_id, lam, mu, std_total, std_within, std_mode)

Note: This relies on internal API (_calc_params), so tests may break if internal API changes.
"""

import numpy as np
import pytest
import pandas as pd

# Import guards
fdhpy = pytest.importorskip("fdhpy", reason="fdhpy not installed – FDHI aleatory tests skipped")

from fdhpy import KuehnEtAl2024

try:
    from openquake.fdha.primary_surf_displ import Kuehn2024PrimaryFD
    from openquake.fdha.primary_surf_displ.kuehn2024.load_data import DATA as DATA_COEFFICIENTS
except ImportError as e:
    pytest.skip(f"pfdha FD models not importable for aleatory tests: {e}", allow_module_level=True)


def _get_pfdha_median_coeffs(style):
    """Extract median coefficients for a given style from pfdha's DATA_COEFFICIENTS."""
    mean_coeffs_data = DATA_COEFFICIENTS[style]['mean']
    if isinstance(mean_coeffs_data, pd.DataFrame):
        if 'median' in mean_coeffs_data.index:
            return mean_coeffs_data.loc['median']
        else:
            return mean_coeffs_data.iloc[0]
    return mean_coeffs_data


class TestKuehn2024Aleatory:
    """
    Aleatory uncertainty tests for Kuehn et al. (2024).
    
    Note: These tests use pfdha's internal _calc_params() method.
    """

    @pytest.mark.parametrize("magnitude", [7.0, 7.5])
    @pytest.mark.parametrize("xl", [0.3, 0.5])
    @pytest.mark.parametrize("style", ["normal", "reverse", "strike-slip"])
    def test_total_sigma(self, magnitude, xl, style):
        """
        Compare total sigma (combined magnitude and x/L aleatory).
        
        fdhpy: sqrt(sigma_mag^2 + sigma_xl_u1^2)
        pfdha: std_total from _calc_params()
        """
        # fdhpy
        fdhpy_model = KuehnEtAl2024(
            style=style,
            magnitude=magnitude,
            xl=xl,
            version="median_coeffs",
            folded=True,
            displ_array=np.array([0.1]),
        )
        
        # Trigger internal calculations
        _ = fdhpy_model.stat_params_info
        fdhpy_sigma_mag = fdhpy_model._calc_sigma_mag()
        fdhpy_sigma_xl_u1, _ = fdhpy_model._calc_sigma_xl()
        fdhpy_total_sigma = np.sqrt(fdhpy_sigma_mag**2 + fdhpy_sigma_xl_u1**2)

        # pfdha
        pfdha_model = Kuehn2024PrimaryFD()
        coeffs = _get_pfdha_median_coeffs(style)
        _, lam, mu, std_total, std_within, std_mode = pfdha_model._calc_params(
            coeffs, magnitude, np.array([xl]), style
        )
        pfdha_total_sigma = float(std_total[0]) if hasattr(std_total, '__len__') else float(std_total)

        np.testing.assert_allclose(
            pfdha_total_sigma, float(fdhpy_total_sigma),
            rtol=1e-6, atol=1e-10,
            err_msg=f"Kuehn 2024 total_sigma: {style}, M={magnitude}, x/L={xl}"
        )

    @pytest.mark.parametrize("magnitude", [7.0, 7.5])
    @pytest.mark.parametrize("xl", [0.3, 0.5])
    @pytest.mark.parametrize("style", ["normal", "reverse", "strike-slip"])
    def test_mean_mu(self, magnitude, xl, style):
        """
        Compare mean mu (in transformed space).
        
        fdhpy: _calc_mean_mu(u_star)
        pfdha: mu from _calc_params()
        """
        # fdhpy
        fdhpy_model = KuehnEtAl2024(
            style=style,
            magnitude=magnitude,
            xl=xl,
            version="median_coeffs",
            folded=True,
            displ_array=np.array([0.1]),
        )
        # Note: fdhpy uses folded xl for u_star
        fdhpy_mu = fdhpy_model._calc_mean_mu(xl)

        # pfdha
        pfdha_model = Kuehn2024PrimaryFD()
        coeffs = _get_pfdha_median_coeffs(style)
        _, lam, mu, std_total, std_within, std_mode = pfdha_model._calc_params(
            coeffs, magnitude, np.array([xl]), style
        )
        pfdha_mu = float(mu[0]) if hasattr(mu, '__len__') else float(mu)

        np.testing.assert_allclose(
            pfdha_mu, float(fdhpy_mu),
            rtol=1e-6, atol=1e-10,
            err_msg=f"Kuehn 2024 mean_mu: {style}, M={magnitude}, x/L={xl}"
        )

    @pytest.mark.parametrize("magnitude", [7.0, 7.5])
    @pytest.mark.parametrize("style", ["normal", "reverse", "strike-slip"])
    def test_sigma_mag_component(self, magnitude, style):
        """
        Compare magnitude-related sigma component.
        
        fdhpy: _calc_sigma_mag()
        pfdha: std_mode from _calc_params()
        """
        xl = 0.5  # Fixed for this test

        # fdhpy
        fdhpy_model = KuehnEtAl2024(
            style=style,
            magnitude=magnitude,
            xl=xl,
            version="median_coeffs",
            folded=True,
            displ_array=np.array([0.1]),
        )
        _ = fdhpy_model.stat_params_info
        fdhpy_sigma_mag = fdhpy_model._calc_sigma_mag()

        # pfdha
        pfdha_model = Kuehn2024PrimaryFD()
        coeffs = _get_pfdha_median_coeffs(style)
        _, lam, mu, std_total, std_within, std_mode = pfdha_model._calc_params(
            coeffs, magnitude, np.array([xl]), style
        )
        pfdha_sigma_mag = float(std_mode) if not hasattr(std_mode, '__len__') else float(std_mode[0])

        np.testing.assert_allclose(
            pfdha_sigma_mag, float(fdhpy_sigma_mag),
            rtol=1e-6, atol=1e-10,
            err_msg=f"Kuehn 2024 sigma_mag: {style}, M={magnitude}"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])

















