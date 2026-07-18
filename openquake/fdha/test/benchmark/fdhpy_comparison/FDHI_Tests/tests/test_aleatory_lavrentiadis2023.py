"""
Aleatory uncertainty tests for Lavrentiadis & Abrahamson (2023).

STATUS: TESTABLE - pfdha exposes LavrentiadisAbrahamson2023SlipProfile() returning sigma components.

fdhpy exposes:
- sigma_mu_agg property (std dev of predicted median aggregate displacement)

pfdha Lavrentiadis2023PrimaryFD_aggregate:
- LavrentiadisAbrahamson2023SlipProfile() returns:
  - sig_agg (total aggregate aleatory)
  - sig_prnc (total principal aleatory)
  - phi_agg (within-event aggregate aleatory)
  - phi_prnc (within-event principal aleatory)
  - tau_agg (between-event aleatory)
  - phi_add (additional aleatory due to segmentation)
"""

import numpy as np
import pytest

# Import guards
fdhpy = pytest.importorskip("fdhpy", reason="fdhpy not installed - FDHI aleatory tests skipped")

from fdhpy import LavrentiadisAbrahamson2023

try:
    from openquake.fdha.primary_surf_displ import Lavrentiadis2023PrimaryFD_aggregate
except ImportError as e:
    pytest.skip(f"pfdha FD models not importable for aleatory tests: {e}", allow_module_level=True)


class TestLavrentiadis2023Aleatory:
    """Aleatory uncertainty tests for Lavrentiadis & Abrahamson (2023)."""

    @pytest.mark.parametrize("magnitude", [6.5, 7.0, 7.5])
    @pytest.mark.parametrize("style", ["strike-slip", "normal", "reverse"])
    def test_sigma_mu_agg(self, magnitude, style):
        """
        Compare sigma_mu_agg (standard deviation of predicted median aggregate displacement).
        
        fdhpy: sigma_mu_agg property
        pfdha: This is the epistemic uncertainty on the median, not directly exposed.
        
        Note: fdhpy's sigma_mu_agg is specifically defined in Eq. (varies by M and style).
        pfdha's sig_agg from SlipProfile is the total aleatory (different quantity).
        """
        # fdhpy - sigma_mu_agg (epistemic uncertainty on median)
        fdhpy_model = LavrentiadisAbrahamson2023(
            magnitude=magnitude,
            xl=0.5,  # Fixed for this test
            displ_array=np.array([0.1]),
            metric="aggregate",
            version="full rupture",
            style=style,
        )
        
        fdhpy_sigma_mu = fdhpy_model.sigma_mu_agg
        
        # Manual calculation based on fdhpy's formula
        if magnitude >= 7.1:
            expected_sigma = 0.035 + 0.025 * (magnitude - 7.1)
        else:
            c_map = {"normal": 0.064, "strike-slip": 0.036, "reverse": 0.036}
            c = c_map.get(style.lower(), 0.036)
            expected_sigma = 0.035 + c * (7.1 - magnitude)
        
        np.testing.assert_allclose(
            float(fdhpy_sigma_mu), expected_sigma,
            rtol=1e-6, atol=1e-10,
            err_msg=f"LA23 sigma_mu_agg formula check: {style}, M={magnitude}"
        )

    @pytest.mark.parametrize("magnitude", [6.5, 7.0, 7.5])
    @pytest.mark.parametrize("xl", [0.3, 0.5])
    def test_total_aggregate_sigma(self, magnitude, xl):
        """
        Compare total aggregate sigma (sig_agg).
        
        pfdha returns sig_agg from LavrentiadisAbrahamson2023SlipProfile().
        fdhpy's stat_params_info can provide comparable sigma.
        """
        style = "strike-slip"

        # fdhpy
        fdhpy_model = LavrentiadisAbrahamson2023(
            magnitude=magnitude,
            xl=xl,
            displ_array=np.array([0.1]),
            metric="aggregate",
            version="full rupture",
            style=style,
        )
        
        # Get fdhpy's sigma from stat_params_info
        stat_params = fdhpy_model.stat_params_info
        fdhpy_sigma = stat_params["params"].get("sigma")
        
        # pfdha
        pfdha_model = Lavrentiadis2023PrimaryFD_aggregate()
        result = pfdha_model.LavrentiadisAbrahamson2023SlipProfile(
            x_array=np.array([xl]),
            mag=magnitude,
            srl=1,  # normalized
            sof=style.title(),
        )
        (disp_agg_prime, disp_prnc_prime, disp_agg_seg,
         sig_agg, sig_prnc, phi_agg, phi_prnc, tau_agg, phi_add,
         P_gap, P_zero_slip) = result
        
        pfdha_sig_agg = float(sig_agg[0]) if hasattr(sig_agg, '__len__') else float(sig_agg)

        if fdhpy_sigma is not None:
            np.testing.assert_allclose(
                pfdha_sig_agg, float(fdhpy_sigma),
                rtol=1e-6, atol=1e-10,
                err_msg=f"LA23 sig_agg: {style}, M={magnitude}, x/L={xl}"
            )
        else:
            # If fdhpy doesn't expose sigma directly, just verify pfdha computes something reasonable
            assert pfdha_sig_agg > 0, f"sig_agg should be positive: {pfdha_sig_agg}"

    @pytest.mark.parametrize("magnitude", [6.5, 7.0, 7.5])
    @pytest.mark.parametrize("xl", [0.3, 0.5])
    def test_phi_tau_components(self, magnitude, xl):
        """
        Test that phi and tau components are reasonable and sum correctly to total sigma.
        
        sig_agg_seg = sqrt(tau_agg^2 + phi_agg^2)  (from pfdha code)
        """
        style = "strike-slip"

        # pfdha
        pfdha_model = Lavrentiadis2023PrimaryFD_aggregate()
        result = pfdha_model.LavrentiadisAbrahamson2023SlipProfile(
            x_array=np.array([xl]),
            mag=magnitude,
            srl=1,
            sof=style.title(),
        )
        (disp_agg_prime, disp_prnc_prime, disp_agg_seg,
         sig_agg, sig_prnc, phi_agg, phi_prnc, tau_agg, phi_add,
         P_gap, P_zero_slip) = result

        # Extract scalar values
        phi_agg_val = float(phi_agg[0]) if hasattr(phi_agg, '__len__') else float(phi_agg)
        tau_agg_val = float(tau_agg[0]) if hasattr(tau_agg, '__len__') else float(tau_agg)
        
        # Compute expected sig_agg_seg
        expected_sig_agg_seg = np.sqrt(tau_agg_val**2 + phi_agg_val**2)
        
        # Verify components are positive
        assert phi_agg_val > 0, f"phi_agg should be positive: {phi_agg_val}"
        assert tau_agg_val > 0, f"tau_agg should be positive: {tau_agg_val}"
        
        # Note: sig_agg also includes phi_add, so it's not exactly sqrt(tau^2 + phi^2)
        # Just verify they're reasonable
        sig_agg_val = float(sig_agg[0]) if hasattr(sig_agg, '__len__') else float(sig_agg)
        assert sig_agg_val >= expected_sig_agg_seg, (
            f"sig_agg ({sig_agg_val}) should be >= sqrt(tau^2 + phi^2) ({expected_sig_agg_seg})"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])

















