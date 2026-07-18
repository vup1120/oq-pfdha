"""
Aleatory uncertainty tests for Petersen et al. (2011).

STATUS: TESTABLE - pfdha exposes calc_params_* methods returning (mu, sd).

fdhpy exposes:
- stat_params_info["params"]["mu"] (natural log space, exp(mu) in cm)
- stat_params_info["params"]["sigma"] (natural log units)

pfdha Petersen2011PrimaryFD:
- calc_params_elliptical(mag, X_L_ratio) -> (mu, sd)
- calc_params_quadratic(mag, X_L_ratio) -> (mu, sd)
- calc_params_bilinear(mag, X_L_ratio) -> (mu, sd)

These are comparable after verifying units (both in ln(cm) space).
"""

import numpy as np
import pytest

# Import guards
fdhpy = pytest.importorskip("fdhpy", reason="fdhpy not installed - FDHI aleatory tests skipped")

from fdhpy import PetersenEtAl2011

try:
    from openquake.fdha.primary_surf_displ import Petersen2011PrimaryFD
except ImportError as e:
    pytest.skip(f"pfdha FD models not importable for aleatory tests: {e}", allow_module_level=True)


class TestPetersen2011Aleatory:
    """Aleatory uncertainty tests for Petersen et al. (2011)."""

    @pytest.mark.parametrize("magnitude", [6.5, 7.0, 7.5])
    @pytest.mark.parametrize("xl", [0.1, 0.3, 0.5])
    def test_mean_mu_elliptical(self, magnitude, xl):
        """Compare mean mu (in ln(cm) space) for elliptical model."""
        # fdhpy
        fdhpy_model = PetersenEtAl2011(
            magnitude=magnitude,
            xl=xl,
            version="elliptical",
            displ_array=np.array([0.1]),  # dummy, needed for initialization
        )
        stat_params = fdhpy_model.stat_params_info
        fdhpy_mu = stat_params["params"]["mu"]

        # pfdha
        pfdha_model = Petersen2011PrimaryFD()
        mu, sd = pfdha_model.calc_params_elliptical(mag=magnitude, X_L_ratio=np.array([xl]))
        pfdha_mu = float(mu[0])

        np.testing.assert_allclose(
            pfdha_mu, fdhpy_mu,
            rtol=1e-6, atol=1e-10,
            err_msg=f"Petersen 2011 elliptical mean: M={magnitude}, x/L={xl}"
        )

    @pytest.mark.parametrize("magnitude", [6.5, 7.0, 7.5])
    @pytest.mark.parametrize("xl", [0.1, 0.3, 0.5])
    def test_sigma_elliptical(self, magnitude, xl):
        """Compare sigma (in ln units) for elliptical model."""
        # fdhpy
        fdhpy_model = PetersenEtAl2011(
            magnitude=magnitude,
            xl=xl,
            version="elliptical",
            displ_array=np.array([0.1]),
        )
        stat_params = fdhpy_model.stat_params_info
        fdhpy_sigma = stat_params["params"]["sigma"]

        # pfdha
        pfdha_model = Petersen2011PrimaryFD()
        mu, sd = pfdha_model.calc_params_elliptical(mag=magnitude, X_L_ratio=np.array([xl]))
        pfdha_sigma = float(sd[0])

        np.testing.assert_allclose(
            pfdha_sigma, fdhpy_sigma,
            rtol=1e-6, atol=1e-10,
            err_msg=f"Petersen 2011 elliptical sigma: M={magnitude}, x/L={xl}"
        )

    @pytest.mark.parametrize("magnitude", [6.5, 7.0, 7.5])
    @pytest.mark.parametrize("xl", [0.1, 0.3, 0.5])
    def test_mean_mu_quadratic(self, magnitude, xl):
        """Compare mean mu (in ln(cm) space) for quadratic model."""
        # fdhpy
        fdhpy_model = PetersenEtAl2011(
            magnitude=magnitude,
            xl=xl,
            version="quadratic",
            displ_array=np.array([0.1]),
        )
        stat_params = fdhpy_model.stat_params_info
        fdhpy_mu = stat_params["params"]["mu"]

        # pfdha
        pfdha_model = Petersen2011PrimaryFD()
        mu, sd = pfdha_model.calc_params_quadratic(mag=magnitude, X_L_ratio=np.array([xl]))
        pfdha_mu = float(mu[0])

        np.testing.assert_allclose(
            pfdha_mu, fdhpy_mu,
            rtol=1e-6, atol=1e-10,
            err_msg=f"Petersen 2011 quadratic mean: M={magnitude}, x/L={xl}"
        )

    @pytest.mark.parametrize("magnitude", [6.5, 7.0, 7.5])
    @pytest.mark.parametrize("xl", [0.1, 0.3, 0.5])
    def test_sigma_quadratic(self, magnitude, xl):
        """Compare sigma (in ln units) for quadratic model."""
        # fdhpy
        fdhpy_model = PetersenEtAl2011(
            magnitude=magnitude,
            xl=xl,
            version="quadratic",
            displ_array=np.array([0.1]),
        )
        stat_params = fdhpy_model.stat_params_info
        fdhpy_sigma = stat_params["params"]["sigma"]

        # pfdha
        pfdha_model = Petersen2011PrimaryFD()
        mu, sd = pfdha_model.calc_params_quadratic(mag=magnitude, X_L_ratio=np.array([xl]))
        pfdha_sigma = float(sd[0])

        np.testing.assert_allclose(
            pfdha_sigma, fdhpy_sigma,
            rtol=1e-6, atol=1e-10,
            err_msg=f"Petersen 2011 quadratic sigma: M={magnitude}, x/L={xl}"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])

















