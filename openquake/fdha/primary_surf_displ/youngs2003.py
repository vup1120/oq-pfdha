"""
Principal fault-displacement model of Youngs et al. (2003), using the
Wells and Coppersmith (1994) magnitude-displacement scaling relations.

References
----------
Youngs, R.R., et al. (2003). A methodology for probabilistic fault
displacement hazard analysis (PFDHA). Earthquake Spectra, 19(1), 191-219.
"""

import numpy as np
from scipy.stats import gamma, norm, beta
from openquake.fdha.primary_surf_displ.base import BasePrimarySurfDispl, BaseSecondarySurfDispl

class Youngs2003PrimaryFD(BasePrimarySurfDispl):
    """
    Model of Youngs et al. (2003) for the probability of exceeding
    principal surface fault displacement thresholds on normal faults.
    
    This model uses Wells & Coppersmith (1994) magnitude-displacement scaling
    relations. Two coefficient sets are available:
    
    - style="all": Uses "All styles" coefficients from WC94 
    - style="normal": Uses "Normal faulting" specific coefficients from WC94
    
    Reference:
        Youngs, R. R., et al. (2003). A methodology for probabilistic fault
        displacement hazard analysis (PFDHA). Earthquake Spectra, 19(1), 191-219.
    """
    
    # Wells & Coppersmith (1994) coefficients for "All styles"
    # (recommended - consistent with Youngs et al. 2003 paper and fdhpy)
    _WC94_ALL = {
        "AD": {"intercept": -4.80, "slope": 0.69, "sigma": 0.36},
        "MD": {"intercept": -5.46, "slope": 0.82, "sigma": 0.42},
    }
    
    # Wells & Coppersmith (1994) coefficients for "Normal faulting"
    _WC94_NORMAL = {
        "AD": {"intercept": -4.45, "slope": 0.63, "sigma": 0.33},
        "MD": {"intercept": -5.90, "slope": 0.89, "sigma": 0.38},
    }
    
    # Integration parameters
    _DZ = 0.1   # Step size in epsilon space

    _ACCEPTED_DISP_TYPES = frozenset(["AD", "MD"])
    _ACCEPTED_STYLES = frozenset(["all", "normal"])

    def __init__(self, n_sigma=6.0):
        """
        :param n_sigma:
            Half-width of the ±σ ε-space integration truncation. Defaults to 6
            (improves accuracy over the historical ±3σ). Overridable from the
            logic tree via ``[Youngs2003PrimaryFD] n_sigma = <value>``.
        """
        super().__init__()
        self._N_EPS = float(n_sigma)  # ±n_sigma truncation in epsilon space
        if self._N_EPS <= 0.0:
            raise ValueError(f"n_sigma must be positive; got {self._N_EPS}")
    
    def _get_wc94_coeffs(self, style, norm_disp_type):
        """
        Get Wells & Coppersmith (1994) coefficients based on style and displacement type.
        
        :param style: Faulting style ("all" or "normal")
        :param norm_disp_type: Normalization type ("AD" or "MD")
        :returns: Dictionary with 'intercept', 'slope', 'sigma' keys
        """
        if style == "all":
            return self._WC94_ALL[norm_disp_type]
        elif style == "normal":  # normal
            return self._WC94_NORMAL[norm_disp_type]
    
    def get_prob(self, d, X_L_ratio, mag, style, norm_disp_type):
        """
        Model of Youngs et al. (2003) for the probability of exceeding
        threshold values of primary displacement [m].

        This method uses numerical integration (epsilon-based weighted sum)
        to convolve the statistical distributions and capture total aleatory
        variability, matching the fdhpy implementation.

        :param d: Target displacement in meters (scalar or array-like, shape (n_displacements,))
        :param X_L_ratio: Ratio of distance from the closest rupture end to the total rupture length
                         (scalar or array-like, shape (n_sites,))
        :param mag: Earthquake magnitude (scalar)
        :param style: Faulting style for WC94 coefficients. Valid options are:
                     - "all": Use WC94 "All styles" coefficients
                     - "normal": Use WC94 "Normal faulting" coefficients
        :param norm_disp_type: Normalization displacement type. Valid options are "AD" or "MD".
        :returns: Probability of exceeding target displacement (m), shape (n_displacements, n_sites).
        """
        # Validate inputs
        style_lower = style.lower() if isinstance(style, str) else str(style).lower()
        if style_lower not in self._ACCEPTED_STYLES:
            raise ValueError(
                f"Invalid style '{style}'. Accepted values are: {', '.join(self._ACCEPTED_STYLES)}"
            )
        
        if norm_disp_type not in self._ACCEPTED_DISP_TYPES:
            raise ValueError(
                f"Invalid displacement type '{norm_disp_type}'. Accepted values are: {', '.join(self._ACCEPTED_DISP_TYPES)}"
            )
        if not np.isscalar(mag):
            raise ValueError("mag must be a scalar value")
        
        # Convert inputs to numpy arrays
        d = np.atleast_1d(d)  # Shape (n_displacements,)
        X_L_ratio = np.atleast_1d(X_L_ratio)  # Shape (n_sites,)
        
        # Fold X/L to [0, 0.5] using modulus reflection, robust to 1±eps
        r = X_L_ratio - np.floor(X_L_ratio)
        X_L_ratio = 0.5 - np.abs(r - 0.5)
        
        # Get WC94 coefficients based on style and displacement type
        coeffs = self._get_wc94_coeffs(style_lower, norm_disp_type)
        mu = coeffs["intercept"] + coeffs["slope"] * mag  # Mean in log10 space
        sigma = coeffs["sigma"]  # Std dev in log10 space
        
        # Get gamma/beta distribution parameters based on x/L
        if norm_disp_type == "AD":
            alpha = np.exp(-0.193 + 1.628 * X_L_ratio)  # Shape (n_sites,)
            beta_param = np.exp(0.009 - 0.476 * X_L_ratio)   # Shape (n_sites,)
        else:  # "MD"
            alpha = np.exp(-0.705 + 1.138 * X_L_ratio)  # Shape (n_sites,)
            beta_param = np.exp(0.421 - 0.257 * X_L_ratio)   # Shape (n_sites,)
        
        # Create epsilon array for integration (matching fdhpy: ±6σ, step 0.1)
        epsilons = np.arange(-self._N_EPS, self._N_EPS + self._DZ, self._DZ)  # Shape (n_eps,)
        prob_eps = norm.pdf(epsilons)  # Shape (n_eps,)
        
        # Compute array of XD (AD or MD) values
        z = np.power(10, mu + epsilons * sigma)  # Shape (n_eps,)
        
        # Compute normalized displacement D/XD for each displacement and XD value
        # d: (n_displacements,), z: (n_eps,)
        # y: (n_displacements, n_eps)
        y = d[:, np.newaxis] / z[np.newaxis, :]
        
        # Initialize output array: (n_displacements, n_sites)
        n_displacements = len(d)
        n_sites = len(X_L_ratio)
        prob_exceeding_d = np.zeros((n_displacements, n_sites))
        
        # For each site, compute the exceedance probability
        for i_site in range(n_sites):
            a_site = alpha[i_site]
            b_site = beta_param[i_site]
            
            if norm_disp_type == "AD":
                # Gamma distribution survival function (1 - CDF = exceedance probability)
                # y: (n_displacements, n_eps)
                sf_matrix = gamma.sf(y, a_site, loc=0, scale=b_site)  # (n_displacements, n_eps)
            else:  # "MD"
                # Beta distribution survival function
                # Truncation to correct for D/MD > 1
                cdf_matrix = beta.cdf(y, a_site, b_site)  # (n_displacements, n_eps)
                cdf_value_at_1 = beta.cdf(1.0, a_site, b_site)
                cdf_matrix = np.where(y > 1, 1.0, cdf_matrix / cdf_value_at_1)
                sf_matrix = 1.0 - cdf_matrix  # (n_displacements, n_eps)
            
            # Compute weighted sum (numerical integration)
            prob_exceeding_d[:, i_site] = np.dot(sf_matrix, prob_eps) * self._DZ
        
        return prob_exceeding_d

    def get_prob_D_AD(self, D_AD, x_L_ratio):
        """
        Model of Youngs et al. (2003) for the probability of normalized displacement
        using the D/AD (average displacement) formulation.

        :param D_AD: Normalized displacement (D/Avg_D), shape (n_displacements, 1) or scalar
        :param x_L_ratio: Ratio of distance from the closest rupture end to the total rupture length,
                         shape (1, n_sites) or scalar
        :returns: Probability of exceeding 'D_AD', shape (n_displacements, n_sites) or scalar
        """
        # Allow tiny numerical noise around bounds
        tol = 1e-10
        if np.any((x_L_ratio < -tol) | (x_L_ratio > 0.5 + tol)):
            raise ValueError("x_L_ratio must be between 0 and 0.5")
        
        # Gamma distribution parameters from Youngs et al. (2003) Eq. 9-10
        a = np.exp(-0.193 + 1.628 * x_L_ratio)  # Shape (1, n_sites) or scalar
        b = np.exp(0.009 - 0.476 * x_L_ratio)   # Shape (1, n_sites) or scalar
        return gamma.sf(D_AD, a, loc=0, scale=b)

    def get_prob_D_MD(self, D_MD, x_L_ratio):
        """
        Model of Youngs et al. (2003) for the probability of normalized displacement
        using the D/MD (maximum displacement) formulation.

        :param D_MD: Normalized displacement (D/Max_D), shape (n_displacements, 1) or scalar
        :param x_L_ratio: Ratio of distance from the closest rupture end to the total rupture length,
                         shape (1, n_sites) or scalar
        :returns: Probability of exceeding 'D_MD', shape (n_displacements, n_sites) or scalar
        """
        tol = 1e-12
        if np.any((x_L_ratio < -tol) | (x_L_ratio > 0.5 + tol)):
            raise ValueError("x_L_ratio must be between 0 and 0.5")
        
        # Beta distribution parameters from Youngs et al. (2003) Eq. 11-12
        a = np.exp(-0.705 + 1.138 * x_L_ratio)  # Shape (1, n_sites) or scalar
        b = np.exp(0.421 - 0.257 * x_L_ratio)   # Shape (1, n_sites) or scalar
        return beta.sf(D_MD, a, b)
