# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#
# Copyright (C) 2012-2024 GEM Foundation
#
# OpenQuake is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""
Moss and Ross (2011) primary surface fault displacement model.
"""

import numpy as np
from scipy.stats import beta
from scipy.stats import gamma
from scipy.stats import norm
from scipy.stats import weibull_min

from openquake.fdha.primary_surf_displ.base import BasePrimarySurfDispl


class MossRoss2011PrimaryFD(BasePrimarySurfDispl):
    """
    Model of Moss and Ross (2011) for exceedance probabilities of primary
    surface fault displacement thresholds.

    The implementation restores the historical model that was removed during
    the May 2025 vectorization cleanup, while returning the current standard
    shape: ``(n_displacements, n_sites)``.
    """

    _ACCEPTED_DISP_TYPES = frozenset(["AD", "MD"])
    _INTEGRATION_DISPLACEMENTS = np.logspace(np.log10(0.001), np.log10(10), 100)

    def get_prob(self, d, X_L_ratio, mag, norm_disp_type):
        """
        Return probability of exceeding primary displacement threshold(s).

        :param d:
            Target displacement in meters, scalar or ``(n_displacements,)``.
        :param X_L_ratio:
            Along-strike position ratio, scalar or ``(n_sites,)``.
        :param mag:
            Earthquake magnitude, scalar.
        :param norm_disp_type:
            Normalization displacement type: ``"AD"`` or ``"MD"``.
        :returns:
            Exceedance probabilities with shape ``(n_displacements, n_sites)``.
        """
        if norm_disp_type not in self._ACCEPTED_DISP_TYPES:
            raise ValueError(
                f"Invalid displacement type '{norm_disp_type}'. Accepted values are: "
                f"{', '.join(sorted(self._ACCEPTED_DISP_TYPES))}"
            )
        if not np.isscalar(mag):
            raise ValueError("mag must be a scalar value")

        d_arr = np.atleast_1d(d).astype(float)
        x_l = np.atleast_1d(X_L_ratio).astype(float)

        if np.any(d_arr <= 0.0):
            raise ValueError("d must be positive")

        r = x_l - np.floor(x_l)
        x_fold = 0.5 - np.abs(r - 0.5)

        integration_displacements = self._INTEGRATION_DISPLACEMENTS
        norm_ratio = d_arr[:, np.newaxis] / integration_displacements[np.newaxis, :]

        if norm_disp_type == "AD":
            mag_weights = self.get_prob_avg_displacement(integration_displacements, mag)
            out = np.zeros((d_arr.size, x_fold.size), dtype=float)
            for idx, x_l_site in enumerate(x_fold):
                out[:, idx] = np.dot(
                    self.get_prob_D_AD(norm_ratio, x_l_site),
                    mag_weights,
                )
            return out

        mag_weights = self.get_prob_max_displacement(integration_displacements, mag)
        out = np.zeros((d_arr.size, x_fold.size), dtype=float)
        for idx, x_l_site in enumerate(x_fold):
            out[:, idx] = np.dot(
                self.get_prob_D_MD(norm_ratio, x_l_site),
                mag_weights,
            )
        return out

    def get_prob_D_AD(self, D_AD, X_L_ratio, variant="gamma"):
        """
        Probability of exceeding normalized displacement D/AD.

        Moss and Ross (2011) present two source-defined distributions on
        D/AD (their Eqs. 6 and 7); both pass goodness-of-fit equally well.

        :param variant: ``"gamma"`` (Eq. 7, default; matches the historical
            implementation) or ``"weibull"`` (Eq. 6).
        """
        self._check_folded_x_l(X_L_ratio)
        v = str(variant).strip().lower()
        if v == "gamma":
            # Moss and Ross (2011) Eq. 7.
            a = np.exp(-30.4 * X_L_ratio**3 + 19.9 * X_L_ratio**2
                       - 2.29 * X_L_ratio + 0.574)
            b = np.exp(50.3 * X_L_ratio**3 - 34.6 * X_L_ratio**2
                       + 6.6 * X_L_ratio - 1.05)
            return gamma.sf(D_AD, a, loc=0, scale=b)
        if v == "weibull":
            # Moss and Ross (2011) Eq. 6.
            k = np.exp(-31.8 * X_L_ratio**3 + 21.5 * X_L_ratio**2
                       - 3.32 * X_L_ratio + 0.431)
            lam = np.exp(17.2 * X_L_ratio**3 - 12.8 * X_L_ratio**2
                         + 3.99 * X_L_ratio - 0.38)
            return weibull_min.sf(D_AD, c=k, scale=lam)
        raise ValueError(
            f"Unknown variant '{variant}'. Accepted: 'gamma', 'weibull'.")

    def get_prob_avg_displacement(self, target_ad, mag):
        """
        Normalized probability mass for average displacement in meters.
        """
        avg_displacement = -2.2192 + 0.3244 * mag
        sigma = 0.17
        return self._normalized_log10_weights(target_ad, avg_displacement, sigma)

    def get_prob_D_MD(self, D_MD, X_L_ratio):
        """
        Probability of exceeding normalized displacement D/MD.
        """
        self._check_folded_x_l(X_L_ratio)

        a = np.exp(0.713 + 0.901 * X_L_ratio)
        b = np.exp(1.74 - 1.86 * X_L_ratio)
        return beta.sf(D_MD, a, b)

    def get_prob_max_displacement(self, target_md, mag):
        """
        Normalized probability mass for maximum displacement in meters.
        """
        max_displacement = -3.1971 + 0.5102 * mag
        sigma = 0.31
        return self._normalized_log10_weights(target_md, max_displacement, sigma)

    @staticmethod
    def _check_folded_x_l(X_L_ratio):
        tol = 1e-12
        if np.any((X_L_ratio < -tol) | (X_L_ratio > 0.5 + tol)):
            raise ValueError("X_L_ratio must be folded between 0 and 0.5")

    def _normalized_log10_weights(self, target_displacement, mean, sigma):
        target_displacement = np.atleast_1d(target_displacement).astype(float)
        if np.any(target_displacement <= 0.0):
            raise ValueError("target displacement values must be positive")

        prob = norm.pdf(np.log10(target_displacement), loc=mean, scale=sigma)
        grid_prob = norm.pdf(
            np.log10(self._INTEGRATION_DISPLACEMENTS),
            loc=mean,
            scale=sigma,
        )
        return prob / np.sum(grid_prob)
