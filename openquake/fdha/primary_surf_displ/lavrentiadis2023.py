# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#
# Copyright (C) 2024-2026 Yen-Shin Chen, OGS
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""
Module :mod:`openquake.fdha.primary_surf_displ.lavrentiadis2023` implements
the model of Lavrentiadis and Abrahamson (2023) in two classes, one per
displacement definition (the class choice IS the definition; the version,
full rupture or individual segment, is the ``output_type`` parameter):

- :class:`Lavrentiadis2023PrimaryFD_aggregate` - the AGGREGATE-definition
  variants, ``output_type`` ``disp_agg_prime`` (full rupture, default) or
  ``disp_agg_seg`` (individual segment);
- :class:`Lavrentiadis2023PrimaryFD_principal` - the sum-of-principal
  variants, ``output_type`` ``disp_prnc_prime`` (full rupture, default) or
  ``disp_prnc_seg`` (individual segment).

Equation, table and page numbers refer to Lavrentiadis and Abrahamson
(2023). The regression terms follow the authors' reference implementation
(NHR3-UCLA/LA23_PFDHA_model, commit c3ba1136, MIT licence), which differs
from the printed paper in two places: the cross term of Eq. 34 (the code's
2*rho*phi_agg*phi_b2; the printed phi_b2**2 is dimensionally inconsistent)
and the predictor of Eq. 32 for the full-rupture model (the code uses
mu'_agg of Eq. 22). Items tagged ``[fdhpy-IMP]`` follow fdhpy 1.0.3
(Sarmiento et al. 2025) where the paper gives no equation; items tagged
``[DERIVED]`` follow from the equations cited with them.
"""

import logging

import numpy as np
from scipy import integrate
from scipy import stats as scipystats

from openquake.fdha.params import check_bool, check_style
from openquake.fdha.primary_surf_displ.base import BasePrimarySurfDispl

_LOGGER = logging.getLogger(__name__)

#: Styles of faulting with LA23 coefficients; the coefficient arrays below
#: follow this order.
STYLES = ("normal", "strike-slip", "reverse")

#: Magnitude range of the model: LA23 p. 24 (M5.0-M8.5 for all three
#: styles); fdhpy ``_CONDITIONS``. Outside it the model warns and
#: extrapolates.
MAGNITUDE_RANGE = (5.0, 8.5)

AGGREGATE_OUTPUT_TYPES = ("disp_agg_prime", "disp_agg_seg")
PRINCIPAL_OUTPUT_TYPES = ("disp_prnc_prime", "disp_prnc_seg")

#: Individual-segment variants: functions of X_seg/L_seg (Eq. 14), not of
#: the full-rupture X/L.
SEGMENT_OUTPUT_TYPES = ("disp_agg_seg", "disp_prnc_seg")

# Model coefficients, copied from the authors' implementation; arrays are
# ordered as STYLES (normal, strike-slip, reverse).
# general aggregate slip, Eq. 8 (Table 3)
c_0 = 0.272
c_1 = 0.913
c_2 = -2.128
c_3 = 1.15
# tapering, Eqs. 9-14 (Table 3)
c_5 = np.array([0.7, 0.7, -0.2])
c_6 = np.array([-0.262, -0.314, 0.0])
c_7 = 0.033
M_1 = np.array([6.25, 6.50, 5.00])
# effect of segmentation on the median, Eqs. 19-21 (Table 4)
c_10 = np.array([-0.406, -0.394, -0.849])
c_11 = np.array([-8.205, -2.950, 8.026])
c_12 = np.array([-165.0, -64.3, 237.2])
c_13 = np.array([-1372, -833, 927])
c_14 = np.array([-3013, -2145, 831])
c_14a = np.array([-0.904, -1.010, -3.555])
c_15 = np.array([-0.440, -0.155, -0.087])
c_16 = np.array([0.0767, 0.0263, 0.0148])
c_17 = np.array([0.0275, 0.0177, 0.0087])
c_17a = np.array([0.00276, 0.00534, 0.00234])
# effect of segmentation on the aleatory variability, Eq. 23 (Table 4)
c_18 = np.array([-0.215, -0.201, -0.165])
c_19 = np.array([0.0430, 0.0364, 0.0304])
c_20 = np.array([0.00574, 0.0102, 0.00773])
# gap probability, Eqs. 25-28 (Table 4)
c_21 = np.array([-0.082, -0.135, -0.151])
c_22 = np.array([0.027, 0.027, 0.033])
c_23 = np.array([-0.0088, 0.0063, 0.0032])
c_24a = np.array([0.020, 0.040, 0.04])
c_24b = np.array([0.00, 0.0050, 0.00])
c_24c = np.array([0.0247, 0.00273, 0.00843])
c_25 = np.array([0.60, 0.41, 0.43])
c_26a = np.array([0.162, 0.153, 0.156])
c_26b = np.array([0.084, 0.017, 0.031])
c_26c = np.array([0.0080, 0.00312, 0.00575])
# principal displacement, Eqs. 32-37 (Table 5)
b_0 = np.array([-3.240, 0.867, 1.65])
b_1 = np.array([9.105, 3.767, 1.349])
b_2 = np.array([-0.097, -0.048, -0.062])
b_3 = np.array([0.286, 0.191, 0.28])
b_4 = np.array([1.195, 1.058, 0.0])
b_5 = np.array([-1.526, -1.418, 0.0])
phi_b2 = 0.11
rho_b2 = -0.15
# maximum aggregate displacement, Eq. 38 (Table 6)
e_1 = 0.326
e_2 = 0.41
e_3 = np.array([-0.0166, -0.0144, -0.0172])
e_4 = np.array([0.0244, -0.0109, 0.0094])

_WARNED = set()


def _warn_once(key, message, *args):
    """Log ``message`` once per ``key`` per process (hazard loops call the
    model once per rupture)."""
    if key not in _WARNED:
        _WARNED.add(key)
        _LOGGER.warning(message, *args)


def _components(X_L_ratio, mag, style):
    """
    Every LA23 term at the locations ``X_L_ratio`` for one magnitude and
    style. Location-dependent terms have the shape of
    ``np.atleast_1d(X_L_ratio)``; the others are scalars.
    """
    i = STYLES.index(style)
    xl = np.atleast_1d(np.asarray(X_L_ratio, dtype=float))
    # fold x/L to [0, 0.5] by modulus reflection, robust to 1 +- eps
    r = xl - np.floor(xl)
    xl = 0.5 - np.abs(r - 0.5)

    # Eq. 8, wavenumber-domain aggregate slip
    d_agg = np.exp(
        c_0 + c_1 * (xl - 0.3) + c_2 * (xl - 0.3) ** 2 + c_3 * (mag - 7.0))
    # Eq. 9, X/L taper; Eq. 13, magnitude taper
    xl1 = 0.15 - 0.10 * np.minimum(np.maximum(mag - 7.0, 0.0), 1)
    T_xl = np.minimum(xl - xl1, 0) / xl1
    T_M = np.maximum(np.minimum((7.0 - mag) / (7.0 - M_1[i]), 1.0), 0.0)
    # Eq. 14, median D_agg**0.3 of an individual segment
    mu_agg_seg = (d_agg * np.exp(c_5[i] * T_xl)) ** 0.3 + c_6[i] * T_M + c_7

    # Eqs. 19-22, median D_agg**0.3 of the full rupture
    xl_offset1 = np.minimum(xl - 0.3, 0)
    xl_offset2 = np.maximum(xl - 0.4, 0)
    f_NDmu = (
        c_11[i] * xl_offset1
        + c_12[i] * xl_offset1**2
        + c_13[i] * xl_offset1**3
        + c_14[i] * xl_offset1**4
        + c_14a[i] * xl_offset2
    )
    f_NDmu = c_10[i] + f_NDmu
    Dmu_max = (c_15[i] + c_16[i] * mag + c_17[i] * (mag - 6.7) ** 2
               + c_17a[i] * (mag - 6.7) ** 3)
    mu_agg_prime = mu_agg_seg + Dmu_max * f_NDmu

    # Eqs. 25-28, probability of being in a gap (full rupture only)
    c_24 = c_24a[i] + c_24b[i] * (mag - 5) + c_24c[i] * (mag - 5) ** 3
    c_26 = c_26a[i] + c_26b[i] * (mag - 5) + c_26c[i] * (mag - 5) ** 3
    P_gap_max = c_21[i] + c_22[i] * mag + c_23[i] * (mag - 6.5) ** 2
    f_NPGap = 20.0 * c_24 * np.minimum(np.maximum(xl - 0.10, 0.0), 0.05)
    f_NPGap += (0.13**-1 * (1.0 - c_24)
                * np.minimum(np.maximum(xl - 0.15, 0.0), 0.13))
    f_NPGap -= (10.0 * (1.0 - c_25[i])
                * np.minimum(np.maximum(xl - 0.30, 0.0), 0.10))
    f_NPGap -= (10.0 * (c_25[i] - c_26)
                * np.minimum(np.maximum(xl - 0.40, 0.0), 0.10))
    p_gap = P_gap_max * f_NPGap

    # Eq. 32, P(D_P = 0 | not in a gap). Full rupture: predictor mu'_agg
    # (Eq. 22) as in the authors' code. Individual segment: predictor mu_agg
    # (Eq. 14), as printed in Eq. 32 [DERIVED, decision B2]
    p_zero_prime = 1 / (1 + np.exp(b_0[i] + b_1[i] * mu_agg_prime))
    p_zero_seg = 1 / (1 + np.exp(b_0[i] + b_1[i] * mu_agg_seg))

    # Eq. 33, median D_P**0.3 of the full rupture; the individual-segment
    # median applies the same b_2 offset to Eq. 14 [fdhpy-IMP; p. 23]
    mu_prnc_prime = np.maximum(mu_agg_prime + b_2[i], 0.0)
    mu_prnc_seg = np.maximum(mu_agg_seg + b_2[i], 0.0)

    # Eqs. 15-16, within- and between-event variability; Eq. 34 (code form,
    # see the module docstring) with tau_P = tau_agg (p. 24); Eq. 23
    phi_agg = np.maximum(np.minimum(0.120 + 0.150 * (mag - 6.0), 0.270), 0.120)
    tau_agg = np.maximum(np.minimum(0.115 + 0.060 * (mag - 6.0), 0.205), 0.115)
    phi_prnc = np.sqrt(phi_agg**2 + phi_b2**2 + 2 * rho_b2 * phi_agg * phi_b2)
    phi_add = c_18[i] + c_19[i] * mag + c_20[i] * (mag - 6.7) ** 2

    return {
        "xl": xl,
        "mu_agg_seg": mu_agg_seg,
        "mu_agg_prime": mu_agg_prime,
        "mu_prnc_seg": mu_prnc_seg,
        "mu_prnc_prime": mu_prnc_prime,
        "p_gap": p_gap,
        "p_zero_prime": p_zero_prime,
        "p_zero_seg": p_zero_seg,
        "phi_agg": phi_agg,
        "tau_agg": tau_agg,
        "phi_prnc": phi_prnc,
        "phi_add": phi_add,
        # p. 16 (below Eq. 16), individual segment
        "sig_agg_seg": np.sqrt(tau_agg**2 + phi_agg**2),
        # Eq. 24, full rupture
        "sig_agg_prime": np.sqrt(tau_agg**2 + phi_agg**2 + phi_add**2),
        # Eq. 34 with Eq. 24 (Delta-adjustments apply to D_P, p. 24)
        "sig_prnc_prime": np.sqrt(tau_agg**2 + phi_prnc**2 + phi_add**2),
        # [fdhpy-IMP] individual segment: no phi_add
        "sig_prnc_seg": np.sqrt(tau_agg**2 + phi_prnc**2),
    }


def _censored_mean(mu, sigma):
    """
    E[max(Y, 0)**(10/3)] for Y ~ N(mu, sigma**2): the mean of the
    power-normal displacement, left-truncated at zero as in Eq. 17 and
    integrated without renormalisation as in Eq. 35 [fdhpy-IMP: quadrature].
    """
    mu, sigma = np.broadcast_arrays(np.asarray(mu, float),
                                    np.asarray(sigma, float))
    out = np.empty(mu.shape)
    for idx in np.ndindex(mu.shape):
        m, s = float(mu[idx]), float(sigma[idx])
        lo, hi = max(0.0, m - 12.0 * s), max(0.0, m + 12.0 * s)
        if hi <= 0.0:
            out[idx] = 0.0
            continue
        out[idx] = integrate.quad(
            lambda y: y ** (10 / 3) * scipystats.norm.pdf(y, m, s), lo, hi,
            points=[m] if lo < m < hi else None, epsabs=0.0, epsrel=1e-12,
            limit=200)[0]
    return out


def _displacement(percentile, mu, sigma, w_nonzero):
    """
    Displacement (m) at ``percentile`` (or the mean for ``-1``) of the
    mixture of a point mass at zero, weight 1 - w_nonzero, and the
    power-normal distribution of Eq. 17 left-truncated at zero.

    The percentile is the root of P(D > d) = 1 - percentile, with P(D > d)
    the zero-scaled exceedance of Eqs. 25 and 31 [DERIVED from Eq. 31,
    decision G3]; the mean is w_nonzero times the censored mean of the
    non-zero part (Eq. 35).
    """
    mu, sigma, w = np.broadcast_arrays(
        np.asarray(mu, float), np.asarray(sigma, float),
        np.asarray(w_nonzero, float))
    if percentile == -1:
        return w * _censored_mean(mu, sigma)
    if not 0.0 < percentile < 1.0:
        raise ValueError(
            f"percentile must be in (0, 1), or -1 for the mean; got {percentile}")
    # P(D <= d) = 1 - w * (1 - Phi(z)) for d >= 0, so Phi(z) = 1 - (1 - p) / w
    with np.errstate(divide="ignore", invalid="ignore"):
        q = np.where(w == 1.0, percentile, 1.0 - (1.0 - percentile) / w)
    q = np.nan_to_num(q, nan=0.0, neginf=0.0)
    y = mu + sigma * scipystats.norm.ppf(np.clip(q, 0.0, 1.0))
    y = np.where(q > 0.0, y, 0.0)
    return np.maximum(y, 0.0) ** (1 / 0.3)


def _squeeze(prob):
    """Normalize the (n_displacements, n_locations) output shape: a single
    site gives (n_displacements,), a single displacement (n_locations,)."""
    if getattr(prob, "ndim", 0) == 2:
        if prob.shape[1] == 1:
            return prob[:, 0]
        if prob.shape[0] == 1:
            return prob[0, :]
    return prob


class Lavrentiadis2023PrimaryFD_aggregate(BasePrimarySurfDispl):
    """Aggregate fault-displacement model of Lavrentiadis and Abrahamson (2023).

    Model of aggregate fault displacement as a function of magnitude,
    normalized along-strike position, and faulting style, run in the
    principal (primary_surf_displ) slot.

    References
    ----------
    Lavrentiadis, G., and Abrahamson, N.A. (2023). Fault-displacement models
    for aggregate and principal displacements. Earthquake Spectra, 41(4),
    2806-2837. https://doi.org/10.1177/87552930231201531

    Model contract: DISPLACEMENT_DEFINITION = "aggregate",
    DISPLACEMENT_COMPONENT = "net" -- STATIC, the class choice IS the
    definition. Lavrentiadis & Abrahamson (2023) develop their model on the
    FDHI aggregate displacement (net slip summed across principal and
    distributed ruptures in the measurement aperture); Sarmiento et al.
    (2025, Earthquake Spectra) Table 1 lists LA23 under the aggregate
    definition. This class serves ONLY the aggregate-definition variants:
    ``output_type = "disp_agg_prime"`` (default, full rupture, the
    "simplified FDM without segmentation") or ``"disp_agg_seg"`` (individual
    segment). The sum-of-principal variants (``disp_prnc_prime``,
    ``disp_prnc_seg``; principal-strand slip WITHOUT distributed ruptures)
    are a different physical quantity and live in
    :class:`Lavrentiadis2023PrimaryFD_principal`; requesting them here raises
    ValueError (and is rejected at logic-tree validation time, FDLT-015).
    As an aggregate model, this class runs single-bucket and forbids
    secondary-slot models (FDLT-013).

    Zero displacement (``include_zero_slip``, default False; fdhpy defaults
    to True): ``disp_agg_prime`` scales the exceedance by 1 - P(Gap)
    (Eq. 25). ``disp_agg_seg`` has no zero term -- a single segment has no
    gap (p. 21) and P(D_P = 0) concerns the principal displacement only --
    so the flag is ignored with a warning.

    Magnitudes outside M 5.0-8.5 log a warning and are extrapolated.
    """

    DISPLACEMENT_DEFINITION = "aggregate"
    DISPLACEMENT_COMPONENT = "net"

    #: ``output_type`` values this class accepts; the first is the default.
    OUTPUT_TYPES = AGGREGATE_OUTPUT_TYPES
    SEGMENT_OUTPUT_TYPES = SEGMENT_OUTPUT_TYPES

    def __init__(self, style=None, output_type=None, include_zero_slip=None):
        """
        :param style: optional faulting style pinned by the logic-tree
            branch ('normal', 'strike-slip' or 'reverse'); ``None`` defers
            to the ``get_prob`` call (legacy default: 'normal').
        :param output_type: optional output-type selector pinned by the
            logic-tree branch (one of :attr:`OUTPUT_TYPES`); ``None``
            defers to the call (default: the first of :attr:`OUTPUT_TYPES`).
        :param include_zero_slip: optional flag pinned by the logic-tree
            branch; ``None`` defers to the call (default: False).
        """
        self.style = check_style(type(self).__name__, style, accepted=STYLES)
        if output_type is not None:
            output_type = str(output_type)
            self.check_output_type(output_type)
        self.output_type = output_type
        self.include_zero_slip = check_bool(
            type(self).__name__, "include_zero_slip", include_zero_slip)

    @classmethod
    def check_output_type(cls, output_type):
        """Raise ValueError unless ``output_type`` is served by this class
        (logic-tree pin, job parameter or call argument)."""
        if output_type in cls.OUTPUT_TYPES:
            return
        if output_type in PRINCIPAL_OUTPUT_TYPES:
            raise ValueError(
                f"{cls.__name__} is the AGGREGATE-definition model "
                f"({' / '.join(cls.OUTPUT_TYPES)}); the sum-of-principal "
                f"metric {output_type} is a different displacement "
                "definition (Sarmiento et al. 2025, Table 1). Select the "
                "Lavrentiadis2023PrimaryFD_principal model class instead.")
        raise ValueError(
            f"Invalid output_type '{output_type}' for {cls.__name__}; "
            f"accepted: {', '.join(cls.OUTPUT_TYPES)}")

    def _resolve(self, mag, style, output_type, include_zero_slip):
        """Fall back to the constructor pins, then the defaults; validate."""
        name = type(self).__name__
        if style is None:
            style = self.style
        if style is None:
            _warn_once((name, "style"),
                       "%s: no style given; using the legacy default "
                       "'normal'", name)
            style = "normal"
        style = check_style(name, style, accepted=STYLES)
        if output_type is None:
            output_type = self.output_type
        if output_type is None:
            output_type = self.OUTPUT_TYPES[0]
        output_type = str(output_type)
        self.check_output_type(output_type)
        if include_zero_slip is None:
            include_zero_slip = self.include_zero_slip
        include_zero_slip = bool(check_bool(
            name, "include_zero_slip", include_zero_slip))
        mag_lo, mag_hi = np.min(mag), np.max(mag)
        if mag_lo < MAGNITUDE_RANGE[0] or mag_hi > MAGNITUDE_RANGE[1]:
            _warn_once((name, "magnitude", round(float(mag_lo), 3),
                        round(float(mag_hi), 3)),
                       "%s: magnitude %s is outside the model range "
                       "M %.1f-%.1f (LA23 p. 24); extrapolating", name,
                       f"{mag_lo:g}" if mag_lo == mag_hi
                       else f"{mag_lo:g}-{mag_hi:g}", *MAGNITUDE_RANGE)
        return style, output_type, include_zero_slip

    def _distribution(self, X_L_ratio, mag, style, output_type,
                      include_zero_slip):
        """
        (mu, sigma, w_nonzero) of D**0.3 at each location: the median and
        standard deviation of the power-normal part (Eq. 17) and the
        probability of a non-zero displacement.
        """
        comp = _components(X_L_ratio, mag, style)
        if output_type == "disp_agg_prime":
            # Eq. 22, Eq. 24; zero term 1 - P(Gap), Eq. 25
            mu, sigma = comp["mu_agg_prime"], comp["sig_agg_prime"]
            w = 1 - comp["p_gap"]
        elif output_type == "disp_agg_seg":
            # Eq. 14, p. 16; no gap for a single segment (p. 21)
            mu, sigma = comp["mu_agg_seg"], comp["sig_agg_seg"]
            w = 1.0
            if include_zero_slip:
                _warn_once(
                    (type(self).__name__, "zero", output_type),
                    "%s: include_zero_slip is ignored for disp_agg_seg: a "
                    "single segment has no gap (LA23 p. 21) and "
                    "P(D_P = 0) (Eq. 32) concerns the principal "
                    "displacement only", type(self).__name__)
        elif output_type == "disp_prnc_prime":
            # Eqs. 31-34
            mu, sigma = comp["mu_prnc_prime"], comp["sig_prnc_prime"]
            w = (1 - comp["p_zero_prime"]) * (1 - comp["p_gap"])
        elif output_type == "disp_prnc_seg":
            # [fdhpy-IMP] median and sigma; zero term 1 - P(D_P = 0) with
            # P(Gap) = 0 for a single segment (p. 21) [DERIVED from
            # Eqs. 31-32, decision B2]
            mu, sigma = comp["mu_prnc_seg"], comp["sig_prnc_seg"]
            w = 1 - comp["p_zero_seg"]
        else:
            raise ValueError(f"Invalid output_type '{output_type}'")
        if not include_zero_slip:
            w = 1.0
        mu, sigma, w = np.broadcast_arrays(mu, sigma, w)
        return mu, sigma, w

    def get_prob(self, d, X_L_ratio, mag, style=None, output_type=None,
                 include_zero_slip=None):
        """
        Probability of exceeding displacement thresholds.

        :param d:
            Target displacement in meters.
        :param X_L_ratio:
            Ratio of distance from the closest rupture end to the total rupture
            length (any value; folded to [0, 0.5]).
        :param mag:
            Earthquake magnitude.
        :param style:
            Style of faulting ("normal", "strike-slip" or "reverse").
        :param output_type:
            Displacement metric to evaluate, one of :attr:`OUTPUT_TYPES`.
        :param include_zero_slip:
            If ``True`` the probability accounts for the zero-displacement
            terms (see the class docstring).
        :returns:
            Probability of exceeding ``d``, shape (n_displacements,
            n_locations) squeezed for a single site or displacement.
        :raises ValueError:
            If ``style`` or ``output_type`` is invalid for this class.
        """
        style, output_type, include_zero_slip = self._resolve(
            mag, style, output_type, include_zero_slip)
        return self._evaluate(
            d, X_L_ratio, mag, style=style, output_type=output_type,
            include_zero_slip=include_zero_slip)

    def _evaluate(self, d, X_L_ratio, mag, style, output_type,
                  include_zero_slip):
        """Shared implementation for all published output_type variants
        (the public classes pin/restrict output_type; this does not)."""
        mu, sigma, w = self._distribution(
            X_L_ratio, mag, style, output_type, include_zero_slip)
        z = np.atleast_1d(d)[:, np.newaxis] ** 0.3
        return _squeeze(w * scipystats.norm.sf(x=z, loc=mu, scale=sigma))

    def get_cdf(self, d, X_L_ratio, mag, style=None, output_type=None,
                include_zero_slip=None):
        """P(D <= d) = 1 - :meth:`get_prob`; same parameters and shape."""
        return 1.0 - self.get_prob(
            d, X_L_ratio, mag, style=style, output_type=output_type,
            include_zero_slip=include_zero_slip)

    def get_displ_site(self, X_L_ratio, mag, percentile=0.5, style=None,
                       output_type=None, include_zero_slip=None):
        """
        Deterministic displacement (m) at each location.

        :param percentile: aleatory quantile in (0, 1), or -1 for the mean.
            With the zero terms the displacement is the percentile of the
            mixture of zero displacement and the power-normal distribution
            (the root of P(D > d) = 1 - percentile), not fdhpy's percentile
            scaled by the non-zero probability.
        :returns: array of shape ``np.atleast_1d(X_L_ratio).shape``
        """
        style, output_type, include_zero_slip = self._resolve(
            mag, style, output_type, include_zero_slip)
        mu, sigma, w = self._distribution(
            X_L_ratio, mag, style, output_type, include_zero_slip)
        return _displacement(percentile, mu, sigma, w)

    def get_displ_profile(self, mag, percentile=0.5, xl_step=0.05,
                          style=None, output_type=None,
                          include_zero_slip=None):
        """
        Displacement profile along the full rupture.

        :param xl_step: x/L increment. The grid always closes at x/L = 1;
            fdhpy returns no profile when the step does not divide 1
            [fdhpy-IMP deviation].
        :returns: ``(xl, displ)`` arrays
        """
        if not 0.0 < xl_step <= 1.0:
            raise ValueError(f"xl_step must be in (0, 1]; got {xl_step}")
        n = int(np.floor(1.0 / xl_step + 1e-9))
        xl = np.arange(n + 1) * xl_step
        if xl[-1] < 1.0 - 1e-9:
            xl = np.append(xl, 1.0)
        return xl, self.get_displ_site(
            xl, mag, percentile=percentile, style=style,
            output_type=output_type, include_zero_slip=include_zero_slip)

    def get_displ_avg(self, mag, style=None):
        """Average displacement: defined for the principal displacement only
        (Eqs. 35-37); see :class:`Lavrentiadis2023PrimaryFD_principal`."""
        raise ValueError(
            f"{type(self).__name__}: LA23 gives the average displacement for "
            "the principal displacement only (Eqs. 35-37); use "
            "Lavrentiadis2023PrimaryFD_principal.get_displ_avg")

    def get_displ_max(self, mag, percentile=0.5, style=None):
        """
        Maximum aggregate displacement (m) of the full rupture, Eqs. 38-39
        and Table 6.

        :param percentile: aleatory quantile in (0, 1), or -1 for the mean.
        """
        style, output_type, _ = self._resolve(mag, style, None, None)
        if output_type != "disp_agg_prime":
            raise ValueError(
                f"{type(self).__name__}: LA23 gives the maximum displacement "
                "for the full rupture only (Eq. 38); this instance is pinned "
                f"to output_type = {output_type}")
        if self.include_zero_slip:
            _warn_once((type(self).__name__, "zero", "max"),
                       "%s: include_zero_slip is ignored for the maximum "
                       "displacement", type(self).__name__)
        mu, sigma = self._max_disp_params(mag, style)
        return float(_displacement(percentile, mu, sigma, 1.0))

    @staticmethod
    def _max_disp_params(mag, style):
        """Median and sigma of MD_agg**0.3, Eqs. 38-39."""
        i = STYLES.index(style)
        mu_agg_prime = _components(0.25, mag, style)["mu_agg_prime"][0]
        dm_pwr = e_1
        dm_pwr += e_2 * np.minimum(np.maximum(mag - 6.0, 0.0), 1.0)
        dm_pwr += (e_3[i] * np.maximum(mag - 7.0, 0.0)
                   + e_4[i] * np.maximum(mag - 7.0, 0.0) ** 2)
        sig_dm = 0.13 + 0.095 * np.minimum(np.maximum(mag - 6.0, 0.0), 1.0)
        sig_dm += 0.050 * np.minimum(np.maximum(mag - 7.0, 0.0), 0.5)
        return mu_agg_prime + dm_pwr, sig_dm

    def get_stat_params(self, X_L_ratio, mag, style=None, output_type=None):
        """
        Parameters of the power-normal distribution of D**0.3 (Eq. 17):
        ``{"mu": ..., "sigma": ..., "distribution": ...}``; the zero terms
        are given by :meth:`get_p_gap` and :meth:`get_p_zero_slip`.
        """
        style, output_type, _ = self._resolve(mag, style, output_type, False)
        mu, sigma, _ = self._distribution(
            X_L_ratio, mag, style, output_type, False)
        return {
            "mu": mu,
            "sigma": sigma,
            "distribution": "D**0.3 ~ Normal(mu, sigma), left-truncated at "
                            "zero (Eq. 17)",
        }

    def get_p_gap(self, X_L_ratio, mag, style=None):
        """P(Gap), the probability of a site in a gap between segments of
        the full rupture, Eq. 25."""
        style, _, _ = self._resolve(mag, style, None, False)
        return _components(X_L_ratio, mag, style)["p_gap"]

    def get_p_zero_slip(self, X_L_ratio, mag, style=None, segment=False):
        """
        P(D_P = 0 | not in a gap), Eq. 32: with the full-rupture predictor
        mu'_agg (authors' code, fdhpy) or, for ``segment=True``, the
        individual-segment predictor mu_agg of Eq. 14.
        """
        style, _, _ = self._resolve(mag, style, None, False)
        comp = _components(X_L_ratio, mag, style)
        return comp["p_zero_seg"] if segment else comp["p_zero_prime"]

    def get_sigma_mu_agg(self, mag, style=None):
        """Epistemic standard deviation of the median D_agg**0.3, Eq. 29."""
        style, _, _ = self._resolve(mag, style, None, False)
        if mag >= 7.1:
            return 0.035 + 0.025 * (mag - 7.1)
        slope = 0.064 if style == "normal" else 0.036
        return 0.035 + slope * (7.1 - mag)

    # Legacy interface of the authors' implementation (x in metres along a
    # rupture of length srl). Kept as thin wrappers.

    def LavrentiadisAbrahamson2023SlipProfile(self, x_array, mag, srl, sof="Strike-Slip"):
        """
        Lavrentiadis and Abrahamson 2023 fault displacement model for
        aggregate and principal displacement.

        Parameters
        ----------
        x_array : np.array()
            Along strike location - unnormalized (m).
        mag : real
            Moment Magnitude.
        srl : real
            Surface Rupture Length (m).
        sof : string, optional
            Style of faulting. The valid options are: Normal, Strike-Slip, and Reverse
            The default is 'Strike-Slip'.

        Returns
        -------
        disp_agg_prime : np.array()
            Aggregate displacements for entire event rupture (m).
        disp_prnc_prime : np.array()
            Principal displacements for the entire event rupture (m).
        disp_agg_seg : np.array()
            Aggregate displacements for single segment (m).
        sig_agg : np.array()
            Total aggregate aleatory variability (tau_agg, phi_agg, phi_add).
        sig_prnc : np.array()
            Total principal aleatory variability (tau_agg, phi_prnc, phi_add).
        phi_agg : np.array()
            Within-event aggregate aleatory variability.
        phi_prnc : np.array()
            Within-event principal aleatory variability.
        tau_agg : np.array()
            Between-event aleatory variability.
        phi_add : np.array()
            Additional aleatory variability due to segmentation.
        P_gap : np.array()
            Segment gap probability.
        P_zero_slip : np.array()
            Zero slip probability.
        """
        style = check_style(type(self).__name__, sof, accepted=STYLES)
        comp = _components(np.asarray(x_array, dtype=float) / srl, mag, style)
        return (
            comp["mu_agg_prime"] ** (1 / 0.3),
            comp["mu_prnc_prime"] ** (1 / 0.3),
            comp["mu_agg_seg"] ** (1 / 0.3),
            comp["sig_agg_prime"],
            comp["sig_prnc_prime"],
            comp["phi_agg"],
            comp["phi_prnc"],
            comp["tau_agg"],
            comp["phi_add"],
            comp["p_gap"],
            comp["p_zero_prime"],
        )

    def LavrentiadisAbrahamson2023SlipProfilePrc(self, x_array, mag, srl, sof="Strike-Slip", prc=0.5):
        """
        Aggregate (full rupture), principal (full rupture) and aggregate
        (individual segment) displacement profiles (m) at percentile
        ``prc``, without the zero terms.
        """
        style = check_style(type(self).__name__, sof, accepted=STYLES)
        xl = np.asarray(x_array, dtype=float) / srl
        return tuple(
            _displacement(prc, *self._distribution(xl, mag, style, ot, False))
            for ot in ("disp_agg_prime", "disp_prnc_prime", "disp_agg_seg"))

    def LavrentiadisAbrahamson2023AvgDisp(self, mag, srl, sof="Strike-Slip"):
        """Median average principal displacement (m) and the ratio AD'_P /
        D'_P(X/L = 0.25), Eqs. 36-37."""
        style = check_style(type(self).__name__, sof, accepted=STYLES)
        return Lavrentiadis2023PrimaryFD_principal._avg_disp(mag, style)

    def LavrentiadisAbrahamson2023MaxDisp(self, mag, srl, sof="Strike-Slip"):
        """Median maximum aggregate displacement (m) and its sigma in
        D**0.3 units, Eqs. 38-39."""
        style = check_style(type(self).__name__, sof, accepted=STYLES)
        mu, sigma = self._max_disp_params(mag, style)
        return float(_displacement(0.5, mu, sigma, 1.0)), float(sigma)


class Lavrentiadis2023PrimaryFD_principal(Lavrentiadis2023PrimaryFD_aggregate):
    """Sum-of-principal variants of Lavrentiadis and Abrahamson (2023).

    Identical regression machinery as
    :class:`Lavrentiadis2023PrimaryFD_aggregate`, for the principal-strand
    slip summed within the measurement aperture, WITHOUT distributed
    ruptures. Following the Petersen2011PrimaryFD_* variant idiom, the class
    choice IS the definition; ``output_type`` selects only the version:
    ``"disp_prnc_prime"`` (default, full rupture, Eqs. 31-34) or
    ``"disp_prnc_seg"`` (individual segment: median max(mu_agg + b_2, 0) with
    mu_agg of Eq. 14 and sigma sqrt(phi_P**2 + tau_agg**2) [fdhpy-IMP]). The
    aggregate values are rejected (also at logic-tree validation time,
    FDLT-015).

    Model contract: DISPLACEMENT_DEFINITION = "sum-of-principal",
    DISPLACEMENT_COMPONENT = "net" -- Lavrentiadis & Abrahamson (2023)
    derive the principal-prime profile from the aggregate one (their b_2
    offset and zero-slip probability); like Chiou et al. (2025)'s D_SP it
    excludes distributed ruptures, so secondary-slot models remain
    legitimate alongside this class (it is NOT aggregate; FDLT-013 does not
    apply).

    Zero displacement (``include_zero_slip``, default False; fdhpy defaults
    to True): ``disp_prnc_prime`` scales the exceedance by
    (1 - P(D_P = 0)) (1 - P(Gap)) (Eq. 31). ``disp_prnc_seg`` scales it by
    1 - P(D_P = 0) with the Eq. 14 predictor of Eq. 32 and P(Gap) = 0 for a
    single segment (p. 21) [DERIVED from Eqs. 31-32]; fdhpy ignores the flag
    for individual segments.
    """

    DISPLACEMENT_DEFINITION = "sum-of-principal"
    DISPLACEMENT_COMPONENT = "net"

    OUTPUT_TYPES = PRINCIPAL_OUTPUT_TYPES

    @classmethod
    def check_output_type(cls, output_type):
        if output_type in cls.OUTPUT_TYPES:
            return
        if output_type in AGGREGATE_OUTPUT_TYPES:
            raise ValueError(
                f"{cls.__name__} evaluates the sum-of-principal metric "
                f"({' / '.join(cls.OUTPUT_TYPES)}); got '{output_type}', an "
                "aggregate-definition metric. Select "
                "Lavrentiadis2023PrimaryFD_aggregate for the aggregate "
                "variants.")
        raise ValueError(
            f"Invalid output_type '{output_type}' for {cls.__name__}; "
            f"accepted: {', '.join(cls.OUTPUT_TYPES)}")

    def get_prob(self, d, X_L_ratio, mag, style=None,
                 include_zero_slip=None, output_type=None):
        """
        Probability of exceeding sum-of-principal displacement thresholds;
        parameters as in
        :meth:`Lavrentiadis2023PrimaryFD_aggregate.get_prob`, with
        ``output_type`` one of :attr:`OUTPUT_TYPES` (default
        ``disp_prnc_prime``).
        """
        return super().get_prob(
            d, X_L_ratio, mag, style=style, output_type=output_type,
            include_zero_slip=include_zero_slip)

    def get_displ_avg(self, mag, style=None):
        """Median average principal displacement (m) of the full rupture,
        Eqs. 36-37 (no aleatory variability is given for it)."""
        style, output_type, _ = self._resolve(mag, style, None, None)
        if output_type != "disp_prnc_prime":
            raise ValueError(
                f"{type(self).__name__}: LA23 gives the average displacement "
                "for the full rupture only (Eqs. 35-37); this instance is "
                f"pinned to output_type = {output_type}")
        if self.include_zero_slip:
            _warn_once((type(self).__name__, "zero", "avg"),
                       "%s: include_zero_slip is ignored for the average "
                       "displacement", type(self).__name__)
        return self._avg_disp(mag, style)[0]

    @staticmethod
    def _avg_disp(mag, style):
        """(AD'_P, AD'_P / D'_P(X/L = 0.25)), Eqs. 36-37."""
        i = STYLES.index(style)
        ratio_ad = float(np.exp(b_3[i] + b_4[i] * np.exp(b_5[i] * (mag - 5.0))))
        mu = _components(0.25, mag, style)["mu_prnc_prime"][0]
        return ratio_ad * float(mu ** (1 / 0.3)), ratio_ad

    def get_displ_max(self, mag, percentile=0.5, style=None):
        """Not available: LA23 gives the maximum displacement for the
        aggregate displacement only (Eqs. 38-39, Table 6)."""
        raise ValueError(
            f"{type(self).__name__}: LA23 gives the maximum displacement for "
            "the aggregate displacement only (Eqs. 38-39, Table 6); use "
            "Lavrentiadis2023PrimaryFD_aggregate.get_displ_max")
