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
Module :mod:`openquake.fdha.secondary_surf_displ.nishizaka2026` implements
the far-field surface-rupture displacement model of Nishizaka et al. (2026)
in :class:`Nishizaka2026SecondaryFD`.

Supported fault styles: strike-slip only. Coefficients are available for the
2016 Kumamoto earthquake dataset only (the 2019 Ridgecrest far-field
displacement data were too sparse to regress; see the paper).

Following Takao et al. (2016) and Equations 3-4 of Nishizaka et al. (2026),
the normalized displacement ``x = SRD / EAD`` (surface rupture displacement
over the average displacement of the earthquake source faults) at a distance
``r1`` from the source fault is gamma distributed:

    x ~ Gamma(shape = a, scale = b(r1)),  b(r1) = (c6 / c8) * exp(c7 * r1)

with coefficients from Table 2 of the paper, regressed for all data and for
the proximal (r2 <= 1 km) and non-proximal (r2 > 1 km) subsets, where ``r2``
is the shortest surface distance from mapped, pre-existing active faults.

Reference
---------
Nishizaka, N., Onishi, K., Ikeda, M., Si, H., Yamamoto, K., & Tsuji, T.
(2026). Characteristics of far-field surface ruptures caused by two recent
strike-slip earthquakes: insights into fault displacement prediction.
Seismological Research Letters. https://doi.org/10.1785/0220250293
"""

import numpy as np
from scipy.stats import gamma
from openquake.fdha.primary_surf_displ.base import BaseSecondarySurfDispl


class Nishizaka2026SecondaryFD(BaseSecondarySurfDispl):
    """
    Far-field displacement exceedance model of Nishizaka et al. (2026) for
    strike-slip earthquakes (2016 Kumamoto regression), conditional on a
    surface rupture occurring at the site.

    Notes
    -----
    - The model predicts displacement normalized by the average displacement
      of the earthquake source faults (EAD, in metres), which must be
      supplied by the user — e.g. from a fault model or an empirical
      magnitude scaling relation. For the 2016 Kumamoto earthquake the paper
      uses EAD = 1.9 m (Asano and Iwata, 2021).
    - ``dataset`` selects the regression subset: ``"proximal"``
      (r2 <= 1 km from mapped, pre-existing active faults),
      ``"nonproximal"`` (r2 > 1 km), or ``"all"``. Alternatively, pass the
      site's ``r2`` distance (km) and the subset is chosen automatically.
    - Like the companion occurrence model, the displacement model is
      independent of earthquake magnitude (magnitude enters only through
      EAD).
    """

    #: Table 2 of Nishizaka et al. (2026): gamma shape ``a`` and regression
    #: coefficients (c6, c7, c8); distances in km.
    COEFFS = {
        "all": {"a": 1.00, "c6": 0.808, "c7": -0.252, "c8": 2.30},
        "proximal": {"a": 1.06, "c6": 0.804, "c7": -0.232, "c8": 2.27},
        "nonproximal": {"a": 1.07, "c6": 0.115, "c7": -0.0903, "c8": 2.27},
    }

    #: r2 threshold (km) separating the proximal and non-proximal subsets.
    R2_THRESHOLD_KM = 1.0

    def _select_dataset(self, dataset, r2):
        if dataset is not None:
            ds = str(dataset).strip().lower().replace("-", "").replace("_", "")
            if ds not in self.COEFFS:
                raise ValueError(
                    f"Invalid dataset '{dataset}'. "
                    f"Accepted: {sorted(self.COEFFS.keys())}")
            return ds
        if r2 is None:
            raise ValueError("Provide either 'dataset' or 'r2'")
        return ("proximal" if float(r2) <= self.R2_THRESHOLD_KM
                else "nonproximal")

    def _gamma_params(self, r, coeffs):
        r1_arr = np.asarray(r, dtype=float)
        if np.any(r1_arr < 0):
            raise ValueError("Distance r must be non-negative")
        scale = (coeffs["c6"] / coeffs["c8"]) * np.exp(coeffs["c7"] * r1_arr)
        return coeffs["a"], scale

    def get_prob(self, d, r, ead, dataset=None, r2=None):
        """
        Probability that the surface-rupture displacement exceeds ``d``,
        conditional on a surface rupture occurring at distance ``r`` from
        the earthquake source fault.

        :param d:
            Target displacement threshold(s) in metres, scalar or
            array-like of shape (n_displacements,).
        :param r:
            Distance r1 from the earthquake source fault in km, scalar or
            array-like of shape (n_sites,).
        :param ead:
            Average displacement of the earthquake source faults in metres
            (scalar, > 0). E.g. 1.9 m for the 2016 Kumamoto earthquake.
        :param dataset:
            ``"proximal"``, ``"nonproximal"``, or ``"all"`` (Table 2
            regression subsets). May be omitted if ``r2`` is given.
        :param r2:
            Optional shortest surface distance (km) from mapped,
            pre-existing active faults; selects ``"proximal"`` when
            r2 <= 1 km, ``"nonproximal"`` otherwise. Ignored when
            ``dataset`` is given.
        :return:
            Exceedance probability with shape (n_sites, n_displacements)
            (squeezed to scalar/1D for scalar inputs).
        """
        ds = self._select_dataset(dataset, r2)
        coeffs = self.COEFFS[ds]

        ead = float(ead)
        if ead <= 0.0:
            raise ValueError(f"ead must be positive, got {ead}")

        d_arr = np.atleast_1d(np.asarray(d, dtype=float))
        if np.any(d_arr < 0):
            raise ValueError("Displacement thresholds must be non-negative")
        r_arr = np.atleast_1d(np.asarray(r, dtype=float))

        shape_a, scale = self._gamma_params(r_arr, coeffs)
        x = d_arr[np.newaxis, :] / ead                # (1, n_d)
        scale_b = scale[:, np.newaxis]                # (n_sites, 1)
        prob = gamma.sf(x, shape_a, scale=scale_b)    # (n_sites, n_d)

        if prob.size == 1:
            return float(prob[0, 0])
        return prob.squeeze()

    def get_percentile_displacement(self, p, r, ead, dataset=None, r2=None):
        """
        Displacement (metres) at the ``p``-th percentile of the gamma
        distribution at distance ``r`` (km) — e.g. ``p=90`` reproduces the
        90th-percentile curves of Figure 6 of the paper when divided
        by ``ead``.

        :param p: Percentile in (0, 100).
        :param r: Distance r1 from the source fault in km (scalar or array).
        :param ead: Average displacement of the source faults in metres.
        :param dataset: See :meth:`get_prob`.
        :param r2: See :meth:`get_prob`.
        :return: Displacement in metres (scalar or array following ``r``).
        """
        p = float(p)
        if not 0.0 < p < 100.0:
            raise ValueError(f"Percentile must be in (0, 100), got {p}")
        ds = self._select_dataset(dataset, r2)
        shape_a, scale = self._gamma_params(r, self.COEFFS[ds])
        x_p = gamma.ppf(p / 100.0, shape_a, scale=scale)
        out = x_p * float(ead)
        return out.item() if np.ndim(out) == 0 else out
