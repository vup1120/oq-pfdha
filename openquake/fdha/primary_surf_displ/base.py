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
Module :mod:`openquake.fdha.primary_surf_displ.base` defines abstract base
classes for :class:`BasePrimarySurfDispl <BasePrimarySurfDispl>` and
 :class:`BaseSecondarySurfDispl <BaseSecondarySurfDispl>`  
"""

import abc


class BasePrimarySurfDispl(metaclass=abc.ABCMeta):
    """Abstract base class for principal (primary) fault-displacement models.

    Subclasses implement :meth:`get_prob`, returning the probability that the
    principal displacement exceeds a given value (in metres).
    """

    #: Reference-line treatment this model needs when the source has no
    #: continuous fault trace (multiFaultSource / kite sections); one of
    #: 'lcp', 'ecs', 'segments'. Declarative, mirroring hazardlib's
    #: REQUIRES_DISTANCES pattern: the FDHA context maker computes the union
    #: of declared requirements once per rupture. Irrelevant for single-strand
    #: sources, whose trace is used directly.
    MULTIFAULT_REFERENCE_LINE = "lcp"

    @abc.abstractmethod
    def get_prob(self):
        """
        Return the probability that the primary displacement will exceed
        a certain value [m]
        """

    def __str__(self):
        """
        Returns the name of the class
        """
        return self.__class__.__name__

    def __repr__(self):
        """
        Returns the name of the class in angular brackets
        """
        return "<%s>" % self.__class__.__name__


class BaseSecondarySurfDispl(metaclass=abc.ABCMeta):
    """Abstract base class for distributed (secondary) fault-displacement models.

    Subclasses implement :meth:`get_prob`, returning the probability that the
    distributed displacement exceeds a given value (in metres).
    """

    #: See BasePrimarySurfDispl.MULTIFAULT_REFERENCE_LINE.
    MULTIFAULT_REFERENCE_LINE = "lcp"

    #: Near-field regularisation for the distributed *displacement* evaluation,
    #: applied at the calc/ adapter boundary so model files stay paper-faithful.
    #: ``None`` (default) means no floor: bounded models (Takao exponential,
    #: Visini's own 5 m exclusion) leave this alone. ``"footprint_half"`` clamps
    #: the distance fed to the displacement regression to ``max(r, z/2)`` with
    #: ``z = site_footprint_m/1000`` km -- a tool regularisation of Petersen
    #: (2011) eq.18's r -> 0 divergence (docs/design/
    #: rupture_location_uncertainty.md, decision D7).
    NEAR_FIELD_FLOOR = None

    @abc.abstractmethod
    def get_prob(self):
        """
        Return the probability that the secondary displacement will exceed
        a certain value [m]
        """

    def __str__(self):
        """
        Returns the name of the class
        """
        return self.__class__.__name__

    def __repr__(self):
        """
        Returns the name of the class in angular brackets
        """
        return "<%s>" % self.__class__.__name__