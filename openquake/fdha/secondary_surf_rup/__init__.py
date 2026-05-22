# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#
# Copyright (C) 2012-2023 GEM Foundation
#
# OpenQuake is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenQuake is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenQuake. If not, see <http://www.gnu.org/licenses/>.

"""
Package :mod:`openquake.fdha.secondary_surf_rup` contains implementations of different
surface rupture prediction models
"""

from openquake.fdha.secondary_surf_rup.youngs2003 import Youngs2003SecondarySR
from openquake.fdha.secondary_surf_rup.petersen2011 import Petersen2011SecondarySR
from openquake.fdha.secondary_surf_rup.petersen2011 import Petersen2011SecondarySR_default
from openquake.fdha.secondary_surf_rup.visini2025 import Visini2025SecondarySR
from openquake.fdha.secondary_surf_rup.takao2014 import Takao2014SecondarySR
from openquake.fdha.secondary_surf_rup.ferrario2021 import FerrarioLivio2021SecondarySR
from openquake.fdha.secondary_surf_rup.takao2013 import Takao2013SecondarySR
from openquake.fdha.secondary_surf_rup.rodriguez2023 import Rodriguez2023SecondarySR
from openquake.fdha.secondary_surf_rup.fixed import FixedSecondarySR
from openquake.fdha.secondary_surf_rup.moss2022 import Moss2022SecondarySR