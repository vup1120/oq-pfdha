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
Package :mod:`openquake.fdha.secondary_surf_displ` contains implementations of different
fault displacement prediction models
"""
from openquake.fdha.secondary_surf_displ.youngs2003 import Youngs2003SecondaryFD
from openquake.fdha.secondary_surf_displ.petersen2011 import Petersen2011SecondaryFD
from openquake.fdha.secondary_surf_displ.visini2025 import Visini2025SecondaryFD
from openquake.fdha.secondary_surf_displ.moss2022 import Moss2022SecondaryFD