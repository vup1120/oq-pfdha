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
Package :mod:`openquake.fdha.primary_surf_displ` contains implementations of different
fault displacement prediction models
"""
from openquake.fdha.primary_surf_displ.youngs2003 import Youngs2003PrimaryFD
from openquake.fdha.primary_surf_displ.petersen2011 import Petersen2011PrimaryFD
from openquake.fdha.primary_surf_displ.petersen2011 import Petersen2011PrimaryFD_bilinear
from openquake.fdha.primary_surf_displ.petersen2011 import Petersen2011PrimaryFD_elliptical
from openquake.fdha.primary_surf_displ.petersen2011 import Petersen2011PrimaryFD_quadratic
from openquake.fdha.primary_surf_displ.moss_ross2011 import MossRoss2011PrimaryFD
from openquake.fdha.primary_surf_displ.moss2024 import Moss2024PrimaryFD
from openquake.fdha.primary_surf_displ.moss2022 import Moss2022PrimaryFD
from openquake.fdha.primary_surf_displ.takao2013 import Takao2013PrimaryFD
from openquake.fdha.primary_surf_displ.lavrentiadis2023 import Lavrentiadis2023PrimaryFD
from openquake.fdha.primary_surf_displ.lavrentiadis2023 import Lavrentiadis2023PrimaryFD_principal
from openquake.fdha.primary_surf_displ.kuehn2024.kuehn2024 import Kuehn2024PrimaryFD
from openquake.fdha.primary_surf_displ.chiou2025 import Chiou2025PrimaryFD
