"""
Calculation utilities: distance calculators, probability reductions and
context re-exports.
"""

from .rupture_distance import (
    RuptureDistanceCalculator,
    VectorizedRuptureDistanceCalculator,
)
from .probability import _reduce_mc, _reduce_over_axes, _to_sites_x_displ

# Context classes need hazardlib; keep the subpackage importable without it.
try:
    from openquake.fdha.calc.contexts import FDHAContext, FDHAContextMaker
except ImportError:
    pass
