"""Compatibility imports for engine-owned scaling relationships."""
from openquake.hazardlib.scalerel.wc1994 import WC1994 as WellsCoppersmith1994
from openquake.hazardlib.scalerel.thingbaijam2017 import Thingbaijam2017
from openquake.hazardlib.scalerel.leonard2010 import Leonard2010

__all__ = ["WellsCoppersmith1994", "Thingbaijam2017", "Leonard2010"]
