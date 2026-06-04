"""
H3 compatibility layer for OpenQuake with h3 v4.x

This module provides a compatibility layer to make h3 v4.x work with
OpenQuake code that expects h3 v3.x API.
"""

import sys
import logging

logger = logging.getLogger(__name__)


def patch_h3_for_openquake():
    """
    Monkey patch h3 v4.x to be compatible with OpenQuake's expected h3 v3.x API.
    Must be called before importing any OpenQuake modules that use h3.
    """
    try:
        import h3
        
        # Check if we have h3 v4.x (new API)
        if hasattr(h3, 'latlng_to_cell'):
            # Create compatibility module structure
            class MockNumpyInt:
                """Mock module to provide h3 v3 API using h3 v4 functions"""
                geo_to_h3 = h3.latlng_to_cell
                h3_to_geo = h3.cell_to_latlng
            
            class MockApi:
                """Mock api module"""
                numpy_int = MockNumpyInt()
            
            # Inject the compatibility layer
            h3.api = MockApi()
            
            # Also make the functions available at the expected import path
            sys.modules['h3.api'] = MockApi()
            sys.modules['h3.api.numpy_int'] = MockNumpyInt()
            
            logger.debug("H3 v4.x compatibility layer activated for OpenQuake")
            return True
        else:
            # h3 v3.x is already installed, no patching needed
            logger.debug("H3 v3.x detected, no compatibility layer needed")
            return False

    except ImportError:
        logger.warning("h3 library not installed")
        return False

# Apply the patch when this module is imported
patch_h3_for_openquake()
