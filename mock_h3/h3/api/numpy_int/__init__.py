"""
Stub for h3.api.numpy_int to provide geo_to_h3 / h3_to_geo when native h3 is missing.
"""

def geo_to_h3(lat, lon, res=5):
    return 0

def h3_to_geo(index):
    return (0.0, 0.0)

