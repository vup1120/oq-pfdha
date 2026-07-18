# LCP validation fixtures

Ground-truth fixtures for `openquake/fdha/calc/utils/lcp.py`, the pure-Python
port of the Thomas/Milliner Least Cost Path reference-line workflow
(`test/usercase/LCP/Thomas_LCP/Python_LCP/pyLCP/*.ipynb`, which uses
GDAL + `skimage.graph.route_through_array(geometric=True,
fully_connected=True)`).

## ridgecrest/

Copied **verbatim** (`cp`, not transcribed) from the Thomas notebook bundle
(2019 Ridgecrest rupture, UTM zone 11N / EPSG:32611):

| file | role |
|---|---|
| `Images_4_LCP_Ridge_main_Ridgecrest_PL_Cost_Raster_4_LCP.tif` | input cost raster, "PL" case (547x435, 87 m pixels, float32 imagery-derived costs) |
| `Images_4_LCP_Ridge_main_Ridgecrest_S2_Cost_Raster_4_LCP.tif` | input cost raster, "S2" case (799x684, 90 m pixels) |
| `Ridgecrest_{PL,S2}_LCP_{start,end}_point*.txt` | picked LCP endpoints (tab-delimited UTM x, y) |
| `LCP_path_{PL,S2}.txt` | **ground truth**: notebook output path, columns = UTMx UTMy lon lat, one row per path pixel (corner coordinates, `Pixel2Map` convention) |

`test/unit/test_lcp_reference_values.py` reads the GeoTIFFs with Pillow
(test-only dependency; the runtime module needs numpy/scipy/pyproj only),
routes the LCP with `lcp.route_through_raster`, and requires a
**node-for-node exact match** against `LCP_path_*.txt` (verified 2026-07-02:
max |dxy| = 0.0 on both cases, identical path cost).

The notebook bundle also contains `LCP_path_{PL,S2}_orig.txt`; those were
produced from an *earlier* version of the cost rasters and do NOT match the
committed rasters (Hausdorff distance km-scale) - deliberately not fixtures.
The Hector Mine "field map" case is not used either: its cost raster is built
by GDAL shapefile rasterization, which is not bit-reproducible without GDAL;
our trace rasterization is validated by synthetic unit tests instead
(`test_lcp_unit.py`).
