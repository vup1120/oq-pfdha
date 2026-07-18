# -*- coding: utf-8 -*-
"""
Correctness: partial Chelungpu NRML - trace horizontal distance r must be
stable across rupture_mesh_spacing (decoupled from 3-D Rrup mesh density).

Uses the real MeshSpacingError partial fault geometry (not the full trace).
"""
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

from openquake.fdha.calc.utils.parsing import parse_source_model_faults
from openquake.fdha.calc.utils.rupture_distance import VectorizedRuptureDistanceCalculator

pytestmark = pytest.mark.unit

_MESH_DIR = Path(__file__).resolve().parents[1] / "fixtures/chelungpu_mesh"
_PARTIAL_XML = _MESH_DIR / "17_Chelungpu_fault_partial section.xml"


def _rupture_surface(spacing: float):
    sources = parse_source_model_faults(
        str(_PARTIAL_XML),
        hdf5path="",
        rupture_mesh_spacing=spacing,
        width_of_mfd_bin=0.1,
    )
    src = sources["17"]
    rup = next(src.iter_ruptures())
    return rup.surface


def _site_collection(lonlats):
    from openquake.hazardlib.site import Site, SiteCollection
    from openquake.hazardlib.geo import Point

    return SiteCollection(
        [Site(Point(float(lon), float(lat), 0.0)) for lon, lat in lonlats]
    )


@pytest.fixture(scope="module")
def partial_xml_exists():
    assert _PARTIAL_XML.is_file(), f"Missing fixture: {_PARTIAL_XML}"


def test_partial_chelungpu_trace_distance_invariant_fine_meshes(partial_xml_exists):
    """Fine rupture meshes (0.25–0.5 km): same trace polyline → r must match tightly.

    Sinuous real traces: spacing ≥~1 km changes the mesh top-row polyline enough that
    horizontal r can drift (geometric approximation), which is expected - not Rrup noise.
    """
    spacings = (0.25, 0.5)
    sites_ll = [
        (120.695, 23.985),
        (120.710, 24.005),
        (120.680, 23.970),
        (120.725, 24.050),
        (120.700, 23.920),
    ]
    sitecol = _site_collection(sites_ll)

    dist_by_spacing = []
    for sp in spacings:
        surf = _rupture_surface(sp)
        calc = VectorizedRuptureDistanceCalculator(sitecol, surf)
        dist_by_spacing.append(calc.calculate_site_to_trace_distances())

    assert_allclose(dist_by_spacing[1], dist_by_spacing[0], rtol=0.02, atol=0.02)
    assert np.all(dist_by_spacing[0] >= 0.0)
    assert np.all(np.isfinite(dist_by_spacing[0]))


def test_partial_chelungpu_1km_mesh_reasonable_vs_fine(partial_xml_exists):
    """1 km spacing vs 0.25 km: bounded drift on this partial geometry (curved trace)."""
    sites_ll = [
        (120.695, 23.985),
        (120.710, 24.005),
        (120.680, 23.970),
        (120.725, 24.050),
        (120.700, 23.920),
    ]
    sitecol = _site_collection(sites_ll)
    r_fine = VectorizedRuptureDistanceCalculator(
        sitecol, _rupture_surface(0.25)
    ).calculate_site_to_trace_distances()
    r_1km = VectorizedRuptureDistanceCalculator(
        sitecol, _rupture_surface(1.0)
    ).calculate_site_to_trace_distances()
    assert_allclose(r_1km, r_fine, rtol=0.08, atol=0.06)


def test_partial_chelungpu_on_trace_site_near_zero_r(partial_xml_exists):
    """A site on the shallow mesh trace should have r ~ 0 for any practical spacing."""
    surf = _rupture_surface(2.0)
    mesh = surf.mesh
    zero = np.isclose(mesh.depths, 0.0)
    row_counts = zero.sum(axis=1)
    col_counts = zero.sum(axis=0)
    if row_counts.max() >= col_counts.max():
        i = int(np.argmax(row_counts))
        lon0 = float(mesh.lons[i, zero[i, :].argmax()])
        lat0 = float(mesh.lats[i, zero[i, :].argmax()])
    else:
        j = int(np.argmax(col_counts))
        lon0 = float(mesh.lons[zero[:, j], j][0])
        lat0 = float(mesh.lats[zero[:, j], j][0])

    sc = _site_collection([(lon0, lat0)])
    for sp in (0.5, 3.0):
        surf2 = _rupture_surface(sp)
        r = VectorizedRuptureDistanceCalculator(sc, surf2).calculate_site_to_trace_distances()[0]
        assert_allclose(r, 0.0, atol=0.5)
