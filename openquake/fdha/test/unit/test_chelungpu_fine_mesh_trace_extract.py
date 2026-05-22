# -*- coding: utf-8 -*-
"""Regression: fine rupture_mesh_spacing must not collapse trace to one point."""
from __future__ import annotations

from pathlib import Path

import pytest

from openquake.fdha.calc.utils.parsing import parse_source_model_faults
from openquake.fdha.calc.utils.rupture_distance import _extract_fault_trace_from_mesh

_MESH_CASE = Path(
    __file__
).resolve().parents[1] / "fixtures/chelungpu_mesh/17_Chelungpu_fault.xml"


@pytest.mark.skipif(not _MESH_CASE.is_file(), reason="Chelungpu XML not in tree")
@pytest.mark.slow
def test_trace_extract_002_km_has_many_vertices_not_collapsed_to_one():
    src = parse_source_model_faults(
        str(_MESH_CASE), rupture_mesh_spacing=0.02, width_of_mfd_bin=0.1
    )
    s = next(iter(src.values()))
    rup = next(s.iter_ruptures())
    tr = _extract_fault_trace_from_mesh(rup.surface)
    assert tr.shape[0] >= 50, tr.shape
    assert tr.shape[1] == 2
    dlon = float(tr[:, 0].max() - tr[:, 0].min())
    dlat = float(tr[:, 1].max() - tr[:, 1].min())
    assert dlon > 0.01 or dlat > 0.01, (dlon, dlat)
