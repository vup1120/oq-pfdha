# -*- coding: utf-8 -*-
"""
C4 applicability advisory: distributed FD models declare their calibrated
distance range (APPLICABILITY_RANGE, in their own metric); at run time sites
evaluated outside the range are reported by ONE logging.warning per model
per run, carrying the offending-site count. No behaviour change.
"""
import logging

import numpy as np
import pytest

from openquake.fdha.calc.contexts import FDHAContext
from openquake.fdha.calc.hazard import ApplicabilityTracker
from openquake.fdha.secondary_surf_displ.petersen2011 import (
    Petersen2011SecondaryFD)
from openquake.fdha.secondary_surf_displ.visini2025 import (
    Visini2025SecondaryFD)

pytestmark = pytest.mark.unit


def _ctx(r_km, rx_km=None, sids=None):
    r = np.asarray(r_km, dtype=float)
    n = len(r)
    return FDHAContext(
        sids=np.arange(n) if sids is None else np.asarray(sids),
        mag=np.full(n, 7.0),
        rake=np.full(n, 90.0),
        dip=np.full(n, 45.0),
        ztor=np.zeros(n),
        occurrence_rate=np.full(n, 1e-4),
        vs30=np.full(n, 760.0),
        r=r,
        rx=r.copy() if rx_km is None else np.asarray(rx_km, dtype=float),
        x_L=np.full(n, 0.5),
        L=np.full(n, 40.0),
    )


def test_petersen_beyond_2km_warns_once_with_site_count(caplog):
    tracker = ApplicabilityTracker(r_threshold_km=0.1, r_sigma_km=0.0)
    model = Petersen2011SecondaryFD()
    # 0.5 km inside range; 3.0 and 5.0 km beyond the 2 km data limit.
    tracker.observe(model, _ctx([0.5, 3.0, 5.0]))
    # Second rupture revisits the same sites: the count must not inflate.
    tracker.observe(model, _ctx([0.4, 2.5, 4.0]))
    with caplog.at_level(logging.WARNING, logger="openquake.fdha.calc.hazard"):
        tracker.emit()
    warnings = [rec for rec in caplog.records
                if "Petersen2011SecondaryFD" in rec.getMessage()]
    assert len(warnings) == 1  # once per model per run
    msg = warnings[0].getMessage()
    assert "2 site(s)" in msg
    assert "applicability" in msg
    assert "2 km" in msg  # the declared source is named


def test_inside_range_no_warning(caplog):
    tracker = ApplicabilityTracker(r_threshold_km=0.1, r_sigma_km=0.0)
    tracker.observe(Petersen2011SecondaryFD(), _ctx([0.3, 1.0, 1.9]))
    with caplog.at_level(logging.WARNING, logger="openquake.fdha.calc.hazard"):
        tracker.emit()
    assert not caplog.records


def test_no_declared_range_no_warning(caplog):
    from openquake.fdha.secondary_surf_displ.moss2022 import (
        Moss2022SecondaryFD)
    tracker = ApplicabilityTracker(r_threshold_km=0.1, r_sigma_km=0.0)
    tracker.observe(Moss2022SecondaryFD(), _ctx([50.0, 100.0]))
    tracker.observe(None, _ctx([50.0]))
    with caplog.at_level(logging.WARNING, logger="openquake.fdha.calc.hazard"):
        tracker.emit()
    assert not caplog.records


def test_visini_hw_fw_outer_edges(caplog):
    """rx >= 0 (hanging wall) is valid to 10 km, rx < 0 (footwall) only to
    8 km: the same 9 km distance offends on the footwall only."""
    tracker = ApplicabilityTracker(r_threshold_km=0.1, r_sigma_km=0.0)
    model = Visini2025SecondaryFD()
    tracker.observe(model, _ctx([9.0, 9.0], rx_km=[9.0, -9.0]))
    with caplog.at_level(logging.WARNING, logger="openquake.fdha.calc.hazard"):
        tracker.emit()
    warnings = [rec for rec in caplog.records
                if "Visini2025SecondaryFD" in rec.getMessage()]
    assert len(warnings) == 1
    assert "1 site(s)" in warnings[0].getMessage()


def test_visini_5m_floor_masked_by_boxcar_at_sigma0(caplog):
    """At sigma = 0 the distributed term is complementary-masked inside
    |r| <= r_threshold_km, so a 1 m site is NOT an extrapolation (the model
    is never evaluated with weight there); at sigma > 0 the additive path
    does evaluate it and the sub-5 m floor must warn."""
    model = Visini2025SecondaryFD()

    tracker0 = ApplicabilityTracker(r_threshold_km=0.1, r_sigma_km=0.0)
    tracker0.observe(model, _ctx([0.001]))  # 1 m, inside the boxcar
    with caplog.at_level(logging.WARNING, logger="openquake.fdha.calc.hazard"):
        tracker0.emit()
    assert not caplog.records

    tracker1 = ApplicabilityTracker(r_threshold_km=0.1, r_sigma_km=0.05)
    tracker1.observe(model, _ctx([0.001]))
    with caplog.at_level(logging.WARNING, logger="openquake.fdha.calc.hazard"):
        tracker1.emit()
    warnings = [rec for rec in caplog.records
                if "Visini2025SecondaryFD" in rec.getMessage()]
    assert len(warnings) == 1
    assert "1 site(s)" in warnings[0].getMessage()


def test_offending_sites_accumulate_across_ruptures(caplog):
    """The count is the union of offending site ids across the whole run."""
    tracker = ApplicabilityTracker(r_threshold_km=0.1, r_sigma_km=0.0)
    model = Petersen2011SecondaryFD()
    tracker.observe(model, _ctx([3.0], sids=[7]))
    tracker.observe(model, _ctx([4.0], sids=[9]))
    tracker.observe(model, _ctx([5.0], sids=[7]))  # same site again
    with caplog.at_level(logging.WARNING, logger="openquake.fdha.calc.hazard"):
        tracker.emit()
    warnings = [rec for rec in caplog.records
                if "Petersen2011SecondaryFD" in rec.getMessage()]
    assert len(warnings) == 1
    assert "2 site(s)" in warnings[0].getMessage()
