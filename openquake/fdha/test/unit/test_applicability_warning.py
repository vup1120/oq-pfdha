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


def test_gaussian_zero_weight_excludes_warning_in_model_metric(caplog):
    """An on-trace canonical distance masks even a far model-specific
    distance; a nearby site with positive G must still be reported."""
    ctx = _ctx([0.0, 0.01])
    ctx.metrics_for = lambda method: (np.array([3.0, 3.0]), ctx.x_L, ctx.L)
    tracker = ApplicabilityTracker(r_threshold_km=0.1, r_sigma_km=0.05)
    tracker.observe(Petersen2011SecondaryFD(), ctx)
    with caplog.at_level(logging.WARNING, logger="openquake.fdha.calc.hazard"):
        tracker.emit()
    assert len(caplog.records) == 1
    assert "1 site(s)" in caplog.records[0].getMessage()


def test_declared_exclusion_keeps_its_own_metric_and_closed_outer_edge():
    from openquake.fdha.calc.hazard import _inside_declared_exclusion
    ctx = _ctx([0.1, 0.1, 0.1])
    ctx.metrics_for = lambda method: (
        np.array([0.0, 0.004, 0.005]), ctx.x_L, ctx.L)
    np.testing.assert_array_equal(
        _inside_declared_exclusion(Visini2025SecondaryFD(), ctx),
        [True, True, False])


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


def test_visini_sub_5m_site_is_excluded_not_extrapolated(caplog):
    """A site inside Visini's declared 5 m exclusion is never an
    *extrapolation*: on either W_p path the distributed term is not
    evaluated there. At sigma = 0 the complementary split masks it; at
    sigma > 0 the Gaussian complement is positive off-trace, so the kernel
    gates it off (``_inside_declared_exclusion``) and the site carries the
    principal contribution only. Neither case may report extrapolation."""
    model = Visini2025SecondaryFD()

    for r_sigma in (0.0, 0.05):
        tracker = ApplicabilityTracker(r_threshold_km=0.1,
                                       r_sigma_km=r_sigma)
        tracker.observe(model, _ctx([0.001]))  # 1 m, inside the 5 m floor
        caplog.clear()
        with caplog.at_level(logging.WARNING,
                             logger="openquake.fdha.calc.hazard"):
            tracker.emit()
        assert not [r for r in caplog.records
                    if "extrapolating" in r.getMessage()], (
            f"sub-5 m site reported as extrapolation at sigma={r_sigma}")


def test_gated_sites_are_reported_as_an_exclusion(caplog):
    """Gating must not be silent: the run says which sites lost their
    distributed term and why."""
    model = Visini2025SecondaryFD()
    tracker = ApplicabilityTracker(r_threshold_km=0.1, r_sigma_km=0.05)
    tracker.note_excluded(model, _ctx([0.001]), np.array([True]))
    with caplog.at_level(logging.WARNING, logger="openquake.fdha.calc.hazard"):
        tracker.emit()
    msgs = [r.getMessage() for r in caplog.records
            if "Visini2025SecondaryFD" in r.getMessage()]
    assert len(msgs) == 1
    assert "1 site(s)" in msgs[0]
    assert "near-trace exclusion" in msgs[0]
    assert "principal contribution only" in msgs[0]


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
