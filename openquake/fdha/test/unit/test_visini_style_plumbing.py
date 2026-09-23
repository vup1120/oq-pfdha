# -*- coding: utf-8 -*-
"""
Regression tests for the faulting-style plumbing of the Visini (2025)
secondary pipeline.

The Visini SR/FD coefficients are style-specific. The hazard pipeline must
resolve the style like every other model path does (explicit model parameter
wins, otherwise derived from the rupture rake) - it must NEVER silently fall
back to 'normal' on a reverse fault: that shifts the FD median by ~1.7x
(the e = -0.5259 style term) and swaps the SR occurrence tables.
"""
import numpy as np
import pytest
from types import SimpleNamespace

from openquake.fdha.calc.contexts import FDHAContext
from openquake.fdha.calc.hazard import _compute_rupture_contribution
from openquake.fdha.calc.visini import VisiniSecondaryCalculator

pytestmark = [pytest.mark.unit, pytest.mark.visini2025]


class RecorderSR:
    """Stub Visini SR model recording the style it is called with."""
    MULTIFAULT_REFERENCE_LINE = "segments"

    def __init__(self):
        self.styles = []

    def calculate_rank2_total_probability_vectorized(
            self, mag, r_array, rx_array, style, **kwargs):
        self.styles.append(style)
        n = len(np.atleast_1d(r_array))
        return {"P_slice": np.full(n, 0.5),
                "P_along_strike": np.full(n, 0.5),
                "P_total": np.full(n, 0.25)}

    def get_prob(self, mag, r, rx, style, pixel_size, combination):
        # combination-C path
        self.styles.append(style)
        return np.full(len(np.atleast_1d(r)), 0.1)


class RecorderFD:
    """Stub Visini FD model recording the style it is called with."""
    MULTIFAULT_REFERENCE_LINE = "segments"

    def __init__(self):
        self.styles = []

    def get_prob(self, d, mag, s, rx, X_L_ratio=None, dip=None,
                 combination=None, style=None, **kwargs):
        self.styles.append(style)
        n = len(np.atleast_1d(s))
        return np.zeros((n, len(np.atleast_1d(d))))


def _ctx(rake):
    n = 1
    return FDHAContext(
        sids=np.array([0]),
        mag=np.full(n, 6.5),
        rake=np.full(n, rake),
        dip=np.full(n, 45.0),
        ztor=np.zeros(n),
        occurrence_rate=np.full(n, 1e-4),
        vs30=np.full(n, 760.0),
        r=np.full(n, 1.0),
        rx=np.full(n, 1.0),
        x_L=np.full(n, 0.5),
        L=np.full(n, 20.0),
        lons=np.full(n, 13.0),
        lats=np.full(n, 42.0),
    )


def _run(rake, model_params=None, case_label="case3"):
    sr, fd = RecorderSR(), RecorderFD()
    params = dict(model_params or {})
    visini_calc = VisiniSecondaryCalculator(
        base_sec_rup_params={k: v for k, v in params.items()},
        base_sec_displ_params={},
        case_label=case_label,
    )
    calculator = SimpleNamespace(
        secondary_surf_rup_model=sr,
        secondary_surf_displ_model=fd,
        get_model_parameters=lambda name: (
            dict(params) if name == 'secondary_surf_rup' else {}),
    )
    _compute_rupture_contribution(
        ctx=_ctx(rake),
        adapters={},
        target_displacements=np.array([0.01, 0.1]),
        p_sr_red_cfg={},
        s_sr_red_cfg={},
        r_threshold_km=0.1,
        use_visini=True,
        visini_calc=visini_calc,
        calculator=calculator,
    )
    return sr, fd


def test_reverse_rake_derives_reverse_style():
    """No explicit style: a reverse-rake rupture must use 'reverse'."""
    sr, fd = _run(rake=90.0)
    assert sr.styles == ["reverse"]
    assert fd.styles == ["reverse"]


def test_normal_rake_derives_normal_style():
    sr, fd = _run(rake=-90.0)
    assert sr.styles == ["normal"]
    assert fd.styles == ["normal"]


def test_explicit_style_parameter_wins_over_rake():
    """A style set on the model branch overrides the rake-derived one."""
    sr, fd = _run(rake=-90.0, model_params={"style": "reverse"})
    assert sr.styles == ["reverse"]
    assert fd.styles == ["reverse"]


def test_combination_c_receives_the_same_style():
    """case1 activates combination C, whose direct get_prob call previously
    dropped the style unless it was an explicit model parameter."""
    sr, fd = _run(rake=90.0, case_label="case1")
    # combos A, B (vectorized) + C (direct get_prob): all 'reverse'
    assert set(sr.styles) == {"reverse"}
    assert len(sr.styles) == 3
    assert set(fd.styles) == {"reverse"}


def test_strike_slip_rake_fails_loudly_with_real_models():
    """Visini is a dip-slip model: a strike-slip rupture without an explicit
    style must raise, not silently compute with 'normal' coefficients."""
    from openquake.pfd.secondary_surf_rup.visini2025 import (
        Visini2025SecondarySR)
    from openquake.pfd.secondary_surf_displ.visini2025 import (
        Visini2025SecondaryFD)

    visini_calc = VisiniSecondaryCalculator(
        base_sec_rup_params={}, base_sec_displ_params={},
        case_label="case3")
    calculator = SimpleNamespace(
        secondary_surf_rup_model=Visini2025SecondarySR(),
        secondary_surf_displ_model=Visini2025SecondaryFD(),
        get_model_parameters=lambda name: {},
    )
    with pytest.raises(ValueError, match="normal.*reverse|reverse.*normal"):
        _compute_rupture_contribution(
            ctx=_ctx(rake=0.0),  # strike-slip
            adapters={},
            target_displacements=np.array([0.1]),
            p_sr_red_cfg={},
            s_sr_red_cfg={},
            r_threshold_km=0.1,
            use_visini=True,
            visini_calc=visini_calc,
            calculator=calculator,
        )
