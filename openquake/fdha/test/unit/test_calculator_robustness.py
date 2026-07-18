# -*- coding: utf-8 -*-
"""
Robustness guarantees of the branch calculator and the model adapter:
unknown names/parameters fail loudly, explicit zeros are honored, and
model errors propagate instead of silently zeroing hazard.
"""
import numpy as np
import pytest

from openquake.fdha.calc.calculators import (
    BaseFaultRuptureCalculator, MODEL_REGISTRY)
from openquake.fdha.calc.contexts import FDHAContext
from openquake.fdha.calc.model_adapter import LegacyModelAdapter

pytestmark = pytest.mark.unit


def test_registry_is_populated():
    assert 'Youngs2003PrimarySR' in MODEL_REGISTRY
    assert 'Petersen2011PrimaryFD' in MODEL_REGISTRY
    assert 'Visini2025SecondaryFD' in MODEL_REGISTRY


def test_unknown_model_name_rejected():
    with pytest.raises(ValueError, match="Unknown FDHA model class"):
        BaseFaultRuptureCalculator._instantiate_model({'type': 'NoSuchModel'})


def test_unknown_parameter_rejected():
    # 'stlye' is a typo of 'style': neither a constructor nor a get_prob
    # parameter, so it must fail instead of being silently dropped.
    with pytest.raises(ValueError, match="unknown parameter"):
        BaseFaultRuptureCalculator._instantiate_model({
            'type': 'Youngs2003PrimarySR',
            'parameters': {'stlye': 'all'},
        })


def test_known_call_time_parameter_accepted():
    model = BaseFaultRuptureCalculator._instantiate_model({
        'type': 'Youngs2003PrimarySR',
        'parameters': {'style': 'all'},
    })
    assert model.__class__.__name__ == 'Youngs2003PrimarySR'


def test_calc_param_honors_explicit_zero():
    calc = BaseFaultRuptureCalculator.__new__(BaseFaultRuptureCalculator)
    calc.config = {
        'calculation': {'r_threshold_km': 0},
        'parameters': {'r_threshold_km': 0.4},
    }
    # An explicit 0 in [calculation] must win over [parameters] and the
    # default (the old `or` chain treated 0 as missing).
    assert calc._calc_param('r_threshold_km', 0.1) == 0
    calc.config = {'calculation': {}, 'parameters': {}}
    assert calc._calc_param('r_threshold_km', 0.1) == 0.1


def _one_site_ctx(vs30=760.0):
    return FDHAContext(
        sids=np.array([0]), mag=np.array([6.5]), rake=np.array([0.0]),
        dip=np.array([90.0]), ztor=np.array([0.0]),
        occurrence_rate=np.array([0.01]), vs30=np.array([vs30]),
        r=np.array([0.0]), rx=np.array([0.0]),
        x_L=np.array([0.5]), L=np.array([20.0]),
    )


def test_model_error_propagates():
    class BoomPrimarySR:
        def get_prob(self, mag):
            raise ValueError('boom')

    adapter = LegacyModelAdapter(BoomPrimarySR(), {})
    with pytest.raises(RuntimeError, match='boom'):
        adapter.compute_primary_sr(_one_site_ctx(), {'method': 'mean'})


def test_model_returning_none_rejected():
    class NonePrimarySR:
        def get_prob(self, mag):
            return None

    adapter = LegacyModelAdapter(NonePrimarySR(), {})
    with pytest.raises(RuntimeError, match='returned'):
        adapter.compute_primary_sr(_one_site_ctx(), {'method': 'mean'})


def test_missing_vs30_becomes_none_for_model():
    """A NaN vs30 (site without a value, no reference_vs30_value) must reach
    a vs30-requiring model as None so its own validation fires, instead of
    NaN silently failing every threshold comparison."""
    seen = {}

    class Vs30PrimarySR:
        def get_prob(self, mag, vs30):
            seen['vs30'] = vs30
            if vs30 is None:
                raise ValueError('requires vs30')
            return 0.5

    adapter = LegacyModelAdapter(Vs30PrimarySR(), {})
    with pytest.raises(RuntimeError, match='requires vs30'):
        adapter.compute_primary_sr(
            _one_site_ctx(vs30=np.nan), {'method': 'mean'})
    assert seen['vs30'] is None
