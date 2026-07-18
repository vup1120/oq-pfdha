# -*- coding: utf-8 -*-
"""
Validation of the [parameters] primary_sr_reduction / secondary_sr_reduction
configs for hazard calculations.

These are undocumented developer knobs: jobs need not (and should not) set
them. Inside the hazard integral only 'mean' (the default; exact -
expectation commutes with the rate sum and the SR x FD product) and
'median' (legacy central-estimate heuristic) are legal. 'percentile' is
rejected loudly: a per-rupture quantile of exceedance probabilities is not
a fractile of any hazard distribution. The percentile branch of
utils.probability._reduce_mc is deliberately kept alive for non-integral
uses (e.g. a future scenario calculator), so these tests also pin that
behaviour.
"""
import numpy as np
import pytest

from openquake.fdha.calc.calculators import BaseFaultRuptureCalculator
from openquake.fdha.calc.utils.probability import _reduce_mc

pytestmark = pytest.mark.unit

validate = BaseFaultRuptureCalculator._validate_hazard_reduction


class TestValidateHazardReduction:
    def test_mean_is_accepted(self):
        assert validate('primary_sr_reduction', {'method': 'mean'}) == \
            {'method': 'mean'}

    def test_median_is_accepted_and_q_passes_through(self):
        cfg = validate('primary_sr_reduction', {'method': 'median', 'q': 50})
        assert cfg == {'method': 'median', 'q': 50}

    def test_method_is_case_normalized(self):
        assert validate('primary_sr_reduction', {'method': 'Mean'}) == \
            {'method': 'mean'}

    def test_bare_string_shorthand_is_accepted(self):
        assert validate('primary_sr_reduction', 'mean') == {'method': 'mean'}

    def test_percentile_is_rejected_with_guidance(self):
        with pytest.raises(ValueError) as exc:
            validate('primary_sr_reduction', {'method': 'percentile', 'q': 84})
        msg = str(exc.value)
        assert 'primary_sr_reduction' in msg
        assert 'does not produce a hazard fractile' in msg.replace('\n', ' ')
        assert '"mean"' in msg
        assert 'logic-tree realizations' in msg

    def test_percentile_is_rejected_for_secondary_too(self):
        with pytest.raises(ValueError, match='secondary_sr_reduction'):
            validate('secondary_sr_reduction', {'method': 'percentile'})

    def test_unknown_method_is_rejected(self):
        with pytest.raises(ValueError, match="unknown reduction method 'avg'"):
            validate('primary_sr_reduction', {'method': 'avg'})

    def test_non_dict_non_string_is_rejected(self):
        with pytest.raises(TypeError, match='primary_sr_reduction'):
            validate('primary_sr_reduction', [50, 84])

    def test_missing_method_defaults_to_mean(self):
        assert validate('primary_sr_reduction', {})['method'] == 'mean'


class TestPercentileMachineryKeptForFutureUse:
    """_reduce_mc must keep its percentile branch (scenario calculators)."""

    def test_reduce_mc_percentile_still_works(self):
        samples = np.linspace(0.0, 1.0, 101)
        assert _reduce_mc(samples, method='percentile', q=84) == \
            pytest.approx(0.84)

    def test_reduce_mc_percentile_fractional_q(self):
        samples = np.linspace(0.0, 1.0, 101)
        assert _reduce_mc(samples, method='percentile', q=0.84) == \
            pytest.approx(0.84)
