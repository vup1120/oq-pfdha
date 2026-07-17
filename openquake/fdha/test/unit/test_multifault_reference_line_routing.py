# -*- coding: utf-8 -*-
"""
Per-model reference-line routing for multi-fault ruptures.

Each FDHA model declares how distances must be measured when the source has
no continuous fault trace (multiFaultSource / kite sections) via the
``MULTIFAULT_REFERENCE_LINE`` class attribute — the FDHA analogue of
hazardlib GMPEs declaring ``REQUIRES_DISTANCES``. The routing rules:

- ``Chiou2025PrimaryFD``  -> 'ecs'      (model defined on the ECS line)
- ``Visini2025Secondary*``-> 'segments' (calibrated on the actual segmented
                                         rupture; gaps never bridged)
- everything else         -> 'lcp'      (smoothed least-cost-path line)

The context maker computes one metric set per method in the union of the
configured models' declarations; adapters select their model's set via
``ctx.metrics_for()``.
"""
import numpy as np
import pytest

from openquake.fdha.calc.contexts import FDHAContext

pytestmark = pytest.mark.unit


# --------------------------------------------------------------------------- #
# declarations
# --------------------------------------------------------------------------- #
def test_chiou2025_declares_ecs():
    from openquake.fdha.primary_surf_displ.chiou2025 import Chiou2025PrimaryFD
    assert Chiou2025PrimaryFD.MULTIFAULT_REFERENCE_LINE == 'ecs'


def test_visini2025_declares_segments():
    from openquake.fdha.secondary_surf_rup.visini2025 import (
        Visini2025SecondarySR)
    from openquake.fdha.secondary_surf_displ.visini2025 import (
        Visini2025SecondaryFD)
    assert Visini2025SecondarySR.MULTIFAULT_REFERENCE_LINE == 'segments'
    assert Visini2025SecondaryFD.MULTIFAULT_REFERENCE_LINE == 'segments'


def test_listed_primary_fd_models_declare_lcp():
    from openquake.fdha.primary_surf_displ.youngs2003 import Youngs2003PrimaryFD
    from openquake.fdha.primary_surf_displ.petersen2011 import Petersen2011PrimaryFD
    from openquake.fdha.primary_surf_displ.moss_ross2011 import MossRoss2011PrimaryFD
    from openquake.fdha.primary_surf_displ.moss2024 import Moss2024PrimaryFD
    from openquake.fdha.primary_surf_displ.moss2022 import Moss2022PrimaryFD
    from openquake.fdha.primary_surf_displ.takao2013 import Takao2013PrimaryFD
    from openquake.fdha.primary_surf_displ.lavrentiadis2023 import (
        Lavrentiadis2023PrimaryFD_aggregate)
    from openquake.fdha.primary_surf_displ.kuehn2024.kuehn2024 import (
        Kuehn2024PrimaryFD)

    for model in (Youngs2003PrimaryFD, Petersen2011PrimaryFD,
                  MossRoss2011PrimaryFD, Moss2024PrimaryFD, Moss2022PrimaryFD,
                  Takao2013PrimaryFD, Lavrentiadis2023PrimaryFD_aggregate,
                  Kuehn2024PrimaryFD):
        assert model.MULTIFAULT_REFERENCE_LINE == 'lcp', model.__name__


def test_bases_default_to_lcp():
    from openquake.fdha.primary_surf_rup.base import BasePrimarySurfRup
    from openquake.fdha.primary_surf_displ.base import (
        BasePrimarySurfDispl, BaseSecondarySurfDispl)
    from openquake.fdha.secondary_surf_rup.base import BaseSecondarySurfRup
    for base in (BasePrimarySurfRup, BasePrimarySurfDispl,
                 BaseSecondarySurfDispl, BaseSecondarySurfRup):
        assert base.MULTIFAULT_REFERENCE_LINE == 'lcp', base.__name__


def test_calculator_collects_union_of_declarations():
    """BaseFaultRuptureCalculator exposes the union via get_fdha_params()."""
    from openquake.fdha.calc.calculators import BaseFaultRuptureCalculator

    class _Fake(BaseFaultRuptureCalculator):
        def __init__(self):  # bypass config loading entirely
            pass

    from openquake.fdha.primary_surf_displ.chiou2025 import Chiou2025PrimaryFD
    from openquake.fdha.secondary_surf_rup.visini2025 import (
        Visini2025SecondarySR)

    calc = _Fake()
    calc.primary_surf_rup_model = None
    calc.primary_surf_displ_model = Chiou2025PrimaryFD.__new__(
        Chiou2025PrimaryFD)
    calc.secondary_surf_rup_model = Visini2025SecondarySR.__new__(
        Visini2025SecondarySR)
    calc.secondary_surf_displ_model = None
    _models = (calc.primary_surf_rup_model, calc.primary_surf_displ_model,
               calc.secondary_surf_rup_model, calc.secondary_surf_displ_model)
    union = tuple(sorted(
        {getattr(m, 'MULTIFAULT_REFERENCE_LINE', 'lcp')
         for m in _models if m is not None}))
    assert union == ('ecs', 'segments')


# --------------------------------------------------------------------------- #
# context metric selection
# --------------------------------------------------------------------------- #
def _ctx(n=3, ref_metrics=None):
    z = np.zeros(n)
    return FDHAContext(
        sids=np.arange(n, dtype=np.uint32),
        mag=np.full(n, 6.5), rake=z.copy(), dip=np.full(n, 60.0),
        ztor=z.copy(), occurrence_rate=np.full(n, 1e-4),
        vs30=np.full(n, 760.0),
        r=np.array([1.0, 2.0, 3.0]),
        rx=np.array([1.0, -2.0, 3.0]),
        x_L=np.array([0.1, 0.5, 0.9]),
        L=np.full(n, 20.0),
        ref_metrics=ref_metrics,
    )


def _metrics(n=3):
    return {
        'lcp': {'r': np.array([0.9, 1.9, 2.9]),
                'x_L': np.array([0.11, 0.51, 0.91]),
                'L': np.full(n, 21.0)},
        'segments': {'r': np.array([1.5, 6.0, 3.5]),
                     'x_L': np.array([0.12, 0.52, 0.92]),
                     'L': np.full(n, 19.0)},
    }


def test_metrics_for_selects_and_falls_back():
    ctx = _ctx(ref_metrics=_metrics())
    r, x_L, L = ctx.metrics_for('segments')
    np.testing.assert_array_equal(r, [1.5, 6.0, 3.5])
    assert L[0] == 19.0
    # single-strand / undeclared method: canonical arrays
    r, x_L, L = ctx.metrics_for('ecs')
    np.testing.assert_array_equal(r, ctx.r)
    ctx_plain = _ctx(ref_metrics=None)
    r, x_L, L = ctx_plain.metrics_for('segments')
    np.testing.assert_array_equal(r, ctx_plain.r)


def test_filter_masks_ref_metrics():
    ctx = _ctx(ref_metrics=_metrics())
    sub = ctx.filter(np.array([True, False, True]))
    assert len(sub) == 2
    np.testing.assert_array_equal(sub.ref_metrics['segments']['r'], [1.5, 3.5])
    np.testing.assert_array_equal(sub.ref_metrics['lcp']['x_L'], [0.11, 0.91])


def test_adapter_selects_model_declared_metrics():
    """LegacyModelAdapter routes each model to its declared metric set."""
    from openquake.fdha.calc.model_adapter import LegacyModelAdapter
    from openquake.fdha.secondary_surf_rup.visini2025 import (
        Visini2025SecondarySR)
    from openquake.fdha.secondary_surf_rup.petersen2011 import (
        Petersen2011SecondarySR)

    ctx = _ctx(ref_metrics=_metrics())
    visini_adapter = LegacyModelAdapter(
        Visini2025SecondarySR.__new__(Visini2025SecondarySR), {})
    r, _x, _L = visini_adapter._ctx_metrics(ctx)
    np.testing.assert_array_equal(r, [1.5, 6.0, 3.5])       # segments

    petersen_adapter = LegacyModelAdapter(Petersen2011SecondarySR(), {})
    r, _x, _L = petersen_adapter._ctx_metrics(ctx)
    np.testing.assert_array_equal(r, [0.9, 1.9, 2.9])       # lcp
