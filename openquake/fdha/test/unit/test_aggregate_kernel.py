# -*- coding: utf-8 -*-
"""
Aggregate single-bucket kernel path (C4, design D8).

When the principal-slot FD model declares DISPLACEMENT_DEFINITION ==
"aggregate", the whole contribution is

    lambda_total = rate * P_sr * P_fd_aggregate * W_p(r)

in the PRINCIPAL bucket, the distributed bucket is exactly zero and the
secondary models are never evaluated -- on BOTH W_p paths (sigma = 0 boxcar
and sigma > 0 Gaussian). Non-aggregate chains must be bit-for-bit unchanged.
"""
import numpy as np
import pytest

from openquake.fdha.calc.contexts import FDHAContext
from openquake.fdha.calc.hazard import _compute_rupture_contribution
from openquake.fdha.calc.location_weight import location_weight

pytestmark = pytest.mark.unit

RATE = 1e-4
P_SR = 0.8
R_KM = np.array([0.0, 0.05, 0.39, 1.2])  # on-trace, inside h, off, far
D0 = np.array([0.01, 0.1, 1.0])
H_KM = 0.1


class _AggregateModel:
    DISPLACEMENT_DEFINITION = "aggregate"
    DISPLACEMENT_COMPONENT = "net"


class _PrincipalModel:
    DISPLACEMENT_DEFINITION = "principal"
    DISPLACEMENT_COMPONENT = "net"


class _SumOfPrincipalModel:
    """Static sum-of-principal contract (e.g. Chiou2025PrimaryFD,
    Lavrentiadis2023PrimaryFD_principal): NOT aggregate, so the two-bucket
    path applies and secondary models stay legitimate. The contract is the
    plain class attribute -- the class choice IS the definition."""
    DISPLACEMENT_DEFINITION = "sum-of-principal"
    DISPLACEMENT_COMPONENT = "net"


class _PrimarySRAdapter:
    model = None
    model_params: dict = {}

    def compute_primary_sr(self, ctx, red_cfg):
        return np.full(len(ctx), P_SR)


class _PrimaryFDAdapter:
    def __init__(self, model, model_params=None):
        self.model = model
        self.model_params = model_params or {}
        # deterministic, distinguishable P_fd matrix
        self.p_fd = None

    def compute_primary_fd(self, ctx, displacements, red_cfg):
        n, d = len(ctx), len(displacements)
        self.p_fd = (np.arange(n * d, dtype=np.float64).reshape(n, d) + 1.0) \
            / (n * d + 1.0)
        return self.p_fd


class _RecordingSecondaryAdapter:
    """Nonzero secondary contribution + call recording: on the aggregate
    path it must never be invoked."""

    def __init__(self):
        self.calls = 0
        self.model = None
        self.model_params: dict = {}

    def compute_secondary_sr(self, ctx, red_cfg):
        self.calls += 1
        return np.full(len(ctx), 0.5)

    def compute_secondary_fd(self, ctx, displacements, red_cfg):
        self.calls += 1
        return np.full((len(ctx), len(displacements)), 0.25)


def _ctx():
    n = len(R_KM)
    return FDHAContext(
        sids=np.arange(n),
        mag=np.full(n, 7.0),
        rake=np.zeros(n),
        dip=np.full(n, 90.0),
        ztor=np.zeros(n),
        occurrence_rate=np.full(n, RATE),
        vs30=np.full(n, 760.0),
        r=R_KM.copy(),
        rx=R_KM.copy(),
        x_L=np.full(n, 0.5),
        L=np.full(n, 40.0),
    )


def _run(primary_fd_adapter, r_sigma_km, with_secondary=True):
    sec_sr = _RecordingSecondaryAdapter()
    sec_fd = _RecordingSecondaryAdapter()
    adapters = {
        "primary_sr": _PrimarySRAdapter(),
        "primary_fd": primary_fd_adapter,
    }
    if with_secondary:
        adapters["secondary_sr"] = sec_sr
        adapters["secondary_fd"] = sec_fd
    principal, distributed = _compute_rupture_contribution(
        ctx=_ctx(),
        adapters=adapters,
        target_displacements=D0,
        p_sr_red_cfg={},
        s_sr_red_cfg={},
        r_threshold_km=H_KM,
        use_visini=False,
        visini_calc=None,
        calculator=None,
        r_sigma_km=r_sigma_km,
    )
    return principal, distributed, sec_sr, sec_fd, primary_fd_adapter


def _expected_principal(p_fd, r_sigma_km):
    w_p = location_weight(R_KM, r_threshold_km=H_KM, r_sigma_km=r_sigma_km)
    p_sr = np.full(len(R_KM), P_SR)
    return RATE * p_sr[:, np.newaxis] * p_fd * w_p[:, np.newaxis]


@pytest.mark.parametrize("r_sigma_km", [0.0, 0.3],
                         ids=["sigma0_boxcar", "sigma_gaussian"])
def test_aggregate_single_bucket_exact(r_sigma_km):
    principal, distributed, sec_sr, sec_fd, pfd = _run(
        _PrimaryFDAdapter(_AggregateModel()), r_sigma_km)

    # lambda = rate * P_sr * P_fd_aggregate * W_p(r), exactly.
    np.testing.assert_array_equal(
        principal, _expected_principal(pfd.p_fd, r_sigma_km))
    # Distributed bucket is exactly zero and secondary models untouched.
    assert np.all(distributed == 0.0)
    assert sec_sr.calls == 0
    assert sec_fd.calls == 0


@pytest.mark.parametrize("model_cls", [_PrincipalModel, _SumOfPrincipalModel],
                         ids=["principal", "sum_of_principal"])
@pytest.mark.parametrize("r_sigma_km", [0.0, 0.3],
                         ids=["sigma0_boxcar", "sigma_gaussian"])
def test_non_aggregate_path_unchanged(r_sigma_km, model_cls):
    """Control: principal- and sum-of-principal-definition models keep the
    legacy two-bucket combination (complementary at sigma=0, additive at
    sigma>0); only 'aggregate' single-buckets."""
    principal, distributed, sec_sr, sec_fd, pfd = _run(
        _PrimaryFDAdapter(model_cls()), r_sigma_km)

    np.testing.assert_array_equal(
        principal, _expected_principal(pfd.p_fd, r_sigma_km))

    w_p = location_weight(R_KM, r_threshold_km=H_KM, r_sigma_km=r_sigma_km)
    g = (1.0 - w_p) if r_sigma_km == 0.0 else np.ones_like(w_p)
    p_sr = np.full(len(R_KM), P_SR)
    p_dist = p_sr * 0.0 + 0.5  # secondary SR
    expected_dist = (RATE * p_sr[:, np.newaxis]
                     * (p_dist[:, np.newaxis] * np.full((len(R_KM), len(D0)),
                                                        0.25))
                     * g[:, np.newaxis])
    np.testing.assert_allclose(distributed, expected_dist, rtol=0, atol=0)
    assert sec_sr.calls == 1
    assert sec_fd.calls == 1
    # sigma=0: distributed masked inside the boxcar but present outside.
    if r_sigma_km == 0.0:
        assert np.all(distributed[R_KM <= H_KM] == 0.0)
        assert np.all(distributed[R_KM > H_KM] > 0.0)
    else:
        assert np.all(distributed > 0.0)


def test_kernel_routing_ignores_model_params():
    """Contracts are STATIC: the kernel routes on the class attribute alone,
    so model parameters (whatever they are) cannot re-route an aggregate
    class to the two-bucket path or vice versa."""
    _, distributed, sec_sr, _, _ = _run(
        _PrimaryFDAdapter(_AggregateModel(),
                          {"some_param": "disp_prnc_prime"}),
        r_sigma_km=0.0)
    assert np.all(distributed == 0.0)
    assert sec_sr.calls == 0

    _, distributed, sec_sr, _, _ = _run(
        _PrimaryFDAdapter(_SumOfPrincipalModel(),
                          {"some_param": "whatever"}),
        r_sigma_km=0.0)
    assert np.any(distributed > 0.0)
    assert sec_sr.calls == 1


def test_aggregate_total_equals_principal_bucket():
    """Output convention: the aggregate contribution flows through the
    principal bucket, so total == principal exactly."""
    principal, distributed, *_ = _run(
        _PrimaryFDAdapter(_AggregateModel()), r_sigma_km=0.3)
    np.testing.assert_array_equal(principal + distributed, principal)


# ---------------------------------------------------------------------------
# Calculator configuration guard (direct [models.*] jobs, no logic tree):
# aggregate primary FD + configured secondary models must fail LOUDLY at
# model initialisation, mirroring validator FDLT-013.
# ---------------------------------------------------------------------------
def _calc_with_models(models_cfg):
    from openquake.fdha.calc.calculators import BaseFaultRuptureCalculator
    calc = BaseFaultRuptureCalculator.__new__(BaseFaultRuptureCalculator)
    calc.config = {"models": models_cfg}
    calc._initialize_models()
    return calc


def test_calculator_rejects_aggregate_plus_secondary_config():
    with pytest.raises(ValueError, match="AGGREGATE"):
        _calc_with_models({
            "primary_surf_displ": {"type": "Kuehn2024PrimaryFD"},
            "secondary_surf_rup": {"type": "FixedSecondarySR",
                                   "parameters": {"value": 0.0}},
            "secondary_surf_displ": {"type": "Youngs2003SecondaryFD"},
        })


def test_calculator_rejects_aggregate_plus_secondary_fd_only():
    with pytest.raises(ValueError, match="double counts"):
        _calc_with_models({
            "primary_surf_displ": {"type": "Lavrentiadis2023PrimaryFD"},
            "secondary_surf_displ": {"type": "Youngs2003SecondaryFD"},
        })


def test_calculator_allows_aggregate_without_secondary():
    calc = _calc_with_models({
        "primary_surf_displ": {"type": "Kuehn2024PrimaryFD"},
    })
    assert calc.secondary_surf_displ_model is None


def test_calculator_allows_lavrentiadis_principal_class_with_secondary():
    """The class choice IS the definition: the sum-of-principal
    Lavrentiadis2023PrimaryFD_principal variant class may carry secondary
    models (IAEA L23 setup, converted to the class-based scheme)."""
    calc = _calc_with_models({
        "primary_surf_displ": {
            "type": "Lavrentiadis2023PrimaryFD_principal",
            "parameters": {"include_zero_slip": True},
        },
        "secondary_surf_rup": {"type": "FixedSecondarySR",
                               "parameters": {"value": 0.0}},
        "secondary_surf_displ": {"type": "Youngs2003SecondaryFD"},
    })
    assert calc.secondary_surf_displ_model is not None


def test_calculator_allows_principal_primary_with_secondary():
    calc = _calc_with_models({
        "primary_surf_displ": {"type": "Youngs2003PrimaryFD"},
        "secondary_surf_displ": {"type": "Youngs2003SecondaryFD"},
    })
    assert calc.secondary_surf_displ_model is not None
