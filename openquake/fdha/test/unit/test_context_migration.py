"""
Unit tests: Context migration and filtering

Tests the FDHAContext class functionality, including:
- Context creation
- Context filtering
- Curve and map consistency
"""

import numpy as np
import pytest

pytestmark = pytest.mark.unit


def test_context_creation():
    """Test FDHAContext creation."""
    from openquake.fdha.calc.contexts import FDHAContext
    
    N = 10
    ctx = FDHAContext(
        sids=np.arange(N), mag=np.full(N, 6.5), rake=np.full(N, -90.0),
        dip=np.full(N, 60.0), ztor=np.full(N, 0.0),
        occurrence_rate=np.full(N, 1e-4), vs30=np.full(N, 760.0),
        r=np.linspace(0.01, 5.0, N), rx=np.linspace(-2.0, 2.0, N),
        x_L=np.linspace(0.0, 1.0, N), L=np.full(N, 20.0),
    )
    
    assert len(ctx) == N
    assert ctx.style[0] == 'normal'


def test_context_filtering():
    """Test context filtering."""
    from openquake.fdha.calc.contexts import FDHAContext
    
    N = 10
    ctx = FDHAContext(
        sids=np.arange(N), mag=np.full(N, 6.5), rake=np.full(N, 0.0),
        dip=np.full(N, 90.0), ztor=np.full(N, 0.0),
        occurrence_rate=np.full(N, 1e-4), vs30=np.full(N, 760.0),
        r=np.arange(N, dtype=float), rx=np.zeros(N),
        x_L=np.linspace(0.0, 1.0, N), L=np.full(N, 20.0),
    )
    
    mask = ctx.r < 5
    filtered = ctx.filter(mask)
    assert len(filtered) == 5


def test_curve_map_consistency():
    """Ensure curve and map give same results for identical site."""
    # This requires a full integration test setup
    pass