"""
Unit tests: Takao2013 x/L folding behavior

The Takao et al. (2013) regression coefficients in get_prob_D_AD /
get_prob_D_MD are defined against the normalized distance from the *closest*
rupture end (range [0, 0.5]). get_prob therefore must fold the raw along-strike
position x/L in [0, 1] to that range, so the displacement profile is symmetric
about the rupture centre (x/L = 0.5) and peaks there. Without the fold the
gamma mean grows monotonically toward x/L = 1, producing an unphysical
along-strike ramp (see the Chelungpu hazard map north-high/south-low artifact).
"""

import numpy as np
import math
import pytest

from openquake.fdha.primary_surf_displ.takao2013 import Takao2013PrimaryFD

pytestmark = pytest.mark.unit


def _exceed(model, x, norm_disp_type="AD", mag=7.0, d=0.5):
    P = model.get_prob(
        d=np.array([d]), X_L_ratio=np.array([x]), mag=mag, norm_disp_type=norm_disp_type
    )
    assert P.shape == (1, 1)
    return P[0, 0]


@pytest.mark.parametrize("norm_disp_type", ["AD", "MD"])
def test_folding_invariant_under_x_to_1_minus_x(norm_disp_type):
    model = Takao2013PrimaryFD()
    for x in [0.0, 0.1, 0.25, 0.49]:
        p1 = _exceed(model, x, norm_disp_type)
        p2 = _exceed(model, 1.0 - x, norm_disp_type)
        assert math.isfinite(p1) and math.isfinite(p2)
        assert abs(p1 - p2) <= 1e-12


def test_profile_peaks_at_centre_not_at_north_end():
    # Exceedance of a fixed displacement should be largest at the rupture
    # centre and symmetric toward both ends, not monotonically increasing
    # toward x/L = 1.
    model = Takao2013PrimaryFD()
    xs = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
    P = model.get_prob(d=np.array([2.0]), X_L_ratio=xs, mag=7.81, norm_disp_type="AD")[0]
    assert P[2] == pytest.approx(P.max())  # centre is the maximum
    assert P[0] == pytest.approx(P[4], abs=1e-12)  # symmetric ends
    assert P[1] == pytest.approx(P[3], abs=1e-12)
    assert P[2] > P[0]  # centre exceeds the ends


def test_folding_handles_values_outside_unit_interval():
    model = Takao2013PrimaryFD()
    eps = 1e-15
    p0 = _exceed(model, 0.0)
    assert abs(_exceed(model, 1.0 + eps) - p0) <= 1e-12
    assert abs(_exceed(model, -eps) - p0) <= 1e-12
