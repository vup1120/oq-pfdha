"""
Validation of the Nishizaka et al. (2026) far-field models against values
stated in the paper (SRL, doi:10.1785/0220250293).

Two kinds of checks are performed:

1. Paper-anchor tests: the Results section states P2s values read at
   r1 = 10 km for r2 = 0.1 and 10 km (both events, both sides; 250 m
   cells), and 90th-percentile normalized displacements at r1 = 5 and
   10 km for the proximal and non-proximal subsets. The paper reports
   these to one significant figure, so a 15-20% relative tolerance is
   used.

2. Behavioural tests: monotonic decay with r1 and r2, probability bounds,
   vectorization, and input validation.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from openquake.fdha.secondary_surf_rup import Nishizaka2026SecondarySR
from openquake.fdha.secondary_surf_displ import Nishizaka2026SecondaryFD

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# P2s anchors: "At r1 = 10 km, P2s varied with r2 (0.1-10 km) from 4e-2 to
# 3e-3, respectively, on the plain side and from 2e-3 to 4e-5, respectively,
# on the mountainside" (2016 Kumamoto, 250 m cells); "At r1 = 10 km, P2s
# varied from 6e-3 to 2e-3 on the plain side and from 2e-3 to 3e-4 on the
# mountainside for r2 values of 0.1-10 km" (2019 Ridgecrest).
# ---------------------------------------------------------------------------
P2S_PAPER_ANCHORS = [
    # region, side, r1 (km), r2 (km), expected P2s
    ("kumamoto", "plain", 10.0, 0.1, 4e-2),
    ("kumamoto", "plain", 10.0, 10.0, 3e-3),
    ("kumamoto", "mountain", 10.0, 0.1, 2e-3),
    ("kumamoto", "mountain", 10.0, 10.0, 4e-5),
    ("ridgecrest", "plain", 10.0, 0.1, 6e-3),
    ("ridgecrest", "plain", 10.0, 10.0, 2e-3),
    ("ridgecrest", "mountain", 10.0, 0.1, 2e-3),
    ("ridgecrest", "mountain", 10.0, 10.0, 3e-4),
]


@pytest.mark.parametrize("region,side,r1,r2,expected", P2S_PAPER_ANCHORS)
def test_p2s_matches_paper_results_section(region, side, r1, r2, expected):
    model = Nishizaka2026SecondarySR()
    got = model.get_prob(r=r1, r2=r2, region=region, side=side)
    # Paper values are quoted to 1 significant figure
    assert_allclose(got, expected, rtol=0.20)


def test_p2s_decreases_with_r1_and_r2():
    model = Nishizaka2026SecondarySR()
    r1 = np.linspace(0.5, 25.0, 60)
    for region in ("kumamoto", "ridgecrest"):
        for side in ("plain", "mountain"):
            p = model.get_prob(r=r1, r2=1.0, region=region, side=side)
            assert np.all(np.diff(p) < 0), (region, side, "not decreasing in r1")
            r2 = np.linspace(0.1, 15.0, 60)
            p2 = model.get_prob(r=10.0, r2=r2, region=region, side=side)
            assert np.all(np.diff(p2) < 0), (region, side, "not decreasing in r2")
            assert np.all((p >= 0) & (p <= 1))
            assert np.all((p2 >= 0) & (p2 <= 1))


def test_p2s_near_source_approaches_one_for_ridgecrest():
    # "P2s approached 1 regardless of r2 near the earthquake source faults"
    # (Ridgecrest; c3 is ~1e-4 km so the logistic blows up near r1 = 0).
    model = Nishizaka2026SecondarySR()
    for side in ("plain", "mountain"):
        p_near = model.get_prob(r=1e-3, r2=5.0, region="ridgecrest", side=side)
        assert p_near > 0.9, (side, p_near)


def test_p2s_input_validation():
    model = Nishizaka2026SecondarySR()
    with pytest.raises(ValueError):
        model.get_prob(r=1.0, r2=1.0, region="tohoku")
    with pytest.raises(ValueError):
        model.get_prob(r=1.0, r2=1.0, side="valley")
    with pytest.raises(ValueError):
        model.get_prob(r=1.0, r2=1.0, cell_size=100)
    with pytest.raises(ValueError):
        model.get_prob(r=-1.0, r2=1.0)


# ---------------------------------------------------------------------------
# Displacement anchors: "[the 90th percentile of D/EAD for r2 <= 1 km] was
# 0.25 and 0.08 at r1 values of 5 and 10 km, respectively. However, the
# values ... for r2 > 1 km were markedly lower, with a gradual attenuation
# from 0.07 at r1 = 5 km to 0.05 at r1 = 10 km."
# ---------------------------------------------------------------------------
FD_PAPER_ANCHORS = [
    # dataset, r1 (km), expected 90th-percentile D/EAD
    ("proximal", 5.0, 0.25),
    ("proximal", 10.0, 0.08),
    ("nonproximal", 5.0, 0.07),
    ("nonproximal", 10.0, 0.05),
]


@pytest.mark.parametrize("dataset,r1,expected", FD_PAPER_ANCHORS)
def test_fd_90th_percentile_matches_paper(dataset, r1, expected):
    model = Nishizaka2026SecondaryFD()
    ead = 1.9  # 2016 Kumamoto EAD (m), Asano and Iwata (2021), as in paper
    d90 = model.get_percentile_displacement(90, r=r1, ead=ead, dataset=dataset)
    assert_allclose(d90 / ead, expected, rtol=0.15)


def test_fd_dataset_autoselection_by_r2():
    model = Nishizaka2026SecondaryFD()
    prox = model.get_percentile_displacement(90, r=5.0, ead=1.9,
                                             dataset="proximal")
    nonprox = model.get_percentile_displacement(90, r=5.0, ead=1.9,
                                                dataset="nonproximal")
    assert model.get_percentile_displacement(90, r=5.0, ead=1.9, r2=0.5) == prox
    assert model.get_percentile_displacement(90, r=5.0, ead=1.9, r2=2.0) == nonprox
    # proximal sites see larger displacement than non-proximal ones
    assert prox > nonprox


def test_fd_exceedance_consistency_with_percentile():
    # P(D > d90) must equal 0.10 by construction
    model = Nishizaka2026SecondaryFD()
    for dataset in ("all", "proximal", "nonproximal"):
        d90 = model.get_percentile_displacement(90, r=7.0, ead=1.9,
                                                dataset=dataset)
        p = model.get_prob(d=d90, r=7.0, ead=1.9, dataset=dataset)
        assert_allclose(p, 0.10, rtol=1e-6)


def test_fd_exceedance_decays_with_distance_and_threshold():
    model = Nishizaka2026SecondaryFD()
    r = np.array([3.0, 5.0, 10.0, 15.0])
    d = np.array([0.05, 0.2, 0.5])
    p = model.get_prob(d=d, r=r, ead=1.9, dataset="proximal")
    assert p.shape == (4, 3)
    assert np.all(np.diff(p, axis=0) < 0)  # farther -> lower probability
    assert np.all(np.diff(p, axis=1) < 0)  # larger threshold -> lower prob
    assert np.all((p >= 0) & (p <= 1))


def test_fd_input_validation():
    model = Nishizaka2026SecondaryFD()
    with pytest.raises(ValueError):
        model.get_prob(d=0.1, r=5.0, ead=1.9)  # neither dataset nor r2
    with pytest.raises(ValueError):
        model.get_prob(d=0.1, r=5.0, ead=1.9, dataset="riverside")
    with pytest.raises(ValueError):
        model.get_prob(d=0.1, r=5.0, ead=0.0, dataset="all")
    with pytest.raises(ValueError):
        model.get_prob(d=0.1, r=-5.0, ead=1.9, dataset="all")
    with pytest.raises(ValueError):
        model.get_percentile_displacement(0, r=5.0, ead=1.9, dataset="all")
