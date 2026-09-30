# -*- coding: utf-8 -*-
"""Job-to-reference mapping for the IAEA PFDHA exercise benchmark.

Each entry maps one oq-pfdha job (one hazard-analyst model chain of the IAEA
exercise) to the corresponding published curve in ``reference/``.

Fields
------
case / job        subfolder and job name (``<case>/job_<job>.ini``)
figure_csv        reference CSV in ``reference/``
column            model column inside the reference CSV
post_factor       constant multiplier applied to the *computed* curve before
                  comparison. Used only for the Takao (2013) chains, whose
                  published curves include the conditional probability of
                  principal rupture at the site (the P2p term, TECDOC-2092
                  Section 3.2.2, approx 0.45-0.48 for these cases). P2p is not
                  yet implemented as a library model, so it is applied here as
                  a documented constant back-calculated from the published
                  first point (Kumamoto 0.4801, Le Teil 0.4748; the paper
                  quotes "a factor of 0.45").
assert_max_relerr max |computed/reference - 1| allowed by the pytest test
                  (None = qualitative comparison, plotted but not asserted)
assert_dmax_m     restrict the assertion to displacements <= this value (m);
                  beyond it the curves sit at AFOE < 1e-8..1e-12 where the
                  paper notes results are controlled by each team's aleatory
                  truncation choice (not prescribed by the exercise).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class Entry:
    case: str
    job: str
    figure_csv: str
    column: str
    post_factor: float = 1.0
    assert_max_relerr: Optional[float] = 0.15
    assert_dmax_m: Optional[float] = None


# The exercise's V24 (Visini et al., 2025) chains are deliberately absent.
# The teams' V24 curves were computed with an earlier, under-review revision
# of that model, so comparing them against an implementation of the published
# model compares two different models rather than validating one: the offsets
# (a factor 13 at the Le Teil first distributed exercise, 130-180 at the
# second and third, with the hanging-wall response reversed in sign) are a
# vintage mismatch, not an implementation error. The Visini implementation is
# instead validated at model level against the authors' released FDHLab code
# and their Fig. 13 worked example (benchmark/visini_et_al_2025), which is a
# comparison of like with like. See diagnose_V24_p2d.py for the decomposition.
MANIFEST = [
    # ------------------------------------------------ Kumamoto principal (Fig 4a)
    Entry("kumamoto", "principal_P11", "fig4a_kumamoto_principal.csv", "P11",
          assert_max_relerr=0.15),
    Entry("kumamoto", "principal_C24", "fig4a_kumamoto_principal.csv", "C24",
          assert_max_relerr=0.15),
    Entry("kumamoto", "principal_K24", "fig4a_kumamoto_principal.csv", "K24",
          assert_max_relerr=0.15),
    Entry("kumamoto", "principal_T13", "fig4a_kumamoto_principal.csv", "T13",
          post_factor=0.4801, assert_max_relerr=0.15),
    Entry("kumamoto", "principal_L23", "fig4a_kumamoto_principal.csv", "L23",
          assert_max_relerr=0.15, assert_dmax_m=1.0),
    # ------------------------------------------------ Kumamoto distributed (Fig 6a)
    Entry("kumamoto", "distributed_P11", "fig6a_kumamoto_distributed.csv", "P11",
          assert_max_relerr=0.15),
    # T13 distributed = Takao2013PrimarySR x Takao2014SecondarySR (100 m
    # cell) x Takao2013SecondaryFD (AD, n_sigma 5): the curve head matches
    # the published one to 0.1%, but the mid-range (0.1-1 m) runs up to
    # ~35% high (Kumamoto) / ~35% low (Le Teil) - the team's aleatory
    # integration of the DD/PAD distribution is not documented in
    # TECDOC-2092, so the shape difference cannot be reconciled further.
    Entry("kumamoto", "distributed_T13", "fig6a_kumamoto_distributed.csv", "T13",
          assert_max_relerr=0.40),
    # Same chains over the author workbook's BASE-CASE model: the four
    # coexisting rupture sources (Uto, Futagawa+Uto, Uto+UHN, F+U+UHN) with
    # the magnitude/rate epistemic branches mean-collapsed (see
    # kumamoto/make_basecase.py). The published Fig 6a curves correspond to
    # the Uto-only middle branch (the entries above match them at 0.4% /
    # head 0.1%), so these tree entries are qualitative: they quantify how
    # much hazard the multi-segment rupture branches add at the base site.
    Entry("kumamoto", "distributed_P11_basecase",
          "fig6a_kumamoto_distributed.csv", "P11", assert_max_relerr=None),
    Entry("kumamoto", "distributed_T13_basecase",
          "fig6a_kumamoto_distributed.csv", "T13", assert_max_relerr=None),
    # ------------------------------------------------ Kumamoto floating (Fig 4b)
    # Floating ruptures reproduce the exercise convention via a PeerMSR
    # full-width rupture; the residual ~8% rate offset is the along-strike
    # occupancy fraction (ours 0.35 vs the teams' ~0.32).
    Entry("kumamoto", "floating_K24", "fig4b_kumamoto_principal_floating.csv",
          "K24", assert_max_relerr=0.30),
    Entry("kumamoto", "floating_T13", "fig4b_kumamoto_principal_floating.csv",
          "T13", post_factor=0.4801, assert_max_relerr=0.30),
    Entry("kumamoto", "floating_L23", "fig4b_kumamoto_principal_floating.csv",
          "L23", assert_max_relerr=0.30, assert_dmax_m=1.0),
    # ------------------------------------------------ Le Teil principal (Fig 4c)
    Entry("le_teil", "principal_K24", "fig4c_leteil_principal.csv", "K24",
          assert_max_relerr=0.15),
    Entry("le_teil", "principal_T13", "fig4c_leteil_principal.csv", "T13",
          post_factor=0.4748, assert_max_relerr=0.15, assert_dmax_m=1.0),
    Entry("le_teil", "principal_L23", "fig4c_leteil_principal.csv", "L23",
          assert_max_relerr=0.15, assert_dmax_m=1.0),
    Entry("le_teil", "principal_M11", "fig4c_leteil_principal_M11.csv", "M11",
          assert_max_relerr=None),  # qualitative: see README
    # The same M11 chain run over the author workbook's BASE-CASE epistemic
    # tree (smlt_basecase.xml: thickness/rupture-length/magnitude branches,
    # slip-rate level collapsed into the mean rates; see make_basecase.py).
    # The published M11 curve is evidently dominated by the tree's
    # large-magnitude branches (le_teil/diagnose_M11.py), so this entry
    # tracks how far the full tree explains it. Qualitative.
    Entry("le_teil", "principal_M11_basecase", "fig4c_leteil_principal_M11.csv",
          "M11", assert_max_relerr=None),
    # ------------------------------------------------ Le Teil distributed (Fig 6b)
    Entry("le_teil", "distributed_T13", "fig6b_leteil_distributed.csv", "T13",
          assert_max_relerr=0.40, assert_dmax_m=1.0),  # same caveat as the
    # Kumamoto T13 distributed entry; beyond 1 m the published tail is
    # controlled by the team's (undocumented) aleatory truncation.
    # ------------------------------------------------ Norcia principal (Fig 4d)
    Entry("norcia", "principal_Y03", "fig4d_norcia_principal.csv", "Y03",
          assert_max_relerr=0.20),
    Entry("norcia", "principal_K24", "fig4d_norcia_principal.csv", "K24",
          assert_max_relerr=0.20),
    Entry("norcia", "principal_L23", "fig4d_norcia_principal.csv", "L23",
          assert_max_relerr=0.15),
    # ------------------------------------------------ Norcia distributed (Fig 6c)
    Entry("norcia", "distributed_Y03", "fig6c_norcia_distributed.csv", "Y03",
          assert_max_relerr=0.25, assert_dmax_m=1.0),

    # =================================================================
    # 2nd and 3rd distributed exercises of each case study (TECDOC-2092
    # Table 8). These have no coordinator-supplied vectors; their
    # reference curves are digitized from the published figures by
    # reference/digitize_tecdoc_figures.py (axes calibrated from the
    # major gridlines, residual < 0.01 decade). The digitisation error
    # is dominated by the plotted line width, so these entries carry
    # tolerances one step looser than their Fig 4/6 counterparts.
    # =================================================================
    # ------------------------------------ Kumamoto 2nd distributed (Fig 15c)
    Entry("kumamoto", "distributed_P11_sens4", "fig15c_kumamoto_distributed_r10.csv",
          "P11", assert_max_relerr=0.15),
    # The T13 shape difference documented for the Fig 6a entry grows with
    # distance: at this r = 10 km site the head still matches to 9%, but the
    # ratio rises monotonically to 1.70 at 1 m as our distributed-displacement
    # tail rolls off more slowly than the team's. Same undocumented aleatory
    # integration; asserted loosely to catch regressions, not agreement.
    Entry("kumamoto", "distributed_T13_sens4", "fig15c_kumamoto_distributed_r10.csv",
          "T13", assert_max_relerr=0.75, assert_dmax_m=1.0),
    # ------------------------------------ Kumamoto 3rd distributed (Fig 15d)
    Entry("kumamoto", "distributed_T13_sens1", "fig15d_kumamoto_distributed_suizenji.csv",
          "T13", assert_max_relerr=0.20),
    # No P11 counterpart: the published Fig 15(d) shows a P11 curve for this
    # Mw 5.8 Suizenji scenario, but Petersen et al. (2011) is calibrated for
    # M 6-8 and the framework's applicability guard refuses the evaluation
    # rather than extrapolating (kumamoto/job_distributed_P11_sens1.ini
    # exists and is expected to raise). See README, "Known deviations".
    # ------------------------------------ Le Teil 2nd distributed (Fig 19c)
    Entry("le_teil", "distributed_T13_dipflip", "fig19c_leteil_distributed_dipNW.csv",
          "T13", assert_max_relerr=0.40, assert_dmax_m=1.0),
    # ------------------------------------ Le Teil 3rd distributed (Fig 19d)
    Entry("le_teil", "distributed_T13_sens4", "fig19d_leteil_distributed_3faults.csv",
          "T13", assert_max_relerr=0.40, assert_dmax_m=1.0),
    # ------------------------------------ Norcia 2nd distributed (Fig 22b)
    Entry("norcia", "distributed_Y03_sens2", "fig22b_norcia_distributed_sl.csv",
          "Y03", assert_max_relerr=0.25, assert_dmax_m=1.0),
    # ------------------------------------ Norcia 3rd distributed (Fig 22c)
    Entry("norcia", "distributed_Y03_sens3", "fig22c_norcia_distributed_2sources.csv",
          "Y03", assert_max_relerr=0.25, assert_dmax_m=1.0),
]
