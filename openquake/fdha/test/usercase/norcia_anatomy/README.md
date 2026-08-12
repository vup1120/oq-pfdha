# Norcia epistemic anatomy (paper application, Section 7 rework)

Builds on the `test/benchmark/norcia_case3_iaea` case (same two-fault
source model, same three IAEA sites) but replaces the two-chain
verification tree with the **full normal-faulting model library**,
evaluated under the rupture-location weight that the IAEA mapping-accuracy
class implies.

## Sites (IAEA TECDOC-2092, Table 13)

| id | name | position | role |
|----|------|----------|------|
| PF | southern Mount Vettore trace | 42.767 N, 13.278 E, `l/L` = 0.05 | on-trace principal site |
| MS | Mount Serra antithetic fault | 42.749 N, 13.188 E | distributed site at the San Benedetto tunnel; `r` = 7.6 km from the MVF hanging wall, 3.0 km from the NF footwall |
| SL | San Lorenzo antithetic fault | 42.853 N, 13.212 E, `l/L` = 0.4 | distributed site at the closest approach to the MVF trace, `r` = 2.4 km |

Both MS and SL are pre-existing antithetic structures of the Mount
Vettore fault system that were reactivated in 2016; the exercise treats
them as *sites*, not as sources.

## Source tree

`source_model_logic_tree.xml` is a **single `sourceModel` branch at
weight 1.0** - there is deliberately no source-model uncertainty, so
every result here is FDHA-model epistemic uncertainty on the fixed IAEA
source model. Adding parametric source branches would be constrained
anyway: both faults are `characteristicFaultSource` with an
`EvenlyDiscretizedMFD`, and of the OpenQuake uncertainty types only
`incrementalMFDAbsolute` applies to them - `simpleFaultDipAbsolute`,
`setLowerSeismDepthAbsolute`, `setMSRAbsolute`, `maxMagGRRelative` and
`recomputeMmax` are all rejected by that source class. Geometry
branches would require alternative sourceModel XML files.

## FDHA logic trees

- `fdha_logic_tree_anatomy.xml` - **60 end branches**:
  P_sr {Pizza2023 w=0.4, Youngs2003 at its three Pezzopane & Dawson
  normal-faulting data sets w=0.2 each} x P_fd {Y03-AD, Y03-MD,
  Lavrentiadis2023-principal} x S_sr {Y03, Visini2025} with the
  published S_fd pairings {Y03 85th/95th | Visini
  WC1994/TB17/Leonard}, both distributed models on a common
  500 x 500 m cell (see below).
  **Wells & Coppersmith (1993) is excluded** from the P_sr node: 276
  worldwide earthquakes of all slip types is a weaker analogue for a
  normal-faulting case than any of the three normal-faulting sets, the
  same scope test applied to Ferrario & Livio. The three Youngs
  branches are nested geographic scopes of one western-US compilation,
  not independent studies - equal weights express how wide an analogue
  is accepted.
  **Ferrario & Livio (2021) is excluded**: calibrated on
  volcano-tectonic normal faulting, out of scope for the purely
  tectonic Norcia sources.
  Every Youngs branch pins `version` explicitly: a bare
  `[Youngs2003PrimarySR]` takes its coefficients from the rupture rake,
  not from the class default, so an unpinned branch is not
  self-describing.
  `fdhaCalcRSigma` is fixed at **0.02689 km on every branch** - the
  Petersen et al. (2011) Table 2 two-sided sigma for the "Accurately
  located" class that TECDOC-2092 Table 13 assigns to these traces - so
  W_p is the Petersen Gaussian throughout and `r_threshold_km` in the
  INI is inert (design doc D2).
  The primary-SR level uses the proxy-pair idiom (each model duplicated
  `_YDEF`/`_LAV`) so the two primary-FD displacement definitions
  (Y03 principal vs LA23 sum-of-principal) stay in separate,
  definition-consistent branch sets (FDLT-014).
  Kuehn et al. (2024) is excluded: aggregate definition, cannot share a
  chain with distributed models (FDLT-013). Mammarella et al. (2024) is
  excluded: requires case-specific hypocentral-depth priors not
  published for this region.
- `fdha_logic_tree_profile.xml` - IAEA-baseline Y03 chain with the same
  Gaussian W_p, used for the cross-trace profile.
- `fdha_logic_tree_cellsize.xml` - the Visini distributed chain only,
  at 100 m and 500 m site cells, primary levels identical to the
  anatomy tree (72 branches). Its 500 m arm reproduces the anatomy-tree
  Visini curve to 1e-9, which `make_figures.py` asserts.

**The distributed cell convention.** Distributed-rupture probability is
a probability *per cell*, and the cell belongs to the model:

| model | cell | exposed as a parameter? |
|-------|------|-------------------------|
| Youngs et al. (2003) | 500 x 500 m | no - baked into the coefficients |
| Visini et al. (2025) | 10/20/50/100/200/500 m | yes (`pixel_size`) |

Youngs is a fixed 500 m model, so the tree pins Visini to 500 m.
Consequence worth noting: the IAEA exercise specifies a 100 x 100 m
site, and no published tectonic normal-faulting distributed model
except Visini can answer at that scale.

The cell effect is **conditional on the site classification**: running
Visini at 100 m instead of 500 m divides its rate of exceeding the 2016
tunnel offset by 8.2 with no pre-existing fault, but by only 1.04 with
one. Only combinations A and B carry the cell-dependent along-strike
Monte Carlo term; combination C sets that probability to 1 by design
(Visini et al. 2025 p. 12), so once C dominates the cell washes out.
See `discussion_notes.md` A4.

**What the models return.** Both give a full exceedance probability
P(D > d), not a point value. Youngs draws normalised displacement D/MD
from a gamma (shape 2.5) anchored so its **85th or 95th percentile**
matches the distance regression, integrated over a WC1994 lognormal
maximum displacement - so "85th/95th" are two published *fits*, not
output percentiles, and both are carried as branches. Visini predicts
the **median** throw with a lognormal residual (sigma_lnY = 1.027,
truncated at +-3 sigma).

**Visini combinations.** `[calculation].case` selects the Visini et al.
(2025) site decision tree case (`calc/decision_tree.py`), and is the only
consumer of that key:

| case | meaning | combinations |
|------|---------|--------------|
| `case1` | a Rank 1.5 structure beneath the site cannot be excluded | A + B + C |
| `case2` | none beneath the site, but one within 1 km | A + B |
| `case3` | excluded beneath the site and within 1 km | A only |

This study declares **`case1`**: MS and SL are pre-existing antithetic
faults of the Mount Vettore system and both ruptured at the surface in
2016, so the sites sit ON Rank 1.5 structures. That is a declaration
from evidence; `rank1p5_traces_file` is only needed when the case has to
be *derived* from trace geometry, which we do not have for the Mount
Serra and San Lorenzo faults.

> Watch the name collision: these INIs were seeded from the IAEA
> benchmark, where `case = case3` means the exercise's own **scenario**
> "Norcia Case 3". Consumed here it silently meant the decision-tree
> case3 and restricted the run to combination A.

## Running

From this directory (~1 min curves, ~15 min map):

```
python -c "from pathlib import Path; from openquake.fdha.logic_tree.driver import FdhaLogicTree; \
  FdhaLogicTree.from_ini(Path('job_anatomy_curve.ini')).run(outdir=Path('out_curve'))"
# likewise job_profile_transect.ini -> out_transect, job_anatomy_map.ini -> out_map
python make_figures.py
```

- `job_anatomy_curve.ini` - PF / MS / SL hazard curves, 60 branches.
- `job_profile_transect.ini` - 101 sites every 20 m on a 2 km profile
  perpendicular to the southern MVF trace through PF
  (`transect_sites.txt`, offsets in `transect_offsets_km.txt`;
  negative offsets = hanging-wall side).
- `job_anatomy_map.ini` - regional map, 0.01 deg grid, 1e5 a return
  period, 60 branches (~7 min). Feeds panel (a) of `fig_ms_tunnel`.
- `job_cellsize_curve.ini` - Visini chain at 100 m vs 500 m site cell,
  three sites, 72 branches (feeds the dashed curve in
  `fig_ms_tunnel` panel (b)).

## Figures (`figures/`, PDF + PNG)

The paper carries **three** figures, and only these three are written
to disk (`PAPER_FIGURES` in `make_figures.py`):

| figure | content |
|---|---|
| `fig_logic_tree` | the five-level tree, 4 x 3 x (2+3) x 1 = 60 end branches, with weights. Node membership, connectivity and weights are all recovered from the branches that were actually run, so the figure cannot drift from the calculation |
| `fig_ms_tunnel` | mean distributed hazard curves at MS and SL against the 2016 observations (0.20 m tunnel offset; ~1 m road offset). Youngs 85th/95th at 500 m, and Visini at 500/100 m cells with and without a pre-existing fault beneath the site |
| `fig_principal_fractiles` | (a) mean hazard and (b) 84th epistemic fractile of the PRINCIPAL displacement at 1e5 a, on-fault sites only; (c) all 60 branch curves at PF with the mean and the 16th-84th band |

Three further analyses still run and print their numbers to the console
- they are quoted in `discussion_notes.md` and in the section text -
but are no longer rendered: the cross-trace profile (section C), the
three-site branch fans (section B), and the Abrahamson & Liou design
basis (section F). To render them again, add their names to
`PAPER_FIGURES`.

Also dropped earlier: `fig_maps` (its 84th-percentile export carries no
principal column) and `fig_strip` (the fault-parallel strip product was
trialled and abandoned).

Headline numbers are collected in `discussion_notes.md`; the draft
section text is `section7_draft.tex`.
