# Headline numbers for the Section 7 rework

From the 60-branch anatomy tree: IAEA source model, Petersen Gaussian
W_p (sigma = 26.89 m) on every branch, both distributed models on a
500 x 500 m cell, **tectonic normal faulting only**. Outputs in
`out_curve`, `out_transect`, `out_map`, with the cell-size sensitivity
in `out_cellsize` (and its `out_cellsize_case3` twin). "RP" = 1/annual
rate. The observed 2016 San Benedetto tunnel offset is 0.20 m vertical
+ 0.13 m left-lateral (net 0.24 m; Galli et al. 2019 WTC proceedings /
Galli et al. 2020); the San Lorenzo road offset is of order 1 m.

All runs use `case = case1` - a pre-existing fault beneath the site
cannot be precluded - which is the evidence-based classification for
MS and SL (section A). Numbers in earlier drafts of this file were
produced under an inherited `case = case3` and understated the
distributed hazard by x137; see A1.

## 0. Scope and conventions (these condition every number below)

**The source tree carries no uncertainty.** `source_model_logic_tree.xml`
is a single `sourceModel` branch at weight 1.0, so everything below is
FDHA-model epistemic uncertainty on the fixed IAEA source model. Note
that the source model would resist parametric branching anyway: both
faults are `characteristicFaultSource` with an explicit
`incrementalMFD`, and of the OpenQuake uncertainty types only
`incrementalMFDAbsolute` applies - `simpleFaultDipAbsolute`,
`setLowerSeismDepthAbsolute`, `setMSRAbsolute`, `maxMagGRRelative` and
`recomputeMmax` are all rejected by that source class (tested).
Geometry branches would need alternative sourceModel XML files.

The two faults share a source CLASS but not an MFD SHAPE, and the
difference drives several results below (verified from the rate lists in
`source_model_norcia_case3.xml`):

| | MVFS | NFS |
|---|---|---|
| Mw range | 6.4-7.0 | 5.5-6.8 |
| shape | symmetric bell, peak Mw 6.7 | successive ratio constant at 0.7943, i.e. **exponential with b = 1.000** |
| total rate | 4.03e-04 /a | 5.41e-03 /a (**13.4x**) |
| rate-weighted mean Mw | 6.70 | 5.83 |

`characteristicFaultSource` is an OpenQuake rupture TYPOLOGY (the whole
fault surface ruptures); it does not imply a characteristic magnitude
distribution. NFS carries 13.4x the rate at a mean magnitude 0.87 units
lower, which is why it dominates the distributed hazard at MS (2.91 km
away) while contributing little principal displacement.

**The P_sr node uses the normal-faulting data sets.** Youngs et al.
(2003) Equation 4 is one logistic curve fitted to four different data
sets; the `version` parameter selects which. The node now carries the
three Pezzopane & Dawson (1996) normal-faulting sets plus the Italian
Pizza et al. (2023) model:

| branch | weight | data set |
|--------|--------|----------|
| Pizza et al. (2023), normal | 0.4 | Italian normal faults - the only regionally appropriate regression |
| Youngs2003, ExtensionalCordillera | 0.2 | 105 earthquakes |
| Youngs2003, NorthernBasinAndRange | 0.2 | 47 earthquakes |
| Youngs2003, GreatBasin | 0.2 | 32 earthquakes |

Wells & Coppersmith (1993) is dropped: 276 worldwide earthquakes of ALL
slip types is a weaker analogue for a normal-faulting case than any of
the three above, the same scope test that excludes Ferrario & Livio.
(It stays reachable as `version = WC93`, which is bit-identical to
`WC1993PrimarySR`.)

CAVEAT worth one sentence in the paper: the three Youngs branches are
progressively wider geographic scopes of ONE western-US extensional
compilation (Great Basin subset of northern Basin & Range subset of
extensional cordillera), not three independent studies. Equal weights
express "how wide an analogue do I accept", not three independent lines
of evidence.

**Ferrario & Livio (2021) is excluded.** It is calibrated on
volcano-tectonic normal faulting and is out of scope for the purely
tectonic Mount Vettore / Norcia sources. This matters: in the earlier
4-model version it was the single outlier that carried most of the
apparent library spread.

**The distributed cell convention.** Distributed-rupture probability is
a probability *per cell*, and the cell belongs to the model:

| model | cell | exposed as a parameter? |
|-------|------|-------------------------|
| Youngs et al. (2003) | 500 x 500 m | no - baked into the coefficients |
| Visini et al. (2025) | 10/20/50/100/200/500 m | yes (`pixel_size`) |

Youngs is a fixed 500 m model, so Visini is pinned to 500 m to keep the
two comparable. Corollary: the IAEA exercise specifies a 100 x 100 m
site, and no published tectonic normal-faulting distributed model
except Visini can answer at that scale.

**What the models actually return.** Both return a full exceedance
probability P(D > d), not a point value:

- *Youngs et al. (2003)*: normalised displacement D/MD follows a gamma
  distribution with shape a = 2.5, whose scale is set so that its
  **85th (or 95th) percentile equals the distance-dependent regression
  value** x(r) = 0.35 exp(-0.091 r) on the hanging wall,
  0.16 exp(-0.137 r) on the footwall. (Verified: the code's scaling
  constants 4.058 and 5.535 are exactly `gamma.ppf(0.85, 2.5)` and
  `gamma.ppf(0.95, 2.5)`.) That gamma is then integrated over a
  Wells & Coppersmith (1994) lognormal distribution of the principal
  maximum displacement MD (+-3 sigma, sigma_log10 = 0.38).
  **The 85th/95th are two published fits to different percentiles of
  the D/MD data scatter - an epistemic choice between regressions, NOT
  an output percentile.** Both are carried as separate branches.
- *Visini et al. (2025)*: the regression predicts the **median** throw
  Y (a fit on ln Y), with a lognormal residual sigma_lnY = 1.0271
  truncated at +-3 sigma; P(Y > d) integrates that lognormal.

So neither model outputs "a median" or "a mean" to the hazard integral
- each supplies a distribution. Their central tendencies are anchored
differently (Youngs at a high percentile of normalised displacement,
Visini at the median of absolute throw), which is worth one sentence in
the paper because it is a common source of confusion.

## A. The two antithetic-fault sites (fig_ms_tunnel)

Rate of exceeding what was actually measured in 2016, per distributed
chain (weighted mean over the primary nodes that gate it). MS: the
San Benedetto tunnel, 0.20 m vertical. SL: a road offset of order 1 m.

| chain | MS rate (1/a) | MS RP (a) | SL rate (1/a) | SL RP (a) |
|-------|--------------|-----------|--------------|-----------|
| Youngs (2003), 85th, 500 m | 5.27e-07 | 1,898,000 | 5.23e-07 | 1,914,000 |
| Youngs (2003), 95th, 500 m | 3.40e-07 | 2,938,000 | 2.88e-07 | 3,471,000 |
| **Visini (2025), 500 m, pre-existing fault** | **4.88e-05** | **20,500** | **7.81e-06** | **128,000** |
| Visini (2025), 100 m, pre-existing fault | 4.70e-05 | 21,300 | 7.49e-06 | 134,000 |
| Visini (2025), 500 m, no pre-existing fault | 3.57e-07 | 2,804,000 | 1.43e-07 | 7,004,000 |
| Visini (2025), 100 m, no pre-existing fault | 4.34e-08 | 23,050,000 | 1.88e-08 | 53,252,000 |

**The dominant control is not the model but the site classification.**
Whether a pre-existing fault is assumed beneath the site changes the
rate by **x137 at MS and x49 at SL** - far more than the choice between
Youngs and Visini, more than Youngs' own 85th/95th choice (x1.5), and
more than Visini's three magnitude-scaling relations (x1.0).

**Only the pre-existing-fault classification reproduces the
observations.** Without one, both sites return return periods of
2.8-53 Ma for displacements that were measured on the ground in 2016.
With one, MS gives 20,500 a and SL 128,000 a - the right order for
structures of this activity. This is the strongest single piece of
evidence in the section, because it is a prediction tested against a
field measurement rather than a model-versus-model comparison.

The classification is evidence-based here, not a default: MS (Mount
Serra) and SL (San Lorenzo) are pre-existing antithetic faults of the
Mount Vettore system and **both ruptured at the surface in 2016**. In
the Visini decision tree (their Figure 12) that is Case 1 - a Rank 1.5
structure beneath the site cannot be precluded - which activates
combinations A, B and C.

### A1. A configuration trap: the `case` key collision

`[calculation].case` has exactly ONE consumer, the Visini decision-tree
selector (`config_loader.py:315` -> `calculators.py:323` ->
`visini.py:56` -> `choose_combinations`). These INIs were seeded from
the IAEA benchmark, where `case = case3` denotes the exercise's own
SCENARIO name ("Norcia Case 3", both faults active). Consumed here it
silently meant decision-tree case3 and restricted every run to
combination A alone - the x137 understatement above.

Corrected to `case = case1` in all four job INIs, with the collision
documented in each. **Declaring the case does not require
`rank1p5_traces_file`**: that file is only the mechanism for DERIVING
the case from trace geometry when it is not known. Here it is known.

### A2. Why the labels say "pre-existing fault", not "A+B+C"

The combination letters are opaque outside the Visini paper, and the
physical content is exactly the decision tree's first question:

| figure label | decision tree | combinations | what it assumes |
|---|---|---|---|
| pre-existing fault | Case 1 | A + B + C | a Rank 1.5 structure beneath the site cannot be precluded |
| no pre-existing fault | Case 3 | A alone | such structures precluded beneath the site and within 1 km |

Combination A covers Rank 2, which the paper calls "simple" DR -
rupture with no pre-earthquake geologic or geomorphic evidence.
Combination C covers Rank 1.5/21/22/3, i.e. structures mappable BEFORE
the earthquake. So the label states what is being assumed about the
site, which is what a reader needs.

### A3. Combination C bypasses the along-strike term - verified correct

Combination C dominates the total, and it does so because it is not
multiplied by the along-strike Monte Carlo probability that suppresses
A and B. At MS (Mw 6.7, hanging wall, 500 m cell):

| | P_slice | x P_along | total |
|---|---|---|---|
| A | 0.1414 | 0.436 | 0.0617 |
| B | 0.0009 | 0.436 | 0.0004 |
| **C** | **0.2994** | **not applied** | **0.2994** |

This was initially flagged as a possible defect - `1 - prod(1 - P_i)`
appearing to combine quantities of different definition. **It is
correct.** Visini et al. (2025) p. 12: "The Monte Carlo procedure is not
applicable to DRs with Ranks 1.5, 21, 22, and 3. *If any of these DRs
are present beneath the site of interest, the probability calculated
based on the along-strike dimension is 1.*" A pre-existing mapped
structure is present beneath the site with certainty; only its
reactivation probability is in question. Their p. 11 states the same
restriction, and the Conclusion warns that Combination C "may dominate
the hazard".

### A4. The cell-size effect is narrower than it first appeared

Re-running the identical Visini chain at the 100 m IAEA site dimension
(its 500 m arm reproduces the tree curve to 1e-9, asserted in
`make_figures.py`):

| classification | MS: 500 m vs 100 m |
|---|---|
| no pre-existing fault (A alone) | **x8.2** |
| pre-existing fault (A+B+C) | **x1.04** |

**The cell convention only matters when combination A carries the
result.** The mechanism is the same one as A3: the along-strike Monte
Carlo term is strongly cell-dependent (P_along = 0.066 at 100 m vs
0.436 at 500 m, x6.6), and only A and B carry it. Combination C's cell
dependence lives solely in its intercept (a = 5.984 at 100 m vs 5.946
at 500 m), which is negligible - so once C dominates, the cell washes
out.

This supersedes the earlier reading of the x8.2 as a general "invisible
convention" result. It is real, but it is a property of the Rank 2
pathway, not of the model as a whole. Youngs' fixed 500 m cell is still
worth stating (`secondary_surf_rup/youngs2003.py` docstring), and the
IAEA exercise's 100 m site still cannot be answered by any published
tectonic normal-faulting distributed model except Visini.

### A5. Applicability ranges are declared too loosely (open item)

Visini gives ranges PER COMBINATION (pp. 14, 31): A ~0-6 km HW /
0-4 km FW; B ~0-1 km; C ~0-10 km HW / 0-8 km FW. The code declares a
single envelope (HW 10 / FW 8 km, i.e. C's) for all three, so no warning
is raised when A is used beyond its own limit.

MS is 7.59 km from MVFS and SL is 8.08 km from NFS - both beyond A's
6 km hanging-wall limit. **Impact on the reported results is
negligible** (A contributes 0.70% of the total at MS, 0.49% at SL,
because C dominates), but the "no pre-existing fault" comparison curves
ARE extrapolated at those distances and are not currently flagged.
Fixing this needs per-combination ranges plus a tracker that knows which
combination is active.

## B. Full-library spread

Design displacement at the 1e5 a return period (weighted q05 / mean /
q95 over the 60 branches; computed and printed, not rendered):

| site | q05    | mean     | q95     |
|------|--------|----------|---------|
| PF   | 2.04 m | 2.53 m   | 3.58 m  |
| MS   | 0      | 0.398 m  | 0.66 m  |
| SL   | 0      | 0.571 m  | 0.922 m |

- Both antithetic sites carry a **decimetre-to-metre** design
  displacement once the pre-existing-fault classification is applied.
  Under the mistaken case3 configuration MS had none at all and SL only
  0.067 m; the correction moves MS from "no design displacement" to
  0.398 m.
- SL's 0.571 m sits within a factor 1.8 of the ~1 m observed in 2016 -
  a second independent check on the classification, at a different site
  and a different displacement level from the MS tunnel.
- q05 is 0 at both sites because a majority-weight subset of branches
  puts the 1e5 a level below the curve entirely; the left clamp is the
  documented fractile convention, not a failure.

## C. Cross-trace profile (Y03 chain)

Profile through PF, D at 1e5 a, Petersen Gaussian W_p. Reported as
numbers only - the paper carries three figures and this is not one of
them; `make_figures.py` still computes it (the run prints these
values), it simply no longer renders:

- Peak on the trace 3.65 m; the corridor with D >= 1 m spans at least
  -40..+40 m at this profile's 20 m sampling, i.e. 80 m as measured and
  bounded above by the +-2 sigma = +-54 m (108 m) truncation width. The
  true value lies between the two; a 10 m sampling of the same profile
  would pin it, and is cheap (101 sites) if the number matters.
- The principal component exists ONLY inside +-2 sigma; outside it the
  total hazard IS the distributed component.
- Distributed HW/FW asymmetry: 0.493 m at -100 m vs 0.092 m at +100 m,
  0.414 m at -300 m vs 0.063 m at +300 m (factor 5-7 across the trace).

## D. Maps: what grid spacing can and cannot show

The principal corridor is 108 m wide (the +-2 sigma truncation of the
Petersen Gaussian). Measured cost of the map job over the
0.30 x 0.30 deg region, and what each spacing resolves:

| spacing | cell (E-W x N-S) | sites | run time | cells across the corridor |
|---------|------------------|-------|----------|---------------------------|
| 0.01 deg | 816 x 1112 m | 961 | ~5 min | 0.1 |
| 0.005 deg | 408 x 556 m | 3,938 | **19.6 min (measured)** | 0.3 |
| 0.002 deg | 163 x 222 m | 22,801 | ~1.9 h | 0.7 |
| **0.001 deg** | **82 x 111 m** | **90,601** | **~7.5 h** | **1.3** |
| 0.0005 deg | 41 x 56 m | 361,201 | ~30 h | 2.6 |
| 0.0002 deg | 16 x 22 m | 2,253,001 | ~8 days | 6.6 |

Cost is strictly linear in site count and independent of how near the
sites are to the fault: the regional 0.005 deg job ran at 6.7 s per
1000 sites per branch and the all-near-field zoom below at 6.9 s, so
the table extrapolates safely. The conclusion is that **0.001 deg is
not the answer**: seven and a half hours of compute buys a corridor
barely one pixel wide, so near-trace pixels would still be aliased and
still unusable as design values.

The corridor needs ~10-20 m cells, which is unaffordable over the
region. The near-trace structure is therefore characterised by the
cross-trace profile (section C) rather than by any map;
a fault-parallel strip product was trialled and dropped.

**The regional map is now PRINCIPAL ONLY** (panels (a)-(b) of
`fig_principal_fractiles`), and `fig_maps` has been dropped. The reason
is section A: the distributed term depends on the Visini site
classification, which is a PER-SITE geological judgement (case1 at MS
and SL, which sit on Rank 1.5 structures) and cannot be applied
uniformly across a grid without over-predicting everywhere away from
the antithetic faults. The principal term carries no such dependence,
so a principal-only map is well posed at every point. It is non-zero
only inside the +-2 sigma rupture-location band - 80 of 80 on-trace
sites, but only 8 of 992 grid sites - which is to say, on the fault.
**A trap worth recording** (met while trialling the strip product): the
obvious route to a fine near-trace grid - set a small `max_distance_km`
so only near-trace sites survive - is wrong. `max_distance_km` is read
once and used for two different things, the map site filter
(`site_builder.py:133`) *and* the rupture integration distance passed to
`FDHAContextMaker` (`calc/hazard.py:145`). Shrinking it silently
truncates the hazard integral. An explicit site list via
`[geometry].sites_csv` (absolute path required) is the safe route.

**A related subtlety, and a correction.** `build_hazard_map_sites` adds
**exact on-trace sites** on top of the filtered grid, resampled along
strike at the grid spacing (`is_trace = 1` in the map CSVs; the 0.01 deg
regional map has 80 of them against 992 grid sites). So the corridor
*peak* is computed exactly even on a coarse grid - it is only the
corridor *width and flanks* that a coarse grid cannot represent. The
4.13 m at the "PF pixel" of the regional map is a correctly computed
on-trace value that the plotting code spreads over a 1112 m cell.

**Displacement levels - and a retraction.** Two of the 80 on-trace
sites originally sat exactly at 10.0 m, the top of the first
`displacement_measure_levels`. That was first read as censoring of a
true value; it was not. Extending the levels to 30 m simply removed a
ceiling that had been MASKING an artefact: those sites then reported up
to 22.6 m of "distributed" displacement sitting on the fault trace,
~570x the principal term at the same point and eight orders of
magnitude above a grid site 59 m away.

The cause was Visini's ln(s) predictor being evaluated at r = 0, inside
the 5 m exclusion the paper itself declares (pp. 11, 20). The kernel now
gates the distributed term off inside a model's declared `r_min_km`
(commit 47cf5719), so those sites carry the principal contribution only.
On-trace displacement at 1e5 a now peaks at **3.64 m** (was 22.58 m),
no site exceeds 10 m (two did), and the distribution tightens to a 2.91 m
median with a 3.63 m 99th percentile. Levels still run to 30 m, and no
site sits near that ceiling.

Note this also affected **PF**, which lies 4.19 m from the trace and is
therefore inside the exclusion: its Visini distributed term is now gated
off. MS (7.6 km) and SL (2.4 km) are unaffected.

## E. Epistemic fractiles of the principal hazard (fig_principal_fractiles)

Principal displacement at 1e5 a on the 80 on-trace sites, as the
weighted MEAN hazard and as the 84th epistemic fractile:

| | median over sites | range |
|---|---|---|
| mean hazard | 2.91 m | 1.98 - 3.63 m |
| epistemic 84th | 3.76 m | 3.14 - 4.85 m |

The exported `rates_fractiles.h5` and the quantile CSVs carry TOTAL
rates only, with no principal column, so these had to be rebuilt from
`rates_principal` in each of the 60 branch files, weighted by branch
weight. The rebuilt mean reproduces the engine's exported
`displacement_map_mean.csv` exactly (2.91 m median, 1.98-3.63 m range),
which is a useful cross-check since the two paths share no code.

At PF specifically: mean 2.53 m, 16th 2.19 m, 84th 3.30 m - a +-1 sigma
epistemic band of -13% / +30%.

**Two different uncertainties, easily conflated.** These fractiles are
taken ACROSS LOGIC-TREE BRANCHES and are epistemic: the aleatory
variability is already integrated inside each individual branch curve.
They are not the aleatory percentiles of the displacement distribution
(section F), and the two must not be read as versions of one another.

Per-fault epistemic spread, measured on the branch curves (the source
tree has one branch, so all of this is FDHA-model uncertainty):

| | q50 | q84 | q84/q50 | CV |
|---|---|---|---|---|
| MVFS | 3.03 m | 3.89 m | 1.35 | 0.163 |
| NFS | 2.16 m | 3.66 m | **1.73** | **0.299** |

NFS carries 1.8x the epistemic dispersion of MVFS. The likely mechanism
is its lower mean magnitude (5.83 vs 6.70, section 0): the P_sr models
agree closely at Mw 6.4-7.0, where all approach 1, and diverge at
Mw 5.5-6.0. NOT yet verified directly - it is an inference from the
MFD shapes plus the measured spread.

## F. Design-displacement basis

Following Abrahamson and Liou (2025, BSSA 115, 2845-2856), for the
on-fault principal hazard (computed and printed by `make_figures.py`,
not rendered):

| quantity | value |
|---|---|
| rate of surface rupture at the site | 2.98e-04 /a (RP 3,351 a) |
| log-log hazard slope at 1e5 a | **-2.40** |
| equal hazard, mean curve, 1e5 a | 2.91 m |
| 85th percentile given rupture (their eq. 9) | 0.91 m |
| median given rupture | 0.20 m |

**Their alternatives do not apply here, and the check is the slope.**
The paper's three alternatives (upper fractile, risk-targeted, median
given rupture) exist for hazard curves that are FLAT at the design
return period, where the mean-hazard value misrepresents the epistemic
centre. Our slope is -2.40 - steeper than their high-activity case
(-2.07), nowhere near their moderate-activity case (-0.10). Their p.
2853 is explicit that the conditional approach "is intended to be used
for cases in which the hazard curve is flat at the selected return
period, **and not to reduce the displacement for cases in which the
hazard curve is steep**". So 2.91 m stands as the design basis, and the
0.91 m must not be substituted for it. The conditional percentiles are
plotted as a diagnostic - what displacement occurs GIVEN that rupture
happens - not as a design value.

The paper's closing recommendation that PFDHA reports should show
epistemic fractiles is met by section E.
