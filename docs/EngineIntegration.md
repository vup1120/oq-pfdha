# oq-engine Integration — Phase 1 Implementation Plan

Status: draft for implementation, 2026-07-20.
Target: incremental pull requests to the OpenQuake engine, staged through the
fork [vup1120/oq-engine](https://github.com/vup1120/oq-engine).

> This document was reconstructed from a fresh analysis of both codebases
> (oq-pfdha at `2f342af`, oq-engine fork at `8481936`); an earlier
> uncommitted draft of `docs/EngineIntegration.md` and the companion
> `docs/GEM_proposal_fdha.md` were lost with a previous working session and
> are not in git history.

## 1. Goal and strategy

Integrate the oq-pfdha fault-displacement hazard tool into the OpenQuake
engine, starting with the parts that fit the engine's existing architecture
with the least friction ("easy part first"):

1. **Model library** — port the FDHA model classes (`primary_surf_rup`,
   `secondary_surf_rup`, `primary_surf_displ`, `secondary_surf_displ`) into
   `openquake/fdha/` inside the engine.
2. **Distance metrics** — expose the FDHA distance parameters (`r`, `x/L`,
   `L`) through the engine's own distance machinery
   (`KNOWN_DISTANCES` / `get_distances` / surface methods), reusing the
   engine's geometry primitives so the two stacks stay coherent.
3. **Magnitude scaling laws** — move the magnitude–geometry–displacement
   scaling relations used by the FDHA models into
   `openquake.hazardlib.scalerel`, following the engine's existing MSR/ASR
   architecture.

Explicitly **not** in Phase 1 (see §7): an FDHA calculator inside the engine
job runner, job.ini plumbing, datastore outputs, logic-tree integration,
multi-fault reference lines (ECS/LCP), and the heavy-weight models
(Kuehn 2024, Chiou 2025, Lavrentiadis 2023, Visini 2025 rank-2 pipeline).

The standalone oq-pfdha tool remains the validation reference throughout:
every ported class gets tests pinned to oq-pfdha outputs, which are
themselves validated against fdhpy and the published-paper benchmarks
(`openquake/fdha/test/benchmark/`).

## 2. What is already in the fork (Phase 0 baseline)

The fork `vup1120/oq-engine` (branch `master`) already carries the seed of
the integration:

| Fork file | Content |
|---|---|
| `openquake/fdha/primary_surf_rup/base.py` | `BasePrimarySurfRup` with `get_prob(ctx)` — ctx-based signature, engine style |
| `openquake/fdha/primary_surf_rup/youngs2003.py` | `Youngs2003PrimarySR_ExC` / `_GB` / `_nBR` logistic P(surface rupture) models |
| `openquake/fdha/tests/primary_surf_rup/test_youngs2003.py` | unit tests |
| `openquake/hazardlib/geo/surface/base.py:243` | `BaseSurface.get_x_l_ratio(mesh)` → `(x/L, L_km)` computed from `_get_tor()` |
| `openquake/hazardlib/tests/geo/surface/base_surface_test.py` | `GetXLRatioTestCase` pinned to the Norcia trace (x/L ≈ 0.96, L ≈ 33.9 km) |

Not yet done: `x_l` is not wired into the distance dispatch, there is no
trace-distance parameter, no scalerel extensions, and only 3 of ~30 model
classes are ported.

**Branch hygiene**: these seeds currently sit on the fork's `master`. Before
opening PRs to `gem/oq-engine`, rebuild them as focused feature branches off
upstream master (one branch per PR in §6) so each PR is independently
reviewable.

## 3. Architecture mapping

How oq-pfdha concepts land in the engine:

| oq-pfdha (standalone) | oq-engine analogue | Phase 1 action |
|---|---|---|
| `BasePrimarySurfDispl` etc. with `get_prob()` | GSIM pattern: metadata-declaring classes with `compute(ctx, ...)` (`gsim/base.py:380`) | port as `openquake.fdha` classes, ctx-based `get_prob(ctx)` |
| `FDHAContext` / `FDHAContextMaker` (`calc/contexts.py`) | `RuptureContext` / `ContextMaker` (`hazardlib/contexts.py`) | reuse engine contexts; models read `ctx.mag`, `ctx.rake`, `ctx.rx`, … |
| `ctx.r`, `ctx.rx`, `ctx.x_L`, `ctx.L` | `KNOWN_DISTANCES` (`contexts.py:63`) + `get_distances` (`calc/filters.py:77`) + surface methods | add `rtor` + `x_l` distances; `L` as rupture parameter (§4.2) |
| `MULTIFAULT_REFERENCE_LINE` ('ecs'/'lcp'/'segments') | `surface.tor` (`Line` / `MultiLine` GC2) | Phase 2; note engine `MultiSurface.tor` ≈ fdha 'segments' |
| `openquake.fdha.scalerel` (WC1994 SRL/RLD/RW/RA + AD/MD, Leonard2010, Thingbaijam2017) | `openquake.hazardlib.scalerel` (BaseMSR/BaseASR; identical ABCs) | extend hazardlib classes (§5) |
| displacement levels in metres | IMT `Disp` (`hazardlib/imt.py:350`, already registered) | reuse as-is |
| model registry (package `__init__` imports) | `scalerel._get_available_class` scan; `valid.SCALEREL` (`valid.py:53`) | same scan pattern for `openquake.fdha` registries |
| FDHA logic tree XML, INI config | engine gsim logic tree, job.ini | Phase 2 |

The fdha ABCs `BaseMSR`/`BaseASR` are byte-identical to the engine's
(`hazardlib/scalerel/base.py`), and `FDHAContextMaker` was written against
the engine's `SiteCollection`/surface API — the tool was designed for this
port, so most classes move with signature changes only.

## 4. Workstream A — distance metrics (item 2)

### 4.1 Engine anatomy

- A GSIM declares `REQUIRES_DISTANCES ⊆ KNOWN_DISTANCES`
  (`hazardlib/contexts.py:63`); the `ContextMaker` computes the union of all
  requested distances per rupture.
- Every named distance is dispatched in ONE place:
  `openquake/hazardlib/calc/filters.py:55` (`get_dparam`) /
  `:77` (`get_distances`), which call surface methods
  (`get_min_distance`, `get_joyner_boore_distance`, `get_rx_distance`, …).
- All surface types expose `tor` (top-of-rupture line); `MultiSurface.tor`
  is a `MultiLine` and `rx`/`ry0` come from GC2 (`geo/multiline.py`).

### 4.2 What FDHA needs, and the mapping

| fdha ctx field | Meaning | Engine plan |
|---|---|---|
| `r` | horizontal distance (km) to the surface-rupture trace | **new** distance `rtor`: point-to-polyline distance to `surface.tor` (NOT `rjb`, which is the distance to the whole surface projection and is 0 above a dipping fault plane; NOT \|`rx`\|, which is the unclamped GC2 T and differs beyond the rupture ends) |
| `rx` | signed perpendicular distance, + = hanging wall | **exists** (`rx`, GC2). Adopt the engine convention in-engine; see §4.4 |
| `x_L` | normalized along-strike position ∈ [0, 1] | **new** distance `x_l`, backed by the already-merged `BaseSurface.get_x_l_ratio` |
| `L` | rupture trace length (km) | **new rupture parameter** `length` (not a site-dependent distance); fill from `surface.tor` length in the ContextMaker, same source as `get_x_l_ratio`'s `L_km` return |
| `mag`, `rake`, `dip`, `ztor` | rupture params | already in `RuptureContext._slots_` (`contexts.py:1727`) |
| `vs30` | site param (Kuehn 2024 only) | already a site parameter; Phase 2 |

### 4.3 Implementation steps

1. **`BaseSurface.get_tor_distance(mesh)`** (name open: `rtor` vs `rtrace`)
   in `geo/surface/base.py`: clipped point-to-polyline distance to
   `self.tor`, computed in the `OrthographicProjection` local km frame.
   Share the projection sweep with `get_x_l_ratio` — both need the same
   per-segment projection; refactor the two to a common private helper that
   returns `(dist, x_ratio, L)` in one pass, mirroring oq-pfdha's
   `RuptureDistanceCalculator._core` memoization
   (`openquake/fdha/calc/utils/rupture_distance.py`).
   `MultiSurface` override: minimum over surface-reaching section `tor`s
   (fdha 'segments' semantics) — can land in Phase 1 if cheap, else Phase 2.
2. **Wire the dispatch**: add `rtor` and `x_l` to `KNOWN_DISTANCES`
   (`contexts.py:63`) and branches in `get_dparam`/`get_distances`
   (`calc/filters.py`).
3. **Rounding guard**: the context-collapsing logic rounds
   `KNOWN_DISTANCES` values to ≥ 1 km (`contexts.py:254-259`); `x_l` ∈
   [0, 1] must be exempted the way `rx` already is (round, don't ceil), or
   excluded from collapsing entirely.
4. **`length` rupture parameter**: accept `'length'` in
   `REQUIRES_RUPTURE_PARAMETERS`; ContextMaker fills it from the surface
   (`tor` length), consistent with `get_x_l_ratio`.
5. **Tests**:
   - extend `GetXLRatioTestCase` with `rtor` cases (site off the trace end,
     hanging wall/footwall, > maximum-distance);
   - a cross-tool parity test (lives in oq-pfdha, §8): engine
     `rtor`/`x_l`/`rx` vs oq-pfdha `RuptureDistanceCalculator` on shared
     fixtures (Norcia trace, IAEA benchmarks) at fine
     `rupture_mesh_spacing`, with a documented tolerance.

### 4.4 Coherence policy (the deliberate divergences)

oq-pfdha computes `r`/`rx`/`x_L` against the **original NRML trace**
retained at parse time, because the engine's resampled mesh top edge
corner-cuts wiggly traces (see the header of
`openquake/fdha/calc/utils/rupture_distance.py`). Inside the engine we do
NOT carry a parallel geometry stack:

- In-engine FDHA consumes the engine's own primitives — `tor`,
  GC2/`MultiLine`, `OrthographicProjection` — which oq-pfdha already
  delegates to for the underlying math. Only the trace source changes
  (original NRML trace → engine `tor`).
- `SimpleFaultSurface.tor` is mesh row 0 through `keep_corners(1.0)`
  (`geo/surface/simple_fault.py:72`), so it is `rupture_mesh_spacing`-
  dependent. The parity tests quantify the difference; at ≤ 1 km spacing it
  is expected within the 1 km corner tolerance.
- If parity shows material drift for production traces, the Phase 2 fix is
  engine-side: retain the original trace on the surface at construction
  (small change in the simple-fault surface builder), not an fdha-side
  workaround.

## 5. Workstream B — magnitude scaling laws (item 3)

### 5.1 What the tool actually uses

Three relations in `openquake/fdha/scalerel/`, all extending the same
BaseMSR/BaseASR ABCs as the engine, plus PFDHA-specific regressions:

| Class | Beyond the engine's version | Consumed by |
|---|---|---|
| `WellsCoppersmith1994` | SRL/RLD/RW forward+inverse; **AD/MD displacement** forward+inverse, per style, with sigmas | `Youngs2003PrimaryFD` convolution (AD & MD); map principal-band site placement (SRL) |
| `Leonard2010` (interplate dip-slip) | bilinear L(m); W from L capped at 20 km; **AD = 1.7e-5 · L** | AD-based normalization options |
| `Thingbaijam2017` | width/length regressions; **AD** | AD-based normalization; `Mammarella2024PrimarySR` W(m) |
| (in-model tables) Leonard 2014 interplate/SCR W(m) | `m = a + 2.5·log10 W` | `Mammarella2024PrimarySR` (MSR codes 0/1/2 in `primary_surf_rup/mammarella2024.py` TAB1) |

### 5.2 Engine architecture and the precedent

- `openquake.hazardlib.scalerel` auto-discovers every `BaseMSR` subclass in
  the package (`scalerel/__init__.py::_get_available_class`), and
  `valid.SCALEREL` / `valid.mag_scale_rel` (`valid.py:53,1089`) expose them
  to job configuration and source models by class name. New classes are
  picked up with zero registry work.
- **The engine already accepts per-class methods beyond the ABC**:
  `ThingbaijamInterface.get_median_length/get_median_width` and the crustal
  classes' `get_std_dev_width` (`scalerel/thingbaijam2017.py:63-129`). This
  is the pattern to extend — no ABC changes needed.

### 5.3 Implementation steps

1. **`wc1994.py`**: add the FDHA regressions to `WC1994` directly
   (additive, non-breaking, Thingbaijam precedent):
   `get_median_srl`, `get_median_rld`, `get_median_width`,
   `get_average_displacement`, `get_maximum_displacement` — each with the
   per-style coefficient tables and `return_sigma` option from
   `openquake/fdha/scalerel/wc1994.py`. If reviewers prefer the upstream
   class untouched, fall back to a `WC1994_FD(WC1994)` subclass in the same
   module — decide in the PR discussion.
2. **`leonard2014.py`**: add `get_median_width` / `get_std_dev_width` to
   the interplate and SCR classes (`log10 W = (m − a)/2.5`; a = 3.63/3.88
   dip-slip/strike-slip interplate, 4.14/4.22 SCR — the TAB1 constants of
   `mammarella2024.py`).
3. **`thingbaijam2017.py`**: add `get_median_width` to the crustal classes
   (normal −0.829/0.323, reverse −1.669/0.435, strike-slip −0.543/0.261 —
   they already have `get_std_dev_width`), and `get_average_displacement`
   where oq-pfdha provides it.
4. **`leonard2010.py`**: port oq-pfdha's interplate dip-slip `Leonard2010`
   (bilinear L(m), W cap, AD) as a new class beside the existing
   `Leonard2010_SCR*` (e.g. `Leonard2010_Interplate`).
5. **Consumers switch to engine idiom**: when `Mammarella2024PrimarySR` is
   ported (§6, PR-4), replace the `MSR ∈ {0, 1, 2}` integer codes with a
   scalerel instance parameter (`msr=Leonard2014_Interplate()`), resolved
   through `valid.mag_scale_rel` — the same way sources consume
   `magScaleRel`. Keep a code→class mapping in the class docstring for
   traceability to the paper. Same for `Youngs2003PrimaryFD`'s
   `scaling_model="WC1994"` string.
6. **Tests**: unit tests pinned to both the paper values and the oq-pfdha
   outputs for every added coefficient set.

## 6. Workstream C — model library expansion (item 1), and PR sequencing

### 6.1 Porting rules

1. Engine copyright header (GEM), as in the already-ported seed files.
2. `get_prob(ctx, ...)`: read `ctx.mag`, `ctx.rake`, `ctx.rx`, `ctx.rtor`,
   `ctx.x_l`, `ctx.length` (duck-typed: engine recarray context or any
   object with those attributes — keeps classes unit-testable without a
   ContextMaker).
3. Preserve the declarative class attributes
   (`DISPLACEMENT_DEFINITION`, `DISPLACEMENT_COMPONENT`,
   `APPLICABILITY_RANGE`, `SECONDARY_PIPELINE`, …) verbatim — they are the
   FDHA analogue of GSIM metadata and the future engine calculator will
   route on them. `MULTIFAULT_REFERENCE_LINE` is kept but inert until
   Phase 2.
4. One test module per model, values pinned to oq-pfdha (which pins to
   fdhpy / papers).
5. Per-subpackage registry built by subclass scan
   (`scalerel._get_available_class` pattern), not hand-maintained imports.

### 6.2 Port order (dependency-driven) and PR slicing

**PR-1 — scalerel extensions** (§5; no dependencies, smallest, immediately
upstreamable). This is the recommended first PR.

**PR-2 — FDHA distances in hazardlib** (§4): `rtor` + `x_l` + `length`
wiring and tests. Includes moving the fork's existing `get_x_l_ratio` work
into this reviewable branch. Independent of PR-1.

**PR-3 — magnitude-only rupture-probability models** (depends on nothing;
can go in parallel with PR-1/2): complete `openquake/fdha/primary_surf_rup/`
with the remaining P(surface rupture | m) models, all trivial ports —
`wells_coppersmith1993`, `takao2013`, `moss2013`, `moss_ross2011`,
`pizza2023`, `yang2021`, `fixed` — plus the subclass-scan registry and
tests. (`mammarella2024` waits for PR-1: it needs the W(m) scalerels.)

**PR-4 — distance- and scaling-dependent models** (after PR-1 + PR-2):
- `primary_surf_rup/mammarella2024.py` (scalerel-based CPSR);
- `secondary_surf_rup/`: `youngs2003`, `petersen2011`, `takao2013`,
  `takao2014`, `moss2022`, `ferrario2021`, `rodriguez2023`, `fixed`
  (all need `rtor`/`rx`);
- `primary_surf_displ/`: `youngs2003` (needs WC1994 AD/MD + `x_l`),
  `petersen2011`, `moss_ross2011`, `moss2022`, `moss2024`, `takao2013`;
- `secondary_surf_displ/`: `youngs2003`, `petersen2011`, `takao2013`,
  `moss2022`.

PR-4 will likely split into 2–3 PRs by slot to stay reviewable (target
< ~1000 added lines each, tests included).

**Deferred to Phase 2** (heavy data files, extra deps, or special
pipelines): `kuehn2024` (data tables), `chiou2025`, `lavrentiadis2023`,
`visini2025` (both slots; rank-2 Monte Carlo pipeline + segments-r
calibration).

## 7. Out of scope for Phase 1 (Phase 2 backlog)

- FDHA hazard calculator inside the engine (job runner, `job.ini`
  parameters, datastore outputs, exports) — Phase 1 delivers the library
  layers the calculator will stand on.
- FDHA logic tree ↔ engine logic tree integration.
- Multi-fault reference lines (ECS / LCP) and per-model
  `MULTIFAULT_REFERENCE_LINE` routing; `MultiSurface` FDHA distances beyond
  the simple 'segments' minimum.
- Original-NRML-trace retention on engine surfaces (only if parity tests
  demand it, §4.4).
- Near-field floor / principal-vs-distributed weighting (`W_p`,
  `r_threshold_km`) — calculator-level concerns.
- The deferred models listed in §6.2.

## 8. Validation strategy

1. **Unit level (in-engine)**: every ported class ships tests pinned to
   oq-pfdha outputs; scalerel coefficients double-pinned to paper values.
2. **Cross-tool parity harness (in oq-pfdha)**: a test module in this
   repository that imports the engine fork's `openquake.fdha` and asserts
   `get_prob` equality (and distance-metric agreement) against the local
   implementations over a grid of (mag, style, r, x/L) and over the
   existing benchmark fixtures (IAEA, Norcia, Mammarella CPSR golden CSV).
   This is the regression net while the two trees coexist.
3. **Benchmarks**: the oq-pfdha benchmark suite
   (`openquake/fdha/test/benchmark/`) remains the source of truth; any
   in-engine deviation must be explained by a documented convention change
   (§4.4) or fixed.

## 9. Immediate next actions

1. Create feature branch `fdha-scalerel` off upstream `gem/oq-engine`
   master in the fork; implement §5 (PR-1).
2. Create feature branch `fdha-distances`; implement §4 (PR-2), folding in
   the existing `get_x_l_ratio` + test.
3. Create feature branch `fdha-surf-rup-models`; implement §6 PR-3.
4. Add the parity-test harness to oq-pfdha (§8.2) so PRs 1–3 are verified
   against the tool from day one.
