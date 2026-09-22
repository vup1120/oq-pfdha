# oq-engine Integration Plan — FDHA

Status: agreed direction, 2026-09-14.
Target: `gem/oq-engine` (checked out at `~/oq-engine`, master `8d5712fe7f`).
Supersedes the earlier Phase 1 draft of this file.

## 0. Agreed end state

**All core FDHA logic lives in the OpenQuake engine.** The standalone
`oq-pfdha` repository becomes a thin consumer that imports the engine and adds
only the standalone CLI, examples, docs and validation material.

```
FINAL OWNERSHIP
  engine (gem/oq-engine)                    oq-pfdha (this repo, consumer)
  ├─ openquake/pfd/             models      ├─ CLI (fdha) + examples
  ├─ openquake/hazardlib/       distances   ├─ benchmark/validation suites
  │    scalerel/, calc/…                    ├─ docs + web GUI
  └─ openquake/calculators/displacement.py          └─ imports openquake.pfd from engine
```

Consequences that must be decided/fixed first:

1. **Namespace collision — resolved.** The engine's library was renamed from
   `openquake.fdha` to **`openquake.pfd`** (engine commit `e28c8ad91c`), so the
   import name `openquake.fdha` is free. oq-pfdha keeps `openquake.fdha` as its
   CLI/consumer package and later imports the engine's `openquake.pfd`; no
   oq-pfdha rename is needed (the old PR-0 is dropped).
2. **Single source of truth for models.** The engine's older FDHA seeds
   (`primary_surf_rup/youngs2003.py` ExC/GB/nBR, `primary_surf_displ/youngs2003.py`
   AD/MD, committed 2026-02..06) are superseded by oq-pfdha's newer library
   (4 categories, 34 model modules, declarative contracts, 2026-07). They are
   **deleted and replaced**, never merged with the oq-pfdha versions.
3. **The engine keeps its own general infrastructure; oq-pfdha keeps none of
   its `calc/`, `logic_tree/`, `main.py` or output writers.** Those exist only
   because oq-pfdha had to be self-contained; in the end state the engine's
   source reader, logic tree, contexts, datastore, stats and exports are used.

## 0.1 Settled decisions (2026-09-14)

All decisions below are confirmed and are implemented by the PRs as written.

| # | Decision | Choice |
|---|---|---|
| D1 | Package names | engine library = `openquake.pfd` (renamed, done `e28c8ad91c`); oq-pfdha keeps `openquake.fdha` (CLI/consumer) and imports `openquake.pfd`; no oq-pfdha rename |
| D2 | Trace geometry source | engine `surface.tor` only; no retained original trace unless parity forces it |
| D3 | Datastore datasets | reuse `hcurves-*` / `hmaps` keyed by IMT `Disp`; `hcurves-rlzs` has the full logic-tree cardinality `R = len(sm_rlzs) * gsim_paths * pfd_paths` including the PFD realizations (PFD branches are realization dimensions, not averaged away) |
| D4 | IMT / storage | reuse IMT `Disp`; store annual rates internally (as the classical path) |
| D5 | FDHA logic-tree schema | keep oq-pfdha's XML with `<logicTreeBranchSet>` directly under `<logicTree>`; parsed by `hazardlib.pfd_lt.PFDLogicTree` (moved out of `gsim_lt` and renamed from `FdhaLogicTree`). The legacy `<logicTreeBranchingLevel>` wrapper is accepted too (oq-pfdha has migrated its files) |
| D6 | CLI ownership | engine exposes the FDHA calculation mode; oq-pfdha keeps a thin `fdha` wrapper. The mode is `displacement` (D13) |
| D7 | Canonical model API | oq-pfdha explicit `get_prob(d, mag, rx, r, …)`; one engine ctx adapter |
| D8 | Numeric parity | 1e-12 rel. on rates for identical geometry; documented tolerance where `tor` differs |
| D9 | Site/grid builder | shared helper in hazardlib/commonlib |
| D10 | Output buckets | keep `principal`/`distributed` components + aggregate single-bucket routing |
| D11 | oq-pfdha library | fully removed at PR-9; no vendored fallback |
| D12 | `get_poes` bypass | confirmed — FDHA kernel writes `MapArray` directly |
| D13 | Mode / function naming | `calculation_mode = displacement` (renamed from the provisional `fdha_classical`), `readinput.get_pfd_lt`, `hazardlib.pfd_lt.PFDLogicTree`, library `openquake.pfd`. The mode requires `use_rates = true` and `disagg_by_src = true` (the kernel is per source and yields annual rates). The PFD logic tree is passed through a dedicated `pfd_logic_tree_file` key (`readinput.get_gsim_lt` returns the trivial one-branch `PFDGMPE` tree instead) and to `FullLogicTree` as `extra_lt` (mutually exclusive with the amplification tree), and displacement levels use `intensity_measure_types_and_levels = {"Disp": [...]}` (no `displacement_measure_levels`) |

## 1. Ground truth (both trees, at plan time)

`~/oq-engine` (gem master `8d5712fe7f`) already has:

- `openquake/pfd/{primary_surf_rup,primary_surf_displ,utils.py}` and
  `openquake/pfd/tests/` — older, ctx/arg-based Youngs 2003 primary SR + FD
  (renamed from `openquake/fdha` in `e28c8ad91c`).
- `openquake/hazardlib/scalerel/wc1994.py` — SRL/RLD/RW (no displacement
  methods); `thingbaijam2017.py` — median/std width; `leonard2010/2014`.
- `openquake/hazardlib/geo/surface/base.py::get_x_l_ratio` (x/L + L).
- IMT `Disp` registered in `hazardlib/imt.py`.
- **Missing:** secondary models, `rtor`/`x_l`/`length` distance wiring, FDHA
  logic tree, FDHA calculator, exports.

`~/oq-pfdha` (`main` `2f342af`, 2026-07-19) has the authoritative model library
and semantics:

- `openquake/fdha/{primary_surf_rup,primary_surf_displ,secondary_surf_rup,secondary_surf_displ}/`
  with `DISPLACEMENT_DEFINITION`, `DISPLACEMENT_COMPONENT`,
  `APPLICABILITY_RANGE`, `SECONDARY_PIPELINE`, `MULTIFAULT_REFERENCE_LINE`.
- `openquake/fdha/calc/hazard.py` — the principal/distributed/`W_p` kernel.
- `openquake/fdha/calc/{contexts,model_adapter,location_weight,visini}.py`.
- `openquake/fdha/logic_tree/` — FDHA LT parsing, enumeration, aggregation.
- 4 model slots + calc-param slot (`fdhaCalcRSigma`), benchmark/parity suites.

## 2. Target architecture (engine)

```
oq-engine
├─ openquake/hazardlib/
│  ├─ geo/surface/base.py          # get_rtor(mesh); shared (rtor, x/L, L) sweep
│  ├─ calc/filters.py              # dispatch rtor, x_l
│  ├─ contexts.py                  # KNOWN_DISTANCES += rtor, x_l; length rup param
│  ├─ scalerel/                    # AD/MD relations + missing widths
│  └─ calc/displacement.py         # NEW: FDHA rate kernel
├─ openquake/pfd/                  # CANONICAL model library (lifted from oq-pfdha)
│  ├─ primary_surf_rup/  primary_surf_displ/
│  ├─ secondary_surf_rup/  secondary_surf_displ/
│  ├─ base.py                      # 4 ABCs + declarative contracts
│  ├─ adapter.py                   # ctx recarray -> model call (engine facade)
│  └─ registry.py                  # subclass-scan registries per slot
└─ openquake/calculators/displacement.py   # @base.calculators.add('displacement')
```

Reused engine layers (no duplication allowed):

| Need | Engine component reused |
|---|---|
| Source model / NRML reading | `hazardlib.source_group.read_csm`, `source_group.SourceGroup` |
| Source-model logic tree | `hazardlib.logictree` (already handles arbitrary `uncertaintyType`, `applyToBranches`, weights) |
| Jobs / params / validation | `commonlib.oqvalidation`, `commonlib.readinput` |
| Context building | `hazardlib.contexts.ContextMaker`, `get_cmakers` |
| Rate maps / storage | `hazardlib.map_array.MapArray`/`RateMap`, `calculators.base._store` |
| Mean/quantiles | `hazardlib.stats.mean_curve` / `quantile_curve` |
| Datastore, exports, views, plots | `calculators.base`, `calculators.export`, `calculators.views` |
| Displacement IMLs | IMT `Disp` (`hazardlib/imt.py`) |

## 3. Workstreams

### A. Library consolidation (PR-1)
- Engine library namespace is `openquake.pfd` (rename already done,
  `e28c8ad91c`); oq-pfdha keeps `openquake.fdha` and imports `openquake.pfd`.
- Lift oq-pfdha's four model packages + base ABCs + contracts into the engine's
  `openquake/pfd/`; **delete** the older seeds.
- **Unify the model API.** oq-pfdha's `get_prob(d, mag, rx, r, …)` is canonical
  (it is the benchmarked one); the engine seed's `get_prob(d, x_l, mag, rake)`
  and `get_prob(ctx)` are dropped. The ctx→call translation lives in
  `openquake/pfd/adapter.py`, written once.
- Add subclass-scan registries per slot (the `scalerel._get_available_class`
  pattern) so logic-tree class names resolve with no hand-maintained imports.

### B. Geometry & distances (`hazardlib`)
- `get_rtor(mesh)`: clipped point-to-polyline distance to `surface.tor` in the
  engine's local km frame. **Not** `rjb`, not `|rx|`.
- Refactor `get_rtor` and the existing `get_x_l_ratio` to share one projection
  sweep returning `(distance, x/L, L)` (mirror oq-pfdha's
  `RuptureDistanceCalculator._core` memoization).
- Add `rtor` and `x_l` to `KNOWN_DISTANCES` + `get_dparam`/`get_distances`;
  add a `length` rupture parameter filled from `surface.tor`.
- Exempt `x_l ∈ [0,1]` from the context-collapse rounding (as `rx` is).
- `MultiSurface`: `rtor` = min over surface-reaching section `tor`s
  (oq-pfdha "segments" semantics).
- **Trace-source policy (D2, decided):** in-engine FDHA consumes
  `surface.tor` only. The divergence from oq-pfdha's original-NRML-trace
  behavior is documented and bounded by the parity harness; retaining original
  traces on engine surfaces is an escalation only if parity proves material.

### C. Magnitude–displacement scaling (`hazardlib/scalerel`)
- Add AD/MD regressions to `WC1994` (additive; Thingbaijam's extra methods are
  precedent) or a `WC1994_FD` subclass.
- Add `Leonard2014` interplate/SCR widths; add interplate dip-slip `Leonard2010`.
- Ensure `Thingbaijam2017` covers `Mammarella2024PrimarySR`.
- Consumers take scalerel instances (resolved via `valid.mag_scale_rel`)
  instead of integer "MSR codes".

### D. FDHA hazard kernel (`hazardlib/calc/displacement.py`)
- Engine-idiom pure function: `(contexts, models, W_p params, imls) -> rmap`.
- Encodes oq-pfdha's `calc/hazard.py` semantics verbatim (ported once):
  - `λ_principal = rate·P_sr·P_fd_primary·W_p`,
    `λ_distributed = rate·P_sr·P_dist_combined·G`;
  - two `W_p` paths (σ=0 boxcar + `G=1−W_p`; σ>0 pinned Gaussian + `G=1`);
  - aggregate-definition single-bucket routing;
  - Visini combined pipeline where declared.
- **Bypass `get_mean_stds`/`get_poes` (D12, confirmed)**: FDHA yields
  exceedance probabilities directly; update a `MapArray` and reuse all
  downstream machinery.

### E. FDHA logic tree + job parameters
- Reuse the engine NRML scaffolding via `hazardlib.gsim_lt.FdhaLogicTree`;
  map `fdhaPrimarySRModel`, `fdhaPrimaryFDModel`,
  `fdhaSecondarySRModel`, `fdhaSecondaryFDModel` and `fdhaCalcRSigma` onto the
  four slots + calc-param slot (mirror oq-pfdha's `logic_tree/types.py`).
- `oqvalidation.py`: add `'displacement'` to `ALL_CALCULATORS`; the PFD
  logic tree is passed through a dedicated `pfd_logic_tree_file`; declare
  `r_threshold_km`, `r_sigma_km`; reject scalar-vs-branch
  `r_sigma` conflicts up front. Displacement levels use
  `intensity_measure_types_and_levels = {"Disp": [...]}` (D3/D4).
- Reuse IMT `Disp` for the IML container, storing annual rates (D3/D4).

### F. Calculator (`openquake/calculators/displacement.py`)
- `@base.calculators.add('displacement')`, building on `base.HazardCalculator` /
  `preclassical` to inherit source reading, csm, realizations, datastore,
  checkpointing.
- Reuse `get_cmakers`/`read_full_lt_by_label` for branch enumeration; feed the
  FDHA kernel instead of GSIM poes.
- Store under the engine's existing `hcurves-rlzs` / `hmaps` dataset names
  keyed by IMT `Disp` (D3), using `MapArray`/`RateMap`, so existing stats,
  exports, views and plots work unchanged.
- Map mode: reuse the engine's curve→map inversion and a hazardlib/util site
  builder (port oq-pfdha's `site_builder` as a helper, not a second grid engine).

### G. Exports, views, plots
- With D3/D4 the existing `hcurves`/`hmaps` exporters, views and plots apply
  unchanged; only FDHA-specific additions (params view, displacement labels)
  need new code, registered through the existing `export` dispatch.

### H. Multi-fault reference lines + heavy models
- `MULTIFAULT_REFERENCE_LINE` routing (ECS/LCP/segments) via `MultiSurface`;
  keep model metadata inert until implemented.
- Kuehn 2024, Chiou 2025, Lavrentiadis 2023, Visini 2025 rank-2 pipeline; data
  tables shipped under `openquake/pfd/**/data/` as in oq-pfdha.

### I. Strip oq-pfdha of duplicated logic (final phase)
- Remove oq-pfdha's model packages, `calc/`, `logic_tree/`, output writers once
  the engine's calculator is the functional path. oq-pfdha keeps: CLI, examples,
  benchmark/parity suites, docs, web GUI — all importing the engine.
- Update `pyproject.toml` to depend on the engine (or gel) and re-point the
  `fdha` console script at the engine calculator.

## 4. Migration map (oq-pfdha → engine)

| oq-pfdha | Engine destination | Action |
|---|---|---|
| `primary_surf_rup/`, `primary_surf_displ/`, `secondary_surf_rup/`, `secondary_surf_displ/` | `openquake/pfd/**` | lift; API canonical |
| model base classes + contracts | `openquake/pfd/base.py` | lift |
| `calc/model_adapter.py` | `openquake/pfd/adapter.py` | lift, engine ctx facade |
| `calc/location_weight.py`, `calc/hazard.py` kernel | `openquake/hazardlib/calc/displacement.py` | port semantics |
| `calc/utils/{rupture_distance,segments,interpolation,probability,lcp,ecs}.py` | `hazardlib/calc/`, `hazardlib/geo/` | fold into engine (distances, map inversion) |
| `calc/{contexts,config_loader,calculators}.py` | — | **discard**; use engine ContextMaker/oqvalidation/calculators |
| `logic_tree/**` | engine `logictree` + `calculators/displacement.py` | **discard** LT engine; keep branch-set mapping |
| `scalerel/**` | `hazardlib/scalerel/` | lift/fold |
| `main.py`, `logic_tree/io.py` writers | engine calculator + exports | **discard** |
| `demo/`, `webgui_demo/`, `test/benchmark` | stays in oq-pfdha | consumer |

## 5. PR sequence

| PR | Scope | Depends on | Acceptance |
|---|---|---|---|
| ~~PR-0~~ | ~~Rename oq-pfdha package~~ — **dropped**; the engine library was renamed to `openquake.pfd` instead (`e28c8ad91c`) | — | done |
| ~~PR-1~~ | ~~Engine library consolidation: lift oq-pfdha models into `openquake/pfd`, delete old seeds, adapter, registries~~ — **done** (`cafbc8dce7`) | — | oq-pfdha model tests pass against in-engine classes |
| ~~PR-2~~ | ~~scalerel AD/MD + widths (Workstream C)~~ — **done** (`4b54e8462c`): folded into `hazardlib/scalerel`, `openquake/pfd/scalerel` deleted, `width_model` scalerel instances | — | pinned to papers + oq-pfdha outputs (`hazardlib/tests/scalerel/pfd_scalerel_test.py`) |
| ~~PR-3~~ | ~~`rtor` + `x_l` + `length` (Workstream B)~~ — **done** (`91b09801a3`) | — | distance parity vs oq-pfdha (Norcia/IAEA) |
| ~~PR-4~~ | ~~Remaining SR/FD models (distance-dependent)~~ — **done** (`c3eaa2e1e5`): adapter wired to the engine `rtor`/`x_l`/`length` context, style from rake, wiring + parity tests | PR-2, PR-3 | model tests + parity |
| **PR-5** | FDHA logic tree + oqparam params (Workstream E) — logic tree **done** (`bfbf32886c`, later moved to `hazardlib.pfd_lt.PFDLogicTree` with h5 serialization and Monte-Carlo sampling): oq-pfdha schema + filters; **oqparam declarations deferred to PR-6** | PR-1, PR-4 | branch enumeration matches oq-pfdha manifest |
| **PR-6** | `hazardlib/calc/displacement.py` + `calculators/displacement.py` curve mode under `hcurves-*` (D3/D4/D12) — kernel **done** (`27687f1bdf`): `location_weight` + rate kernel, parity-verified; oqparam surface + logic-tree entry point **done**: `readinput.get_pfd_lt` (reads `pfd_logic_tree_file`) dispatched from `get_gsim_lt` when `calculation_mode == 'displacement'`; **calculator remaining** | PR-3, PR-5 | reproduces `hazard_curve_minimal` within tolerance |
| **PR-7** | map mode + exports/views/plots | PR-6 | reproduces `hazard_map_minimal` |
| **PR-8** | heavy models + multi-fault | PR-6 | benchmarks pass |
| **PR-9** | Strip oq-pfdha duplicated logic (Workstream I) | PR-6..8 | oq-pfdha runs on engine imports only |

## 6. Validation strategy

1. **In-engine unit tests** for every ported model, pinned to oq-pfdha outputs
   (themselves pinned to fdhpy and published papers).
2. **Parity harness in oq-pfdha** (extend
   `test/integration/test_engine_parity.py`): import the engine's
   `openquake.pfd`, assert `get_prob`/distance equality over a grid and over
   the benchmark fixtures — the regression net while both trees coexist.
3. **End-to-end goldens**: the two public examples must reproduce oq-pfdha's
   `aggregate_hazard.csv` / displacement maps — byte-exact where the geometry
   source is identical, else within a documented tolerance.
4. **Benchmarks**: oq-pfdha's `test/benchmark/` stays source of truth; any
   deviation must be explained by a documented convention change.

## 7. Residual risks / open questions

Decisions D1–D13 are settled (§0.1); the items below remain engineering risks.

- **Logic-tree reuse limits** — verify four-slot chains with `applyToBranches`
  enumerate identically to oq-pfdha; add a thin mapping only if engine LT
  cannot express them.
- **`r_sigma` two-path fidelity** — preserve bit-for-bit against oq-pfdha's
  frozen `curve_explicit` / `map_default` fixtures.
- **`hcurves`/`hmaps` compatibility** — confirm no existing consumer hard-codes
  a non-`Disp` IMT assumption that the D3 reuse would violate.
- **`tor` divergence magnitude** — if the parity harness shows more than corner
  tolerance, escalate per D2.

## 8. Immediate next action

PR-0..PR-5 are done on the engine's `oq-integration` branch (which is
currently identical to `master`) and **PR-6** is half-done:
`openquake/hazardlib/calc/displacement.py` holds the FDHA rate kernel
(`location_weight` verified bit-identical to oq-pfdha; the
principal/distributed rate core and the aggregate single-bucket routing), and
`openquake/hazardlib/pfd_lt.py::PFDLogicTree` (moved out of `gsim_lt.py` and
renamed from `FdhaLogicTree`, with h5 serialization and Monte-Carlo sampling)
enumerates the oq-pfdha end branches.

The logic-tree entry point is now wired: `readinput.get_pfd_lt(oqparam)` reads
the dedicated `pfd_logic_tree_file` input and builds the `PFDLogicTree`,
while `readinput.get_gsim_lt` returns the trivial one-branch `PFDGMPE` logic
tree when `calculation_mode == 'displacement'` (D13). The oqparam
surface is also done: `'displacement'` in `ALL_CALCULATORS`, `r_threshold_km`
(default 0.1) and `r_sigma_km` (default 0.0), displacement levels via
`intensity_measure_types_and_levels = {"Disp": [...]}`, a default
`maximum_distance = 10` km, and mandatory `use_rates = true` and
`disagg_by_src = true`. `check_gsim_lt` reads `pfd_logic_tree_file`, and
`is_valid_maximum_distance` / `get_input_files` have `displacement` handling.

**`get_full_lt` is now attacked.** FDHA has no GSIMs, so `FullLogicTree` is
built with a trivial one-branch GSIM logic tree holding one no-op `PFDGMPE`
per source TRT (`openquake/pfd/gsim.py`); the dummy declares the FDHA context
requirements (`REQUIRES_DISTANCES = {rtor, x_l, rx}`,
`REQUIRES_RUPTURE_PARAMETERS = {mag, dip, rake, length}`) so `ContextMaker`
builds contexts with exactly those fields. The PFD logic tree is passed to
`FullLogicTree` as `extra_lt` (the same slot used by the amplification logic
tree; the two are mutually exclusive) and serialized generically via its own
`__toh5__`/`__fromh5__`. This also fixed two engine issues: `set_distances`
treated `x_l` as a paired attribute (because of the underscore) and
`PFDGMPE` needed `rx` for the secondary models.

A fresh in-engine run of the `hazard_curve_minimal` example reproduces
`oq-pfdha`'s `aggregate_hazard.csv` to within **3.75% at the smallest
displacement (0.1 mm), converging to ~6e-4 at 10 m**. The gap is exactly the
`tor`/`rx` divergence of D2: oq-pfdha measures `r` on the original NRML trace
and uses a signed `calculate_signed_site_to_trace_distances()`, while the
engine uses `surface.tor` and `get_rx_distance`.

Proceed with the **second half of PR-6**: `openquake/calculators/displacement.py`
(the `displacement` curve mode storing under `hcurves-rlzs` keyed by
IMT `Disp`). The calculator can now lean on `get_full_lt`/`csm.get_cmakers()`
for source-model realizations and contexts, and on `full_lt.extra_lt` for the
PFD model chains (the `extra_lt` slot is serialized by `FullLogicTree`, so it
is already available in the workers). The end-to-end acceptance is
reproducing `examples/hazard_curve_minimal` (`~/oq-pfdha`, with its INI
migrated to the engine `Disp` format) within the documented `tor` tolerance.

The oq-pfdha FDHA logic-tree files have already been migrated from
`<logicTreeBranchingLevel>` to branch sets directly under `<logicTree>`
(oq-pfdha commit `4453c72`); `PFDLogicTree` accepts both forms.
