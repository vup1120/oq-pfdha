# DISCOVERY.md — r_threshold_km epistemic-uncertainty task, Phase 0

Read-only discovery for the "r_threshold_km as logic-tree branches" task.
No code was modified. All line numbers refer to the current working tree
(branch `percentile-displacement`, commit f7070be3 + untracked LCP/segments files).

---

## 1. `r_threshold_km` end-to-end trace

### 1.1 INI parsing

- `openquake/fdha/calc/config_loader.py:91` — `load_config()` is the single
  canonical INI loader. It performs **no special handling** of
  `r_threshold_km`: the key is parsed generically by `_parse_ini_value()`
  (`config_loader.py:548`) into `config['calculation']['r_threshold_km']`
  (float). Note the contrast with `near_far_threshold_km`, which *is*
  explicitly copied from `[calculation]` to `[parameters]` during
  normalisation (`config_loader.py:355-360`). `r_threshold_km` has no such
  copy step; consumers read both sections themselves.

### 1.2 Existing implementation default

The default **0.1 km** is hardcoded at each consumption point (there is no
single source of truth):

| File:line | Expression |
|---|---|
| `openquake/fdha/calc/calculators.py:115-118` | `float(config['calculation'].get('r_threshold_km') or config['parameters'].get('r_threshold_km', 0.1))` |
| `openquake/fdha/calc/hazard_map.py:38` | `float(calc.get('r_threshold_km', para.get('r_threshold_km', 0.1)))` |
| `openquake/fdha/calc/contexts.py:283` | `fdha_params.get('r_threshold_km', 0.1)` |

**Quirk (pre-existing, must NOT be changed per HARD RULES):**
`calculators.py:116` uses `or`, so an explicit `r_threshold_km = 0` (falsy)
in `[calculation]` silently falls through to `[parameters]` / the 0.1
default. MODE A bit-identity requires leaving this as-is.

### 1.3 Runtime config field

- `BaseFaultRuptureCalculator._initialize_calculation_params()` →
  `self.r_threshold_km` (`calculators.py:115`). This is **the** runtime field.
- `BaseFaultRuptureCalculator.get_fdha_params()` (`calculators.py:158-173`)
  re-exports it as `{'r_threshold_km': ..., 'near_far_threshold_km': ...}`
  to the context maker.

### 1.4 Consumption points

**Hazard-curve AND logic-tree hazard-map path (the live path):**

1. `calc/hazard.py:103` — `calculate_fdha_hazard()` reads
   `calculator.r_threshold_km`.
2. `calc/hazard.py:141` — passed into `_compute_rupture_contribution()`.
3. `calc/hazard.py:396-397` — **the actual routing decision**:
   `mask_principal = np.abs(ctx.r) <= r_threshold_km`;
   `mask_distributed = ~mask_principal`. Principal contribution at
   `hazard.py:401-403`, distributed at `hazard.py:406-416`.

Both logic-tree modes converge here: curve mode via
`driver._run_single()` → `FaultRuptureProbabilityCalculator` →
`calc.run()` → `calculate_fdha_hazard`; map mode via
`driver._run_map()` (`driver.py:290-335`) which builds
`BaseFaultRuptureCalculator(tmp_ini, ...)` and calls
`calculate_fdha_hazard` directly. **One consumption point:
`calculators.py:115` reading the (branch) INI, feeding `hazard.py:396`.**

**Secondary / dormant paths (verified not on the logic-tree route):**

- `calc/hazard_map.py:38` — legacy single-run map entry
  `compute_hazard_map()` reads it into `principal_dist`, which appears
  exactly once in the file (the assignment) — a **dead read**; that path
  routes through `calculate_fdha_hazard` anyway. Reached from
  `calc/run.py:41` (legacy CLI) and `webgui_demo/app.py`, not from
  `FdhaLogicTree`.
- `calc/contexts.py:283` — `FDHAContextMaker` stores `self.r_threshold_km`;
  `FDHAContext.get_principal_mask()/get_distributed_mask()`
  (`contexts.py:203-225`) exist but have **no production callers** (grep:
  only contexts.py itself and tests). Informational only.

### 1.5 Branch-INI flow (how the value reaches per-branch runs today)

`FdhaLogicTree._write_branch_ini()` (`driver.py:1014-1033`) →
`config_builder.build_config()` (`config_builder.py:12-21`) does
`copy.deepcopy(base_config)`, so the job-level
`[calculation].r_threshold_km` scalar is copied verbatim into every
materialised `branch_configs/branch_XXXX.ini`, which
`calculators.py:115` then reads. This is exactly the seam MODE B needs:
overwrite `cfg['calculation']['r_threshold_km']` per end-branch and the
runtime consumes it with **zero new routing logic** (HARD RULE 3-2
satisfied by construction).

Materialised branch INIs bypass public-INI legacy rejection via
`_is_materialized_branch_ini()` (`config_loader.py:166-168`, keyed on the
`branch_configs` path component) — important: per-branch scalar injection
will not trip public-INI validation, and the CONFLICT RULE check must run
on the **job** INI only, never on materialised branch INIs.

---

## 2. FDHA logic-tree machinery map

All under `openquake/fdha/logic_tree/`:

| Concern | Location |
|---|---|
| uncertaintyType registry | `types.py:7-12` `FDHA_UNCERTAINTY_TYPES` (4 model types); `types.py:14-19` `FDHA_SLOTS_BY_UTYPE` maps type → model slot name |
| XML → spec parsing | `nrml_reader.py:24-77` `parse()` — generic; any `uncertaintyType` string is accepted at parse time and carried on `BranchSet.uncertainty_type` |
| `<uncertaintyModel>` text parsing | `param_parser.py:8-28` `parse_uncertainty_model()` — supports plain class name or `[ClassName]\nkey = val` INI block. A line like `r_threshold_km = 0.2` would currently be (mis)read as a *class name*; a type-specific parse path is needed |
| Branch-set validation | `validators.py:44-181` `validate_spec()` — **FDLT-005** rejects any type not in `FDHA_UNCERTAINTY_TYPES`; **FDLT-001** weight sum with tolerance `1e-6` (`validators.py:136`) — this is the validator tolerance to reuse; FDLT-006 checks class registration (must be skipped for a parameter-type branch set); FDLT-007 placeholder weights |
| End-branch enumeration | `enumerator.py:26-71` `enumerate_end_branches()` — Cartesian product per source across branching levels; weights multiply (`enumerator.py:59`). Branch sets whose type has no slot are **silently passed through** (`enumerator.py:39-42`) — i.e. today an unknown type is dropped by the enumerator but blocked earlier by FDLT-005 |
| End-branch shape | `types.py:62-67` `EndBranch(source_id, style, weight, selections: dict[slot, ModelChoice])`; `ModelChoice` (`types.py:54-59`) carries `class_name, params, branch_id, weight` |
| Dedup | `driver.py:1215-1251` `_dedup_end_branches()` — fingerprint over `selections` (`driver.py:1198-1207`, includes params) so distinct threshold values will *not* collapse; identical ones across same-style sources will, correctly |
| Branch-INI materialisation | `driver.py:1014-1033` `_write_branch_ini()` + `config_builder.py:12-21` `build_config()` (writes selections into `cfg['models'][slot]`) + `config_builder.py:24-55` `dump_config_to_ini()` |
| Per-branch param flow XML→INI→runtime | **Exists only for model slots**: `ModelChoice.params` → `[models.<slot>.parameters]` in branch INI → `LegacyModelAdapter`. There is **no existing mechanism** for a branch to set a `[calculation]` key; `build_config()` must be extended for MODE B |
| Manifest | `driver.py:1254-1291` `_build_manifest()` (single-SMLT map/curve), `driver.py:952-1011` `_build_multi_source_manifest()`, map-multi at `driver.py:867-908`. Branch IDs surface via `_fdha_branch_path()` (`driver.py:1186-1195`, joins `selections[slot].branch_id` over sorted slots) and `fingerprint` |
| Aggregation | `aggregation.py:6-14` `weighted_mean` (arithmetic mean of rates), `aggregation.py:17-65` `weighted_fractiles` (weighted empirical quantiles, linear-in-CDF). Multi-source grouping in `driver.py:1354-1419`. **Untouched by this task** — threshold branches are ordinary realizations |
| Job orchestration | `driver.py:62-133` `FdhaLogicTree.__init__/run` (INI-level checks live here and in `config_loader.validate_public_logic_tree_ini_file`, `config_loader.py:231-266`); curve: `_run_curve_or_map_multi_source` (`driver.py:488`); map: `_run_map_multi_source` (`driver.py:648`) → `_run_map` (`driver.py:217`) |
| CLI entry | `main.py:123-159` `_run_logic_tree()` → `FdhaLogicTree.from_ini().run()` |
| Model metadata | `logic_tree/model_metadata.yaml` — **only** `primary_fd_models` with `displacement_definition` + citation. **No distance-applicability metadata exists anywhere** (grepped `secondary_surf_displ/`, `secondary_surf_rup/` for applicability/valid_range/DEFINED_FOR: nothing). ⇒ Phase 2's cross-check must take the single "applicability-range metadata unavailable" WARNING path |

Weight-combination convention already in place (MODE B slots straight in):
combined weight = `sm_branch.weight * eb.weight` (`driver.py:564`,
`driver.py:763-765`); normalised at `driver.py:609-614`.

---

## 3. Existing jobs & tests exercising `r_threshold_km`

**Example jobs (public):**
- `examples/hazard_curve_minimal.ini:15` — sets `r_threshold_km = 0.1`
  explicitly in `[calculation]`; canonical curve job (SMLT + FDHA LT XMLs
  alongside).
- `examples/hazard_map_minimal.ini` — map job; does **not** set
  `r_threshold_km` (exercises the implementation default). Good pair for
  the Phase 1 "explicit + omitted" requirement.

**Tests touching the key:**
- `openquake/fdha/test/unit/test_config_loader.py:107` — INI containing the
  key parses.
- `openquake/fdha/test/unit/test_multi_site_calculator.py:120` — calculator
  built from INI with `r_threshold_km = 0.1`.
- `openquake/fdha/test/unit/test_multifault_source.py:235` — uses
  `r_threshold_km = 1.0`.

**Logic-tree test suites (patterns to mirror in Phases 2–3):**
- Unit: `test/unit/logic_tree/` — `test_validators.py`,
  `test_enumerator.py`, `test_param_parser.py`, `test_config_builder.py`,
  `test_nrml_reader.py`, `test_aggregation.py`.
- Integration: `test/integration/logic_tree/` — notably
  `test_single_branch_equivalence.py` (single-branch ≡ direct run — the
  template for the Phase 3a EQUIVALENCE test),
  `test_two_branch_weighted_sum.py` (template for Phase 3b ANALYTIC
  AGGREGATION), `test_manifest_completeness.py` (Phase 3d),
  `test_fdlt001_weights_halt.py` (weight-sum halt pattern).

**Usercase/benchmark jobs** set the scalar too (read-only context):
`test/usercase/hugo/hazard_map_calabria*/job.ini` (0.01),
`test/usercase/gao/*`, `test/benchmark/norcia_*`, `test/benchmark/visini_*`.
Their committed `out/.../branch_configs/*.ini` files confirm the deepcopy
flow described in §1.5.

---

## 4. Flags: where the task design meets the actual code

1. **`EndBranch.selections` is model-slot-shaped.** The natural
   implementation is a pseudo-slot (e.g. `calc_r_threshold`) whose
   `ModelChoice` carries `params={'r_threshold_km': V}` and a blank/marker
   `class_name`, with `build_config()` special-casing it into
   `cfg['calculation']['r_threshold_km']` instead of `cfg['models']`.
   This automatically gives: Cartesian combination + weight multiplication
   (enumerator), non-collapsing dedup fingerprints, and branch-ID
   visibility in `manifest.json` via `_fdha_branch_path`. The alternative —
   a new parallel enumeration axis outside `selections` — would touch far
   more of `driver.py`. To be decided at the Phase 2 schema gate.
2. **Validator assumes every branch is a model class.** FDLT-006
   (`validators.py:114-123`) would flag `r_threshold_km = 0.2` as an
   unregistered class; validation must branch on `uncertainty_type` before
   the class check. FDLT-005 registry (`types.py:7`) must gain
   `fdhaCalcRThreshold`.
3. **Enumerator silently ignores unknown types** (`enumerator.py:39-42`).
   Once FDLT-005 admits the new type, forgetting to give it a slot would
   *silently drop* the branch set — the Phase 2/3 tests must pin this.
4. **CONFLICT RULE placement.** The scalar lives in
   `self.base_config['calculation']` (`driver.py:66`); the branch set is
   known after `parse_nrml` inside `_enumerate_and_validate()`
   (`driver.py:146-151`). The check must (a) run once per job on the *job*
   INI (never on `branch_configs/` INIs, §1.5), and (b) fire before any
   branch runs. Note `_enumerate_and_validate` is re-run per SMLT branch
   (`driver.py:542`, `driver.py:704`) — the conflict error must not be
   duplicated/masked by that loop.
5. **The `or`-based falsy default** (`calculators.py:116`, §1.2) means a
   MODE B branch value materialised as `0` (positive-float validation in
   Phase 2 should prevent this anyway) would silently become 0.1. Phase 2's
   "must parse as positive float" rule closes this hole at the XML layer;
   the runtime line itself stays untouched.
6. **No distance-applicability metadata exists** (§2, last row) ⇒ Phase 2
   must implement the single "applicability-range metadata unavailable"
   WARNING, not per-model range checks. `model_metadata.yaml` covers only
   primary-FD displacement definitions.
7. **Fixture path.** Task says `tests/fixtures/r_threshold_baseline/`; the
   repo convention is `openquake/fdha/test/fixtures/` (existing:
   `chelungpu_mesh/`, `ecs/`, `lcp/`, `examples_archive/`). Proposal:
   `openquake/fdha/test/fixtures/r_threshold_baseline/` — confirm at the
   Phase 1 gate.
8. **No CHANGELOG file exists** in the repo (only `README.md`; error
   messages cite "CHANGELOG section vX.Y" aspirationally,
   `config_loader.py:181-201`). Phase 4 "CHANGELOG entry" needs a decision:
   create `CHANGELOG.md`, or note the change elsewhere — confirm at the
   Phase 4 gate.
9. **Weight-sum tolerance to reuse:** `1e-6` from FDLT-001
   (`validators.py:136`). No new tolerance will be defined.
10. **Two default-bearing side paths** (`hazard_map.py:38` legacy map entry;
    `contexts.py:283` context maker) are off the logic-tree route (§1.4).
    MODE B needs to touch neither; MODE A leaves both bit-identical.

---

*Phase 0 complete. No edits made outside this report. Awaiting confirmation
before Phase 1 (regression baseline).*
