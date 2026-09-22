# oq-pfdha — Code Overview

A practical tour of what this repository is and how the code is organized.
For the user-facing manual see [`docs/UserManual_Enhanced/`](UserManual_Enhanced/index.md).

## 1. What problem it solves

**Probabilistic Fault Displacement Hazard Analysis (PFDHA)** answers: *at a
given site, what is the annual rate (frequency per year) of earthquake-induced
surface fault displacement exceeding a level `d₀`?*

It is the displacement analogue of standard ground-shaking PSHA. Instead of
ground motion, the "intensity measure" is fault displacement, and it comes from
two places:

- **Principal / primary displacement** — on the fault trace
  (`|r| ≤ r_threshold`): the main fault ruptures the surface at your site.
- **Distributed / secondary displacement** — off-fault
  (`|r| > r_threshold`): splays, shears, and nearby structures rupture.

The result is a **hazard curve** (rate vs. displacement at a point) or a
**hazard map** (displacement at a return period over a grid).

## 2. Inputs and top-level flow

Three inputs feed the `fdha` CLI:

| Input | Role |
|---|---|
| `job.ini` | sites/region, displacement levels, mesh spacing, output options |
| source-model logic tree (NRML XML) | references the OpenQuake fault source model(s) |
| FDHA logic tree (NRML XML) | selects the four scientific models per branch |

Flow (`openquake/fdha/main.py` → `logic_tree/driver.py`):

1. **Parse & validate** the INI, SMLT, FDHA-LT; build the site set
   (`sites` → curve, `region` → map).
2. **Enumerate end branches** — one model chain per combination of logic-tree
   choices, with combined weights.
3. **Run the hazard kernel per branch** (`calc/hazard.py`).
4. **Aggregate** branch rates into a weighted mean + fractiles.
5. **Write outputs** (`out/` next to the INI): `manifest.json`, per-branch
   CSVs/HDF5, aggregate CSVs/HDF5, optional PNG.

## 3. The scientific core

The four model categories form the conditional PFDHA integral
(`docs/UserManual_Enhanced/06-Models.md`):

1. `*PrimarySR` — P(surface rupture on principal fault)
2. `*PrimaryFD` — P(displacement > d │ principal rupture)
3. `*SecondarySR` — P(distributed rupture at distance r)
4. `*SecondaryFD` — P(displacement > d │ distributed rupture)

The kernel (`calc/hazard.py::_compute_rupture_contribution`) evaluates, per
rupture:

```
λ_principal(site,d₀)   = rate · P_sr · P_fd_primary       · W_p(r)
λ_distributed(site,d₀) = rate · P_sr · P_dist_combined(r) · G(r)
λ_total = λ_principal + λ_distributed        # summed over all ruptures
```

`P_dist_combined = P_secondarySR · P_secondaryFD` (or a combined Visini
pipeline).

**Two separate rupture-location-weight paths** (`calc/location_weight.py`), a
key design decision:

- `r_sigma_km == 0` → **boxcar** `W_p = 1{|r| ≤ h}`, with `G = 1 - W_p`
  (complementary split: inside the trace only principal counts, outside only
  distributed).
- `r_sigma_km > 0` → **Petersen's Gaussian** `W_p = exp(-r²/2σ²)`, pinned at 1
  on the trace, truncated at ±2σ, with `G = 1` (principal + distributed
  **summed**).

**Aggregate-definition models** (Kuehn 2024, Lavrentiadis aggregate) already
include the distributed term, so the kernel routes them through a single
bucket: `rate · P_sr · P_fd_aggregate · W_p`, distributed = 0.

The kernel is a **ContextMaker** pattern borrowed from OpenQuake: ruptures are
converted to flat `FDHAContext` arrays (`calc/contexts.py`) of site-rupture
pairs, models receive vectorized inputs, and contributions accumulate by site
id.

## 4. Logic tree and aggregation

`logic_tree/enumerator.py` expands branching levels (model selections +
`fdhaCalcRSigma` scalar branches) with `applyToStyle` / `applyToSources` /
`applyToBranches` filters into weight-carrying `EndBranch`s.
`logic_tree/driver.py` runs each branch, then `logic_tree/aggregation.py`
computes:

- weighted **mean** via `hazardlib.stats.mean_curve`
- weighted **fractiles** via `hazardlib.stats.quantile_curve`

The driver handles nested source-model × FDHA trees, per-source aggregation for
independent faults (means are summed across source groups, not pooled), and
writes a full `manifest.json` describing every realization.

## 5. Calculation modes

- **Hazard curve** (`sites`): direct rate curves per site.
- **Hazard map** (`region` + `region_grid_spacing`): builds a grid site
  collection (active grid + on-trace sites resampled along the fault), computes
  rate cubes, then **inverts** each site's curve at `1/return_period`
  (`calc/utils/interpolation.py::get_map_from_curves`) into a displacement map.

## 6. Module map

| Path | Purpose |
|---|---|
| `main.py` | CLI, dispatch, plotting |
| `calc/config_loader.py`, `params.py` | INI parsing, defaults, validation |
| `calc/calculators.py` | builds the four-model chain, registry, site collection |
| `calc/hazard.py` | the unified hazard kernel |
| `calc/contexts.py` | `FDHAContext` / `FDHAContextMaker`, rake→style, distances |
| `calc/model_adapter.py` | bridges heterogeneous model signatures to the kernel |
| `calc/utils/*` | distances, projection, probability, interpolation, ECS/LCP multi-fault routing |
| `logic_tree/*` | NRML parsing, LT validation, enumeration, branch configs, aggregation, output writers |
| `primary_surf_rup/`, `primary_surf_displ/`, `secondary_surf_rup/`, `secondary_surf_displ/` | the published-model library (e.g. Youngs 2003, Petersen 2011, Takao 2013, Moss 2024, Kuehn 2024, Lavrentiadis 2023, Chiou 2025, Visini 2025…) |
| `demo/` | runnable demos (recent models, r_sigma propagation) |
| `test/` | unit/integration/regression + benchmark reproductions |
| `webgui_demo/` | Streamlit front-end that assembles a job and runs the real engine |

## 7. Design themes worth noting

- **Coherence with OpenQuake**: NRML inputs, hazardlib geometry
  (`OrthographicProjection`), `mean_curve`/`quantile_curve`,
  rates-not-probabilities.
- **Fail loudly**: model errors, site-grid drift across branches, mixed
  displacement definitions, and `r_sigma` scalar-vs-branch conflicts abort the
  run rather than silently producing wrong numbers.
- **Advisory, not silent**: e.g. `ApplicabilityTracker` emits one warning per
  model when sites fall outside a regression's calibrated distance range.
- **Reproducibility**: deterministic output plus frozen fixtures pin the two
  public example jobs byte-for-byte.

## 8. Argument reference: `calc/hazard.py`

`hazard.py` has one public entry point, one per-rupture kernel, and three
helpers. This section documents every argument, with type, default, provenance,
and effect.

### `calculate_fdha_hazard(calculator, sitecol=None, show_progress=True)`

The function the driver calls once per end-branch (and the map path once per
branch/grid).

| Argument | Type / default | What it is and does |
|---|---|---|
| `calculator` | `BaseFaultRuptureCalculator` (required) | The fully configured branch calculator. Everything else is pulled from it: the four `adapters` (model chain), `target_displacements` (the `d₀` grid), `fault_sources` (parsed ruptures), `r_threshold_km`, `r_sigma_km`, the `p_sr_red_cfg`/`s_sr_red_cfg` reductions, its own `sitecol`, and the INI `config` (for `max_distance_km`). Built from a materialized branch INI, so one call = one logic-tree branch. |
| `sitecol` | `Optional[SiteCollection] = None` | Site override. For **hazard curves** it is left `None` and the calculator's own site collection (from `[geometry].sites`) is used. For **hazard maps** the driver passes the pre-built grid `SiteCollection` (active grid + on-trace sites), because the calculator itself does not know the grid. Only the site count and lon/lat are read from it here. |
| `show_progress` | `bool = True` | Wraps the fault-source loop in a `tqdm` bar. The driver turns it off for map runs (many branches) to avoid console spam; the CLI leaves it on. |

**Returns** a dict of numpy arrays: `imls` (`n_displ,`),
`rates`/`rate_principal`/`rate_distributed` (`n_sites, n_displ`, annual
exceedance rates), `site_lons`, `site_lats`, `n_sites`, `n_displ`.

Values derived internally from the calculator (not arguments) matter when
reading the code:

- `max_dist = get_max_distance_km(calculator.config, default=10.0)` — passed
  to `FDHAContextMaker` as the rupture integration cutoff.
- `p_sr_red_cfg` / `s_sr_red_cfg` — how a model's internal MC/epistemic sample
  dimension is collapsed (default `{"method": "mean"}`).
- `r_threshold_km`, `r_sigma_km` — forwarded to the kernel.

### `_compute_rupture_contribution(...)` — the per-rupture kernel

Evaluates one rupture's contribution and returns
`(principal_contrib, distributed_contrib)`, each shape `(N_ctx, n_displ)`.

| Argument | Type / default | Meaning |
|---|---|---|
| `ctx` | `FDHAContext` | The flattened site-rupture context: per-pair `r` (across-strike distance), `rx` (signed hanging-wall/footwall distance), `mag`, `dip`, `rake`, `style`, `occurrence_rate`, `sids`, site coords, plus `metrics_for(...)` for multi-fault reference lines. `N_ctx = len(ctx)`. |
| `adapters` | `Dict[str, LegacyModelAdapter]` | Keyed `'primary_sr'`, `'primary_fd'`, `'secondary_sr'`, `'secondary_fd'`. Each wraps a model instance + its branch parameters and calls its heterogeneous signature. A missing key means "model slot not configured" and the kernel substitutes a neutral value (`P_sr = 1`, or zeros for the others). |
| `target_displacements` | `np.ndarray (n_displ,)` | The `d₀` levels in metres (from `displacement_measure_levels`). Defines the output curve's x-axis. |
| `p_sr_red_cfg` | `Dict[str, Any]` | Reduction config for the **primary** surface-rupture model's internal samples (e.g. `{"method": "mean"}`). Passed straight to `compute_primary_sr`. |
| `s_sr_red_cfg` | `Dict[str, Any]` | Same for the **secondary** models (`compute_secondary_sr`/`compute_secondary_fd`). Defaults to `p_sr_red_cfg`. |
| `r_threshold_km` | `float` | Boxcar half-width `h`, the principal/distributed split. Used **only** on the `sigma == 0` path, where it also sets `G = 1 - W_p`. Ignored when `r_sigma_km > 0`. |
| `use_visini` | `bool` | Switches to the combined Visini SR×FD pipeline instead of the standard `secondary_sr × secondary_fd` product (Visini's regressions are fitted jointly and cannot be split). |
| `visini_calc` | `Optional[VisiniSecondaryCalculator]` | The pre-built Visini calculator (needed only when `use_visini` is true); holds pixel size, along-strike MC settings, and rank-1.5 traces. |
| `calculator` | `BaseFaultRuptureCalculator` | Used in the Visini branch to read style/case/parameters (`get_model_parameters`, `case_label`) and the secondary model's declared `MULTIFAULT_REFERENCE_LINE`. |
| `r_sigma_km` | `float = 0.0` | Two-sided mapping-accuracy σ. `0` → boxcar path (`W_p`, `G = 1-W_p`, complementary). `> 0` → Petersen Gaussian (`W_p`, `G = 1`, additive). Default 0 preserves legacy behaviour. |

What it does with them:

```
W_p = location_weight(ctx.r, r_threshold_km, r_sigma_km)
G   = 1 - W_p   if r_sigma_km == 0 else 1
principal   = rate · P_sr · P_fd_primary      · W_p
distributed = rate · P_sr · P_dist_combined   · G
```

with `rate = ctx.occurrence_rate[0]`. There is one early return: if the primary
FD model declares `DISPLACEMENT_DEFINITION == "aggregate"`, the distributed
bucket is zero and everything flows through `principal` (single-bucket path).

### `ApplicabilityTracker.__init__(r_threshold_km, r_sigma_km)`

- `r_threshold_km` — boxcar half-width; used in `observe` to decide whether a
  site's distributed term actually contributes on the σ=0 path (inside the
  boxcar the complementary split masks it).
- `r_sigma_km` — if `0`, the same masking rule applies; if `> 0`, the
  distributed term is additive everywhere, so every out-of-range site counts.

`observe(model, ctx)` takes the distributed FD model (to read its
`APPLICABILITY_RANGE` and `MULTIFAULT_REFERENCE_LINE`) and the context, and
records offending site ids; `emit()` logs one warning per model. Warning only —
no behaviour change.

### `_resolve_investigation_time(src, ini_time)`

- `src` — a hazardlib source object; its NRML-header `investigation_time` (if
  any) is authoritative for converting a non-parametric source's `probs_occur`
  into annual rates.
- `ini_time` — the `[calculation].investigation_time` fallback. If both exist
  and disagree, it raises (silently rescaling every non-parametric rate is the
  failure it prevents). Parametric sources ignore the value. If neither is
  present, it returns `1.0`.

### `_setup_visini_calculator(calculator)`

Takes only the branch `calculator` and returns `(use_visini, visini_calc)`.
It reads `SECONDARY_PIPELINE` from the configured secondary SR/FD models; if
either declares `'visini'`, it attaches any `rank1p5` XML surfaces and builds a
`VisiniSecondaryCalculator` from the models' parameters (`pixel_size`,
`along_strike_width`, `near_far_threshold_km`, `rupture_traces`,
`rank1p5_traces`, `segment_sampling`, `distribution_type`, `case_label`).
Otherwise it returns `(False, None)` and the kernel uses the standard
`SR × FD` path.

`get_fdha_params()` on the calculator supplies `r_threshold_km`,
`near_far_threshold_km`, `multifault_reference_lines`, `r_sigma_km`,
`surface_rupture_depth_tolerance_km`, and `reference_vs30_value` to the context
maker — none of these are `hazard.py` arguments; they reach the computation
through the calculator.
