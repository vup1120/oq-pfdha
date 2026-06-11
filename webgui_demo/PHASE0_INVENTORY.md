# Phase 0 — Read-Only Discovery Inventory

Repository: `oq-pfdha` (working copy at `/home/user/oq-pfdha`).
All paths below are relative to the repository root. Line numbers refer to the
current checkout (branch `claude/inspiring-tesla-rmtxfq`, derived from `main`).

---

## 1. Model registry

**Registry mechanism.** There is no separate registry file. The effective
registry is the public namespace of the four model packages: the logic-tree
validator accepts a model class iff `hasattr(<package>, class_name)` for one
of the four packages — see `_class_is_registered()` in
`openquake/fdha/logic_tree/validators.py:184-195`. The mapping from
logic-tree slot to package is `FDHA_SLOTS_BY_UTYPE` in
`openquake/fdha/logic_tree/types.py:14-19`.

The GUI must therefore populate its four dropdowns from exactly the classes
exported by each package `__init__.py`, listed below.

### 1.1 Principal (primary) surface-rupture probability — `openquake.fdha.primary_surf_rup`
Exports: `openquake/fdha/primary_surf_rup/__init__.py:24-33`.

| Class name | Defined at |
|---|---|
| `Youngs2003PrimarySR` | `openquake/fdha/primary_surf_rup/youngs2003.py:27` |
| `WC1993PrimarySR` | `openquake/fdha/primary_surf_rup/wells_coppersmith1993.py:27` |
| `MossRoss2011PrimarySR` | `openquake/fdha/primary_surf_rup/moss_ross2011.py:28` |
| `Takao2013PrimarySR` | `openquake/fdha/primary_surf_rup/takao2013.py:29` |
| `Moss2013PrimarySR` | `openquake/fdha/primary_surf_rup/moss2013.py:27` |
| `Yang2021PrimarySR` | `openquake/fdha/primary_surf_rup/yang2021.py:34` |
| `Pizza2023PrimarySR` | `openquake/fdha/primary_surf_rup/pizza2023.py:27` |
| `Mammarella2024PrimarySR` | `openquake/fdha/primary_surf_rup/mammarella2024.py:177` |
| `MammarellaEtAl2024PrimarySR` (alias subclass) | `openquake/fdha/primary_surf_rup/mammarella2024.py:395` |
| `FixedPrimarySR` (utility, constant value) | `openquake/fdha/primary_surf_rup/fixed.py:27` |

### 1.2 Principal (primary) fault displacement — `openquake.fdha.primary_surf_displ`
Exports: `openquake/fdha/primary_surf_displ/__init__.py:23-35`.

| Class name | Defined at |
|---|---|
| `Youngs2003PrimaryFD` | `openquake/fdha/primary_surf_displ/youngs2003.py:15` |
| `Petersen2011PrimaryFD` | `openquake/fdha/primary_surf_displ/petersen2011.py:34` |
| `Petersen2011PrimaryFD_bilinear` | `openquake/fdha/primary_surf_displ/petersen2011.py:190` |
| `Petersen2011PrimaryFD_elliptical` | `openquake/fdha/primary_surf_displ/petersen2011.py:202` |
| `Petersen2011PrimaryFD_quadratic` | `openquake/fdha/primary_surf_displ/petersen2011.py:214` |
| `MossRoss2011PrimaryFD` | `openquake/fdha/primary_surf_displ/moss_ross2011.py:24` |
| `Moss2022PrimaryFD` | `openquake/fdha/primary_surf_displ/moss2022.py:47` |
| `Moss2024PrimaryFD` | `openquake/fdha/primary_surf_displ/moss2024.py:14` |
| `Takao2013PrimaryFD` | `openquake/fdha/primary_surf_displ/takao2013.py:31` |
| `Lavrentiadis2023PrimaryFD` | `openquake/fdha/primary_surf_displ/lavrentiadis2023.py:28` |
| `Kuehn2024PrimaryFD` | `openquake/fdha/primary_surf_displ/kuehn2024/kuehn2024.py:39` |
| `Chiou2025PrimaryFD` | `openquake/fdha/primary_surf_displ/chiou2025.py:47` |

### 1.3 Distributed (secondary) surface-rupture probability — `openquake.fdha.secondary_surf_rup`
Exports: `openquake/fdha/secondary_surf_rup/__init__.py:24-33`.

| Class name | Defined at |
|---|---|
| `Youngs2003SecondarySR` | `openquake/fdha/secondary_surf_rup/youngs2003.py:27` |
| `Petersen2011SecondarySR` | `openquake/fdha/secondary_surf_rup/petersen2011.py:27` |
| `Petersen2011SecondarySR_default` | `openquake/fdha/secondary_surf_rup/petersen2011.py:152` |
| `Visini2025SecondarySR` | `openquake/fdha/secondary_surf_rup/visini2025.py:51` |
| `Takao2014SecondarySR` | `openquake/fdha/secondary_surf_rup/takao2014.py:40` |
| `Takao2013SecondarySR` | `openquake/fdha/secondary_surf_rup/takao2013.py:36` |
| `FerrarioLivio2021SecondarySR` | `openquake/fdha/secondary_surf_rup/ferrario2021.py:35` |
| `Rodriguez2023SecondarySR` | `openquake/fdha/secondary_surf_rup/rodriguez2023.py:35` |
| `Moss2022SecondarySR` | `openquake/fdha/secondary_surf_rup/moss2022.py:51` |
| `FixedSecondarySR` (utility, constant value) | `openquake/fdha/secondary_surf_rup/fixed.py:27` |

### 1.4 Distributed (secondary) fault displacement — `openquake.fdha.secondary_surf_displ`
Exports: `openquake/fdha/secondary_surf_displ/__init__.py:23-26`.

| Class name | Defined at |
|---|---|
| `Youngs2003SecondaryFD` | `openquake/fdha/secondary_surf_displ/youngs2003.py:29` |
| `Petersen2011SecondaryFD` | `openquake/fdha/secondary_surf_displ/petersen2011.py:29` |
| `Visini2025SecondaryFD` | `openquake/fdha/secondary_surf_displ/visini2025.py:44` |
| `Moss2022SecondaryFD` | `openquake/fdha/secondary_surf_displ/moss2022.py:37` |

**GUI sourcing recommendation:** import the four packages at app start and
enumerate exported classes (mirrors `_class_is_registered`), rather than
hard-coding the tables above.

---

## 2. Norcia Case 3 benchmark (IAEA TECDOC-2092)

Directory: `openquake/fdha/test/benchmark/norcia_case3_iaea/`
Reference: IAEA-TECDOC-2092 (2025), DOI 10.61092/iaea.74us-dn4n — see
`REFERENCE.md` in the directory.

### 2.1 Input files (tracked in git)

| File | Role |
|---|---|
| `job_norcia_case3_iaea_curve.ini` | Hazard-curve job, 3 sites (PF `13.278 42.767`, MS `13.188 42.749`, SL `13.212 42.853`) |
| `job_norcia_case3_iaea_map.ini` | Hazard-map job (region `13.05–13.35 E, 42.65–42.95 N`, grid 0.005°, `return_period = 100000`) |
| `job_norcia_case3_iaea_ms_curve.ini` | Single-site (MS) curve variant |
| `source_model_norcia_case3.xml` | NRML source model, sources `MVFS` and `NFS` (normal faulting, rake −90) |
| `config_norcia_case3_iaea_source_model_logic_tree.xml` | Source-model logic tree (1 branch, weight 1.0) |
| `config_norcia_case3_iaea_fdha_logic_tree.xml` | FDHA logic tree: 2 end-branches (Y03 chain: Youngs2003 throughout; V24 chain: Youngs2003 principal + `Visini2025SecondarySR`/`Visini2025SecondaryFD` distributed), weights 0.5/0.5 at the first branching level |

### 2.2 Output files — ⚠️ status caveat

The output directories `out_curve/` and `out_map/` are **gitignored**
(`.gitignore:24` `**/out_curve*/`; outputs in general `.gitignore:22`
`**/out/`). They exist in this working copy only because the benchmark
runner `run_test.py` was executed in this session. A fresh clone does NOT
contain them. **Decision needed at the Phase 0 gate** (see §5, Open
question A): I propose copying the needed CSVs into `webgui_demo/data/`
so the demo is self-contained, which also respects the "no modifications
outside `webgui_demo/`" rule.

Files currently present (produced by `run_test.py`):

**Curve job → `out_curve/`**
- `aggregate_hazard.csv` — columns: `site_id, lon, lat, D0, mean, quantile-0.05, quantile-0.16, quantile-0.5, quantile-0.84, quantile-0.95`; 3 sites × 20 displacement levels. **This is the file for the Results page hazard-curve plot (mean + fractiles).**
- `hazard_curves/branch_0000.csv`, `branch_0001.csv` — per-branch curves, columns `site_id, lon, lat, D0, annual_rate`.
- `manifest.json` — branch metadata and weights.
- `source_model_branches/.../validator_report.txt`, `branch_configs/branch_000N.ini`.

**Map job → `out_map/`**
- `aggregate/displacement_map_mean.csv` — header comment `# return_period = 100000.0`; columns `site_id, lon, lat, is_trace, displ_mean`. **Map outputs DO exist → Page 3 can render a real map.**
- `aggregate/displacement_map_quantile-{0.05,0.16,0.5,0.84,0.95}.csv` — same schema, fractile maps.
- `aggregate/rates_mean.h5`, `aggregate/rates_fractiles.h5` — HDF5 rate grids.
- `component_maps.npz`, `manifest.json`, per-source-model-branch subtree with per-branch `branches/branch_000N.h5`.

Note: rates in the curve CSVs are annual exceedance **rates** (1/yr), not
probabilities (column name `annual_rate`; aggregate column `mean`). Axis
labels on the Results page should say "annual rate of exceedance (1/yr)" to
match the codebase's own plot labels (`openquake/fdha/calc/io_utils.py`).

---

## 3. INI parameters (as parsed by the codebase)

Canonical loader: `load_config()` / `_load_ini_config()` in
`openquake/fdha/calc/config_loader.py:91-139` (v5 canonical INI only; TOML
explicitly rejected at `config_loader.py:116-121`). Section/key
normalisation (OpenQuake-style aliases) in `_normalize_config`,
`config_loader.py:262-360`. Calculation type is auto-detected:
`hazard_map` iff `[geometry]` contains `region`, else `hazard_curve`
(`_detect_calculation_type`, `openquake/fdha/logic_tree/driver.py:1057-1063`).

| Section / key | Default (traceable) | Source of default | Notes |
|---|---|---|---|
| `[general] description` | — (free text) | — | informational |
| `[general] calculation_mode` | `fdha_classical` by convention | used in all shipped jobs (e.g. `job_norcia_case3_iaea_curve.ini:3`) | **not read programmatically** — no parser reference found; keep fixed/disabled in GUI |
| `[geometry] sites` | — required for curves | e.g. `job_norcia_case3_iaea_curve.ini:6` | `lon lat, lon lat, ...` |
| `[geometry] region` | — required for maps | `job_norcia_case3_iaea_map.ini:6` | corner list; presence switches mode (driver.py:1058) |
| `[geometry] region_grid_spacing` | — required for maps | `job_norcia_case3_iaea_map.ini:7` | degrees |
| `[geometry] max_distance_km` | — optional | alias from `[calculation]` handled at `config_loader.py:323-328` | site filtering |
| `[site_params] reference_vs30_value` | — optional | mapped to `site_location.vs30` at `config_loader.py:307-311` | required by `Moss2013PrimarySR` (Vs30-dependent) |
| `[erf] rupture_mesh_spacing` | `0.5` | `ERFConfig`, `config_loader.py:35` | km |
| `[erf] width_of_mfd_bin` | `0.1` | `ERFConfig`, `config_loader.py:36` | magnitude units |
| `[calculation] investigation_time` | `1.0` | `CalculationConfig`, `config_loader.py:50` | years |
| `[calculation] displacement_measure_levels` | — required | JSON parsing noted `config_loader.py:133`; normalised to `parameters.target_displacement` at `config_loader.py:313-318` | JSON: `{"FD": [...]}` |
| `[calculation] r_threshold_km` | `0.1` | `openquake/fdha/calc/calculators.py:115-117`; map path `openquake/fdha/calc/hazard_map.py:38` | principal vs distributed split |
| `[calculation] near_far_threshold_km` | `0.2` | `calculators.py:123-124`; normalisation `config_loader.py:355-360` | Visini near/far regime |
| `[calculation] return_period` | `100000` | `hazard_map.py:37`; logic-tree path `driver.py:251-252` | **hazard maps only** |
| `[calculation] source_model_logic_tree_file` | — required | read at `config_loader.py:225` | NRML XML |
| `[calculation] fdha_logic_tree_file` | — required (singular) | read at `config_loader.py:214,256` | single NRML XML file |
| `[output] mean` | `True` | `config_loader.py:643-652` | emit weighted-mean |
| `[output] quantiles` | `(0.05, 0.16, 0.5, 0.84, 0.95)` | `FRACTILE_QS`, `openquake/fdha/logic_tree/io.py:200` | space-separated floats |

### 3.1 ⚠️ Discrepancies vs. the task brief

1. **`truncation_level` is NOT an INI parameter anywhere in the codebase**
   (repo-wide search for `truncation` matches only per-model internals and
   docs). Distribution truncation is a **per-model parameter `n_sigma`**
   passed inside the logic-tree `<uncertaintyModel>` block, e.g.
   `Youngs2003PrimaryFD` default `n_sigma = 6`
   (`openquake/fdha/primary_surf_displ/youngs2003.py:53-58`),
   `Takao2013PrimaryFD` default `n_sigma = 3`
   (`openquake/fdha/primary_surf_displ/takao2013.py:43-44`).
   → Per constraint #3 ("INI names must match the codebase exactly") the GUI
   will **omit** a `truncation_level` field. See Open question B.
2. **`principal_distance_km`** — confirmed deprecated/inert: it appears only
   inside example/test INI text (e.g.
   `openquake/fdha/test/integration/logic_tree/test_source_model_logic_tree_map.py:210`)
   and is never read by any calculator. Excluded from the GUI, as instructed.

---

## 4. NRML logic-tree schema (as consumed by the parser)

Parser: `openquake/fdha/logic_tree/nrml_reader.py:24-56`; data model:
`openquake/fdha/logic_tree/types.py`.

- Namespace: `http://openquake.org/xmlns/nrml/0.4`.
- Element tree: `<nrml>` → `<logicTree logicTreeID>` →
  `<logicTreeBranchingLevel branchingLevelID>` →
  `<logicTreeBranchSet branchSetID uncertaintyType [applyToSources]
  [applyToBranches] [applyToStyle]>` →
  `<logicTreeBranch branchID>` → `<uncertaintyModel>` (CDATA:
  `[ClassName]` + optional newline-separated `key = value` parameters) +
  `<uncertaintyWeight>` (attribute names read at `nrml_reader.py:33-49`).
- Valid `uncertaintyType` values (`types.py:7-12`):
  `fdhaPrimarySRModel`, `fdhaPrimaryFDModel`, `fdhaSecondarySRModel`,
  `fdhaSecondaryFDModel` → slots `primary_surf_rup`, `primary_surf_displ`,
  `secondary_surf_rup`, `secondary_surf_displ` (`types.py:14-19`).
- `applyToStyle` values (`ALLOWED_STYLES`, `types.py:21`):
  `strike-slip`, `reverse`, `normal`.
- Validation rules the GUI must mirror (`validators.py`):
  - **FDLT-001** (`validators.py:133-143`): weights **per branch set** must
    sum to 1.0 ± 1e-6 — exactly the constraint required by the task.
  - **FDLT-002** (`validators.py:166-179`): `applyToBranches` must reference
    declared branch IDs.
  - Unregistered class error (`validators.py:116-121` + `:184-195`).
  - **FDLT-101** advisory (`validators.py:145-164`): mixed displacement
    definitions (AD vs MD) within one Primary-FD branch set.
- **Single file**: the job key is `fdha_logic_tree_file` (singular,
  `config_loader.py:214`); one XML can hold all four branching levels —
  the Norcia tree (`config_norcia_case3_iaea_fdha_logic_tree.xml`) is the
  canonical template, including `applyToBranches` chaining and CDATA
  parameter blocks.

A concrete serialization template for the GUI is the Norcia FDHA logic tree
itself (2 chains across 4 branching levels, weights 0.5/0.5 at level 1 and
1.0 within each dependent branch set).

---

## 5. Open questions for the Phase 0 gate

**A. Pre-computed outputs are gitignored.** The task assumes Norcia Case 3
outputs are "already present in the repository", but `out_curve/`/`out_map/`
are ignored by `.gitignore` and exist only in this working copy (generated
this session via the repo's own `run_test.py`). Proposal: copy
`aggregate_hazard.csv`, the per-branch curve CSVs, `manifest.json`,
`displacement_map_mean.csv` and the five fractile-map CSVs into
`webgui_demo/data/norcia_case3/` (≈ a few hundred kB, all CSV/JSON) and have
the app read only from there. Confirm or choose an alternative (e.g.
documented "run `run_test.py` first" prerequisite).

**B. `truncation_level` does not exist in the codebase** (§3.1). Proposal:
expose instead the real per-model `n_sigma` parameter inside the logic-tree
builder for models that accept it (defaults traceable per model), or omit
truncation entirely from the demo form. Confirm preference.

**C. `Fixed*` and alias classes in dropdowns.** `FixedPrimarySR` /
`FixedSecondarySR` are utility constants, `MammarellaEtAl2024PrimarySR` and
`Petersen2011SecondarySR_default` are aliases. Include them verbatim (registry
fidelity) or filter utilities/aliases out for clarity? Default if no
preference: include everything the registry exports, since constraint #1
says "populated ONLY from the actual model registry" — filtering is a
curation decision I'd like confirmed.
