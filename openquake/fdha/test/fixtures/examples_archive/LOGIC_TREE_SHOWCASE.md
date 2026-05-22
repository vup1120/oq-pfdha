# FDHA Logic Tree — Slide / Report Pack

Use this file as a **speaker outline**, **figure checklist**, and **runbook** for demonstrating the NRML-based FDHA logic tree (`FdhaLogicTree`), validation demos, and automated tests.

**Repository root** in commands: `/path/to/pfdha` (replace with your clone).

---

## 1. One-minute elevator pitch

- **What:** Epistemic uncertainty in FDHA (primary/secondary surface rupture and displacement models) is expressed as an **OpenQuake-style NRML logic tree** in INI-driven jobs.
- **How:** `[calculation].fdha_logic_tree_file` (singular; deprecated `fdha_logic_tree_files` still accepted in a transition window) routes the CLI to `openquake.fdha.logic_tree.driver.FdhaLogicTree`, which enumerates end-branches, runs the **existing** `BaseFaultRuptureCalculator` / map kernel **unchanged**, then **aggregates annual exceedance rates** (weighted mean + empirical fractiles on rates — not on return periods). NRML from multiple files is merged internally; combined validation uses the merged tree.
- **Why trust it:** Phase-G-style **integration tests** plus **runnable examples** prove single-branch equivalence, closed-form 50/50 blends, style filtering, map-mode arithmetic, and (separately) pixel-wise map equivalence to legacy single-model runs.

---

## 2. Suggested slide deck (≈12 slides)

| # | Title | Key bullets | Suggested figure(s) |
|---|--------|-------------|---------------------|
| 1 | Motivation | FDHA model choice is uncertain; need reproducible, weighted combinations | — |
| 2 | User-facing entry | `fdha job.ini`; `[calculation].fdha_logic_tree_file = tree.xml` (+ `source_model_logic_tree_file`); legacy plural key optional | Snippet from `job_SS_only.ini` or Taiwan INI |
| 3 | NRML concepts | `logicTreeBranchSet`, `uncertaintyModel`, `uncertaintyWeight`; four FDHA `uncertaintyType` values | Tiny XML excerpt |
| 4 | `applyToStyle` | Rake → style (`strike-slip`, `reverse`, `normal`); one INI + merged trees for multi-style sources | Taiwan README table |
| 5 | Aggregation | Mean rate = Σ_branch w·λ; fractiles on **per-branch rates** at fixed D₀; maps from aggregated rates | Formula only |
| 6 | Outputs — curves | `out/hazard_curves/`, `aggregate_hazard.csv`, `manifest.json`, `validator_report.txt` | — |
| 7 | Outputs — maps (v4) | Per-branch `branches/*.h5`; aggregate `rates_mean.h5`, `displacement_map_*.csv` | Taiwan mean map |
| 8 | CLI plotting | `fdha job.ini --plot out.png` — curve and **logic-tree map** quick views (`openquake/fdha/main.py`) | `lt_map.png` or Taiwan PNG |
| 9 | Demo: curve validation | Strike-slip Petersen bilinear vs elliptical; 50/50 matches manual blend | `logic_tree_validation/hazard_curves_comparison.png` |
| 10 | Demo: reverse | Youngs2003 vs Moss2024; same arithmetic check | `logic_tree_validation_reverse/hazard_curves_comparison.png` |
| 11 | Demo: map validation | Map-mode 50/50 vs post-processed blend; displacement checks | `logic_tree_validation/map_mode/hazard_maps_comparison.png` |
| 12 | Demo: Taiwan | 292 end-branches, 3 styles, one INI; multi-source mean hazard map | `logic_tree_validation_taiwan/taiwan_hazard_map_mean.png` |
| 13 | CI / tests | Unit + integration under `openquake/fdha/test/**/logic_tree/` | Pytest command slide |
| 14 | Limitations | No `sourceModel` switching; Moss2013 needs explicit `vs30` in XML when used; metadata-driven advisories (e.g. FDLT-101) | Short bullet list |

Adjust depth: technical audiences can expand slides 3–5; executive audiences can skip to 8–12.

---

## 3. Core features (talking points)

- **Uncertainty slots (1:1 with calculator):**  
  `fdhaPrimarySRModel`, `fdhaPrimaryFDModel`, `fdhaSecondarySRModel`, `fdhaSecondaryFDModel`.
- **Parameter syntax:** `[ClassName]` + INI-style `key = value` inside `<uncertaintyModel>` (reuses existing config parsing).
- **Conditioning:** `applyToBranches` ties dependent branch sets (e.g. primary FD to primary SR).
- **Style filter:** `applyToStyle` limits a branch set to sources classified from rake.
- **Optional source filter:** `applyToSources` still supported where declared source IDs exist.
- **Validators:** FDLT-001 … FDLT-007 (weights, branch IDs, styles, types, registered classes, no placeholder weights); advisory FDLT-101 from `model_metadata.yaml`.
- **Modes:** `hazard_curve` and `hazard_map` both supported under the logic-tree driver.

---

## 4. Public / reference examples (paths relative to repo root)

| Folder | Role | Quick command |
|--------|------|----------------|
| `examples/logic_tree/CharacteristicFaultSourceCase2ClassicalPSHA/` | OQ-style case folder; **runnable** `job_SS_only.ini` (SS tree only); `job.ini` template halts on placeholder NM weights (pedagogical) | `fdha examples/logic_tree/CharacteristicFaultSourceCase2ClassicalPSHA/job_SS_only.ini` |
| `examples/logic_tree_validation/` | Curve: single-branch + 50/50 blend vs analytic mean | `cd examples/logic_tree_validation && python verify.py` |
| `examples/logic_tree_validation_reverse/` | Curve: reverse faulting analogue | `cd examples/logic_tree_validation_reverse && python verify.py` |
| `examples/logic_tree_validation/map_mode/` | Map: same weights as curve demo; `verify_map.py` | `cd examples/logic_tree_validation/map_mode && python verify_map.py` |
| `examples/logic_tree_validation_taiwan/` | Full **multi-style** map demo (292 branches); one INI | `fdha examples/logic_tree_validation_taiwan/fdha_map_taiwan.ini` (long) or `python examples/logic_tree_validation_taiwan/run_taiwan_hazard_map.py` |
| `openquake/fdha/demo/hazard_map/logic_tree_equivalence/` | LT (weight 1.0) **vs** fresh legacy `compute_hazard_map`; 3-panel PNG per scenario | `cd openquake/fdha/demo/hazard_map/logic_tree_equivalence && python compare_lt_vs_legacy.py` |

---

## 5. Committed figures (good for slides)

Copy these into your deck as **evidence** (paths relative to repo root):

| Figure | Description |
|--------|-------------|
| `examples/logic_tree_validation/hazard_curves_comparison.png` | Strike-slip: LT branches + 50/50 mean vs components |
| `examples/logic_tree_validation_reverse/hazard_curves_comparison.png` | Reverse: Youngs2003 vs Moss2024 + blend |
| `examples/logic_tree_validation/map_mode/hazard_maps_comparison.png` | Map-mode blend validation |
| `examples/logic_tree_validation_taiwan/taiwan_hazard_map_mean.png` | Taiwan representative faults — mean displacement map (scripted) |
| `examples/logic_tree_validation_taiwan/lt_map.png` | Same demo — quick **CLI** `--plot` scatter (mean displacement) |

**Regenerated locally (often gitignored):**  
`openquake/fdha/demo/hazard_map/logic_tree_equivalence/*_vs_legacy.png` — run `compare_lt_vs_legacy.py` to produce three-panel legacy vs LT vs difference plots.

---

## 6. Validation story (what to say)

1. **Single-branch equivalence (curve):** Weight 1.0 on each leaf matches the pre-logic-tree calculator on the same virtual INI — `test_single_branch_equivalence.py`.
2. **Two-branch weighted sum:** Analytic `0.5·λ_A + 0.5·λ_B` at every D₀ — `test_two_branch_weighted_sum.py` + `examples/logic_tree_validation/verify.py`.
3. **Style filter:** Rake selects applicable branch sets — `test_style_filter.py`.
4. **Validators halt bad XML:** e.g. weights ≠ 1 — `test_fdlt001_weights_halt.py`; placeholder weights — `test_pedagogical_fdha_ini_halt.py`.
5. **Map mode:** Per-site rate aggregation + inversion checks — `test_map_per_site_arithmetic.py`, `test_map_inversion_consistency.py`.
6. **Map vs legacy:** Pixel-wise match for single-branch trees mirroring shipped map demos — `demo/hazard_map/logic_tree_equivalence/compare_lt_vs_legacy.py` + `summary.json`.
7. **Taiwan:** End-to-end multi-source, multi-style enumeration and **sum of per-source LT means** for total hazard — see `examples/logic_tree_validation_taiwan/README.md`.

---

## 7. Automated tests (copy-paste)

Run from repo root with your Python env activated.

**All logic-tree tests (unit + integration):**

```bash
pytest openquake/fdha/test/unit/logic_tree/ openquake/fdha/test/integration/logic_tree/ -q
```

**Representative integration slices:**

```bash
pytest openquake/fdha/test/integration/logic_tree/test_single_branch_equivalence.py -q
pytest openquake/fdha/test/integration/logic_tree/test_two_branch_weighted_sum.py -q
pytest openquake/fdha/test/integration/logic_tree/test_style_filter.py -q
pytest openquake/fdha/test/integration/logic_tree/test_map_per_site_arithmetic.py -q
```

**Module-level unit tests:**

```bash
pytest openquake/fdha/test/unit/logic_tree/test_validators.py -q
pytest openquake/fdha/test/unit/logic_tree/test_aggregation.py -q
```

---

## 8. CLI one-liners (demo live)

```bash
# Small curve job (logic tree)
fdha examples/logic_tree_validation/single_bilinear/job.ini

# Small map job + save matplotlib figure
fdha examples/logic_tree_validation/map_mode/single_bilinear_map/job.ini --plot /tmp/lt_map_smoke.png

# Verbose Taiwan (long-running)
fdha examples/logic_tree_validation_taiwan/fdha_map_taiwan.ini --plot /tmp/taiwan_lt.png --verbose
```

---

## 9. References inside the repo

- Framework: `openquake/fdha/logic_tree/`
- CLI wiring & LT plotting: `openquake/fdha/main.py`
- Curve validation README: `examples/logic_tree_validation/README.md`
- Map validation README: `examples/logic_tree_validation/map_mode/README.md`
- Reverse demo README: `examples/logic_tree_validation_reverse/README.md`
- Taiwan demo README: `examples/logic_tree_validation_taiwan/README.md`
- Map equivalence README: `openquake/fdha/demo/hazard_map/logic_tree_equivalence/README.md`

---

## 10. Optional “report” PDF workflow

1. Export this Markdown to PDF (Pandoc, VS Code, or GitHub preview print).
2. Embed the PNGs in §5 at full width.
3. Attach `summary.json` from `logic_tree_equivalence/` after a fresh `compare_lt_vs_legacy.py` run as an appendix.

---

*Last updated to match the repository layout as of the logic-tree showcase pack.*
