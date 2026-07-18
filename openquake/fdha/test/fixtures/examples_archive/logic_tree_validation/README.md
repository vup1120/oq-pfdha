# Logic Tree Validation Demo

Self-contained demo that verifies `FdhaLogicTree` aggregates per-branch hazard
results correctly. Matches the design memo (`Logic Tree Validation Demo -
Design v1`) and is the runnable companion to Phase G integration tests 1 and 2
of the v3 FDHA Logic Tree Integration prompt.

What it proves, by construction:

1. **Regression.** Two single-branch trees (weight = 1.0) on
   `Petersen2011PrimaryFD` (version=bilinear) and `Petersen2011PrimaryFD` (version=elliptical)
   exercise the "weight-1.0 tree ≡ pre-logic-tree calculator" path.
2. **Aggregation arithmetic.** A 50/50 blend of those two branches produces
   `λ_mean(D₀) = 0.5·λ_A(D₀) + 0.5·λ_B(D₀)` at every D₀, within 1e-12.

## Scenario

- Single strike-slip characteristic fault, `rake = 0°` (see
  `source_model.xml`). The rake override makes the source unambiguously
  `strike-slip` under `openquake.fdha.calc.contexts.classify_style`.
- Single site on the fault trace (`r = 0`, principal case).
- Primary-SR held constant at `Pizza2023PrimarySR` with `style = "all"`,
  weight 1.0.
- Primary-FD is the only branching axis. Both Petersen aliases share the same
  displacement definition (`D_p_L`, per Sarmiento et al. 2025 Table 1 → PEA11)
  so the advisory validator `FDLT-101` does not fire. This keeps the
  comparison physically meaningful.
- Secondary SR/FD: intentionally omitted.

## Layout

```
examples/logic_tree_validation/
├── README.md                    (this file)
├── source_model.xml             (canonical NRML for OpenQuake-style `job.ini`; same
│                                physical model as `source_model_ss.xml`)
├── source_model_ss.xml          (legacy filename kept for tooling that still
│                                references it, e.g. map_mode plots)
├── source_model_logic_tree.xml  (single-branch SMLT, weight=1.0 on `source_model.xml`)
├── verify.py                    (runs the three ``job.ini`` files and checks invariants)
├── single_bilinear/
│   ├── job.ini
│   └── fdha_logic_tree.xml      (1 Primary-FD branch, weight=1.0)
├── single_elliptical/
│   ├── job.ini
│   └── fdha_logic_tree.xml      (1 Primary-FD branch, weight=1.0)
├── blend_50_50/
│   ├── job.ini
│   └── fdha_logic_tree.xml      (2 Primary-FD branches, 0.5 + 0.5)
└── expected/                    (optional; add girs_reference.csv to enable
                                  regression check 5 in verify.py)
```

Each run writes outputs under `<run>/out/`:

```
<run>/out/
├── hazard_curves/branch_XXXX.csv   (D0, annual_rate)
├── aggregate_hazard.csv            (D0, mean, p05, p16, p50, p84, p95)
├── manifest.json
├── validator_report.txt
└── branch_configs/
```

## Runbook

```bash
cd examples/logic_tree_validation/

fdha single_bilinear/job.ini
fdha single_elliptical/job.ini
fdha blend_50_50/job.ini

python verify.py
```

`verify.py` also auto-runs any INI whose `out/aggregate_hazard.csv` is missing,
so `python verify.py` alone is enough for a fresh checkout.

## Checks

| # | What | Criterion |
|---|---|---|
| 1 | D₀ grid consistency | identical displacement levels (from `[calculation].displacement_measure_levels`) across runs |
| 2 | Aggregation arithmetic | `max|X.mean − 0.5·A.mean − 0.5·B.mean| < 1e-12` |
| 3 | Fractile reduction (2-branch 50/50, linear-interp convention) | `p05 = p16 = p50 = min`; `p84 = min + 0.68·(max − min)`; `p95 = min + 0.90·(max − min)` |
| 4 | `blend_50_50` manifest | exactly two branches, weights `[0.5, 0.5]`, both Petersen aliases |
| 5 | GIRS regression (optional) | `max|single_bilinear.mean − girs_reference| < GIRS_TOL` (default 1e-6) if `expected/girs_reference.csv` is present; else `[SKIP]` |

### Fractile convention (design memo §9)

`openquake/fdha/logic_tree/aggregation.py::weighted_fractiles` uses linear
interpolation of the weighted empirical CDF, not a step function. For two
equal-weight branches the CDF is `[0.5, 1.0]`, so:

- `q ≤ 0.5` returns the leftmost value → `min(A, B)`.
- `q > 0.5` interpolates between `(0.5, min)` and `(1.0, max)`:
  `min + (q − 0.5)/0.5 · (max − min)`.

The design memo (§9) explicitly permits updating Check 3 to match the
implementation when linear interpolation is in use; that is what `verify.py`
does.

## Notes

- `verify.py` does not invent a GIRS reference. If `expected/girs_reference.csv`
  is absent it prints `[SKIP]` and proceeds; do not fabricate one.
- If any of the three Petersen aliases is ever de-registered, substitute only
  with another Petersen-family variant that shares the same displacement
  definition (`D_p_L`). Do **not** swap in models with a different definition
  (e.g. `Chiou2025PrimaryFD` which uses `D_sp_Nstar`).
