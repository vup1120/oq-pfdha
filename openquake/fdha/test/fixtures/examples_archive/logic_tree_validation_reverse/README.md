# Logic Tree Validation Demo — Reverse faulting

Style-parallel of `examples/logic_tree_validation/` for a **reverse** source.
Same workflow (three INIs, same-structure LT XML, `verify.py`, per-branch
CSVs, plot script); only the faulting style, Primary-FD model choices, and
per-model parameter blocks change.

What it proves, by construction:

1. **Regression.** Two single-branch trees (weight = 1.0) on
   `Youngs2003PrimaryFD` and `Moss2024PrimaryFD` exercise the
   "weight-1.0 tree ≡ pre-logic-tree calculator" path.
2. **Aggregation arithmetic.** A 50/50 blend of those two branches produces
   `λ_mean(D₀) = 0.5·λ_A(D₀) + 0.5·λ_B(D₀)` at every D₀, within 1e-12.

## Scenario

- Single reverse characteristic fault (`rake = 90°`, `dip = 45°`) — see
  `source_model_rv.xml`. Classified as `reverse` under
  `openquake.fdha.calc.contexts.classify_style` (rake ∈ [30°, 150°]).
- Single site on the fault trace (`16.16878333, 39.66247618`), which is the
  midpoint of the first trace segment. Same site as the strike-slip sibling
  so the two demos can be compared side-by-side.
- Primary-SR held constant at `Pizza2023PrimarySR` with `style = "all"`,
  weight 1.0 (style-agnostic, reverse-applicable).
- Primary-FD is the only branching axis: `Youngs2003PrimaryFD` vs
  `Moss2024PrimaryFD`.
- Secondary SR/FD: intentionally omitted.

## Scientific caveats (important!)

1. **`Youngs2003PrimaryFD` parametrised with `style = "all"`, not
   `style = "reverse"`.** The Youngs 2003 formulation was derived for
   normal faulting and the implementation only exposes
   `style ∈ {"all", "normal"}`; the calculator raises `ValueError` on
   `style = "reverse"`. The aggregate fit (`"all"`) is the conventional
   engineering workaround when applying Youngs 2003 outside the
   normal-faulting regime. Use with judgement. The XML comment in
   `single_youngs2003/fdha_logic_tree.xml` carries the same warning.
2. **`Moss2024PrimaryFD` displacement-definition metadata is absent.**
   `Youngs2003PrimaryFD` is annotated in
   `openquake/fdha/logic_tree/model_metadata.yaml` with
   `displacement_definition = D_p_V` (cited from Sarmiento et al. 2025,
   Table 1 YEA03). `Moss2024PrimaryFD` is not in that file, so the
   validator FDLT-101 treats it as UNKNOWN and skips the mixed-definition
   advisory (per v3 §F). If a primary-source citation is later added and
   its definition differs from `D_p_V`, FDLT-101 will begin warning on the
   `blend_50_50` branch set. Do NOT invent a citation to silence it.
3. Both FD models are queried with an **average-displacement**
   conditioning (`norm_disp_type = "AD"` for Youngs,
   `version = "AD"` for Moss) and `source = "EQS"`, `completeness = "all"`
   for Moss 2024 — these keep the two models on a broadly comparable
   conditioning basis, but they are not the same physical quantity.

## Layout

```
examples/logic_tree_validation_reverse/
├── README.md                    (this file)
├── source_model_rv.xml          (shared, rake=90, dip=45; _rv suffix
│                                 avoids collision with examples/source_model.xml)
├── verify.py                    (runs the three INIs and checks invariants)
├── plot_hazard_curves.py        (produces hazard_curves_comparison.png)
├── single_youngs2003/
│   ├── fdha.ini
│   └── fdha_logic_tree.xml      (1 Primary-FD branch, weight=1.0)
├── single_moss2024/
│   ├── fdha.ini
│   └── fdha_logic_tree.xml      (1 Primary-FD branch, weight=1.0)
└── blend_50_50/
    ├── fdha.ini
    └── fdha_logic_tree.xml      (2 Primary-FD branches, 0.5 + 0.5)
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
cd examples/logic_tree_validation_reverse/

fdha single_youngs2003/fdha.ini
fdha single_moss2024/fdha.ini
fdha blend_50_50/fdha.ini

python verify.py
python plot_hazard_curves.py   # -> hazard_curves_comparison.png
```

`verify.py` also auto-runs any INI whose `out/aggregate_hazard.csv` is
missing, so `python verify.py` alone is enough for a fresh checkout.

## Checks

| # | What | Criterion |
|---|---|---|
| 1 | D₀ grid consistency | identical `target_displacement` across runs |
| 2 | Aggregation arithmetic | `max|X.mean − 0.5·A.mean − 0.5·B.mean| < 1e-12` |
| 3 | Fractile reduction (2-branch 50/50, linear-interp convention) | `p05 = p16 = p50 = min`; `p84 = min + 0.68·(max − min)`; `p95 = min + 0.90·(max − min)` |
| 4 | `blend_50_50` manifest | exactly two branches, weights `[0.5, 0.5]`, both expected Primary-FD classes, `style = reverse` |

The fractile convention explanation is identical to the strike-slip sibling
README (§"Fractile convention (design memo §9)"); linear interpolation of
the weighted empirical CDF is used, so for two equal-weight branches the
5–50% fractiles collapse to `min` and the 84/95% fractiles sit at the
interpolated values above.

## Observed numbers (illustrative)

| D₀ (m) | single_youngs2003 | single_moss2024 | blend (LT) | 0.5·A + 0.5·B |
|---:|---:|---:|---:|---:|
| 0.01 | 1.065e-02 | 1.070e-02 | 1.068e-02 | 1.068e-02 |
| 0.10 | 9.843e-03 | 1.058e-02 | 1.021e-02 | 1.021e-02 |
| 1.00 | 4.618e-03 | 4.593e-03 | 4.606e-03 | 4.606e-03 |
| 2.00 | 2.556e-03 | 1.689e-03 | 2.122e-03 | 2.122e-03 |
| 5.00 | 7.809e-04 | 1.738e-04 | 4.773e-04 | 4.773e-04 |

`max|LT mean − 0.5·A − 0.5·B| = 0.000e+00` on this grid.

Notice the crossover around D₀ ≈ 1 m: Moss 2024 predicts a higher rate at
small displacements but a lower rate at large displacements than Youngs
2003, so the tails tell different stories. The logic-tree mean rides
exactly on the arithmetic midpoint by construction, which is the whole
point of the check.
