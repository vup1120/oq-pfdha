# Demo: r_threshold_km as epistemic logic-tree branches (fdhaCalcRThreshold)

Builds on the canonical **`examples/hazard_curve_minimal`** job (Youngs2003
primary/secondary models, same site, lon 16.16573727 ~0.39 km from the fault),
but swaps the source for a demo-local **single M7.0 characteristic event** at a
low annual rate (`source_model_M7.xml`, 1e-4/yr). That keeps the hazard curves
at a realistic level (principal plateau ~1e-4/yr) — comparable to a published
single-scenario figure — instead of reflecting the very active
Gutenberg-Richter fault (~0.02/yr) in the stock example.

`r_threshold_km` is the distance cutoff that decides, **per rupture**, whether
a site is treated as **principal** (essentially on the fault → primary
on-fault displacement, metres of slip) or **distributed** (off-fault → only
small secondary displacement). It is a hard-step simplification of the
continuous rupture-location term *f(r)* of Petersen et al. (2011). The single
M7.0 rupture sits ~0.39 km from the site, so the two demo thresholds produce a
**clean binary flip**:

- `0.1 km` → 0.39 km > 0.1 km → the rupture is **distributed** (plateau ~4e-6/yr).
- `0.5 km` → 0.39 km < 0.5 km → the rupture is **principal** (plateau ~9e-5/yr).

Both values are physically plausible; they, the magnitude, and the rate are
demo inputs, **not** recommendations. The ~20× gap between the two curves is
intrinsic to the principal/distributed dichotomy (on-fault slip is metres;
off-fault displacement is small and rarer), which is exactly why the threshold
choice carries epistemic weight.

Four runs:

| Run      | Mechanism | r_threshold_km                    |
|----------|-----------|-----------------------------------|
| baseline | MODE A    | INI scalar `0.1`                  |
| case 1   | MODE B    | one branch `0.1`, weight 1.0      |
| case 2   | MODE B    | one branch `0.5`, weight 1.0      |
| case 3   | MODE B    | `0.1` (w=0.5) + `0.5` (w=0.5)     |

MODE B variants remove the INI scalar (the scalar and a branch set are
mutually exclusive) and append a `fdhaCalcRThreshold` branch set to the FDHA
logic tree.

Checks asserted by the script:

1. case 1 == baseline **bit-for-bit** (branch value 0.1 == the scalar; MODE A
   and MODE B share one consumption point).
2. case 3 mean == 0.5·case1 + 0.5·case2 at every displacement level, to
   machine precision.
3. case 1 and case 2 must genuinely differ (the threshold changes the routing).

Run:

```
python openquake/fdha/demo/r_threshold_epistemic/run_demo.py
```

Outputs (job variants, per-case results, and the three-panel PNG
`out/r_threshold_epistemic_demo.png`: fault-and-site map | hazard curves |
consistency checks) are written to `out/` next to this file.
