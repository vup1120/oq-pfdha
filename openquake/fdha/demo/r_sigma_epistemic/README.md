# Demo: r_sigma_km (mapping accuracy) as epistemic logic-tree branches (fdhaCalcRSigma)

Builds on the canonical **`examples/hazard_curve_minimal`** job (Youngs2003
primary/secondary models), computes **two sites bracketing the 50 m boxcar
edge** (site A at 0.04 km, site B at 0.08 km from the trace), and swaps the
source for a demo-local **single M7.0 characteristic event** at a
low annual rate (`source_model_M7.xml`, 1e-4/yr). That keeps the hazard curves
at a realistic level (plateau ~1e-4/yr) - comparable to a published
single-scenario figure - instead of reflecting the very active
Gutenberg-Richter fault (~0.02/yr) in the stock example.

`r_sigma_km` is the **two-sided mapping-accuracy sigma** of the
rupture-location weight W_p (Petersen et al. 2011, Tables 2-3). It selects one
of two separate W_p paths (`docs/design/rupture_location_uncertainty.md`):

- `σ = 0` → **boxcar** `1{|r| ≤ r_threshold_km}` with the historical
  **complementary** split: inside the half-width **only principal**, outside
  **only distributed** (Youngs 2003 / Takao 2013 either/or).
- `σ > 0` → **Petersen's pure Gaussian** `exp(−r²/2σ²)`, pinned to 1 on the
  trace, truncated at ±2σ, **summed** with the full distributed term
  (Petersen eq. 1 + eq. 2; `r_threshold_km` plays no role on this path).

**Both branches carry the same 50 m location knowledge, treated
differently** - hard cliff vs soft tail - and the two sites show BOTH sides
of that trade:

- **site A, r = 0.04 km (inside the boxcar):** σ=0 gives the **full
  principal and nothing else** (either/or, plateau 8.7e-5 /yr); the four
  classes give W_p = 0.33/0.66/0.83/0.86 × principal + distributed
  (plateaus 3.4e-5 → 8.0e-5 /yr) → the hard edge sits **above** every
  Gaussian class inside the band.
- **site B, r = 0.08 km (outside the boxcar):** σ=0 gives **distributed
  only** (5.4e-6 /yr); Accurate (2σ = 54 m < 80 m) coincides with it
  exactly, while Approximate/Concealed/Inferred give W_p =
  0.19/0.48/0.55 × principal + distributed (up to 5.3e-5 /yr) → the hard
  edge **understates** the hazard just outside the cliff.

σ = 0.05 km sits inside Petersen's own two-sided mapping classes
(0.027–0.116 km); treating the class choice as weighted branches follows
Petersen et al. (2011, p. 811). All values are demo inputs, **not**
recommendations.

Six runs:

| Run      | Mechanism | r_sigma_km                              |
|----------|-----------|-----------------------------------------|
| baseline | MODE A    | absent (default 0 = boxcar)             |
| case 1   | MODE B    | one branch `0.0`, weight 1.0            |
| classes  | MODE B    | one branch each: Accurate `0.02689`, Approximate `0.04382`, Concealed `0.06552`, Inferred `0.07269` (Petersen Tables 2–3, Fig.-9c-style greys) |

MODE B variants append a `fdhaCalcRSigma` branch set to the FDHA logic tree;
the INI never carries an `r_sigma_km` scalar (the scalar and a branch set are
mutually exclusive). The boxcar half-width `r_threshold_km = 0.05` stays in
the INI - it is a fixed calculation parameter, not part of the epistemic
tree, and it only matters on the σ = 0 branch.

Checks asserted by the script:

1. case 1 == baseline **bit-for-bit** (branch value 0 == the σ=0 default;
   MODE A and MODE B share one consumption point).
2. the widest class (Inferred) must genuinely differ from σ=0 at both sites.

(The weighted-mean linearity anchor lives in
`test/integration/logic_tree/test_r_sigma_epistemic.py`, V7 of the design
doc.)

Run:

```
python openquake/fdha/demo/r_sigma_epistemic/run_demo.py
```

Outputs (job variants, per-case results, and the four-panel PNG
`out/r_sigma_epistemic_demo.png`: fault-and-sites map | site A curves |
site B curves | consistency checks) are written to `out/` next to this
file.
