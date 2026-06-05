# Logic Tree Validation — map mode

Map-mode analogue of the sibling curve-mode demo (`../`). Goal: **verify
that the logic-tree framework produces the same hazard map as running
each constituent model individually and post-processing the weighted
average ourselves**, on a single strike-slip fault.

The demo mirrors the curve-mode workflow one-for-one:

| Step | Curve mode (`../`) | Map mode (here) |
|------|---|---|
| 1. Individual-model runs | `single_bilinear/`, `single_elliptical/` (weight 1.0) | `single_bilinear_map/`, `single_elliptical_map/` (weight 1.0) |
| 2. Logic-tree blend run | `blend_50_50/` | `blend_50_50_map/` |
| 3. Numerical check | `verify.py` (rate-space 50/50 match) | `verify_map.py` (rate-space AND displacement-space checks) |
| 4. Visual check | `plot_hazard_curves.py` → `hazard_curves_comparison.png` | `plot_hazard_maps.py` → `hazard_maps_comparison.png` |

Every file in this folder exists exclusively to make that equivalence
observable and testable.

## Inputs

| File | Role |
|------|------|
| `../source_model_ss.xml` | Shared strike-slip source (one fault, rake = 0°). The same source the curve-mode demo uses. |
| `geometry_4x4.ini` | Reference `[geometry]` block; duplicated verbatim in the three scenario INIs so the grid, D0 axis, and site ordering match exactly. |
| `single_bilinear_map/job.ini` | Weight-1.0 LT on `Petersen2011PrimaryFD_bilinear` (reuses `../single_bilinear/fdha_logic_tree.xml`). |
| `single_elliptical_map/job.ini` | Weight-1.0 LT on `Petersen2011PrimaryFD_elliptical` (reuses `../single_elliptical/fdha_logic_tree.xml`). |
| `blend_50_50_map/job.ini` | 50/50 LT over both models (reuses `../blend_50_50/fdha_logic_tree.xml`). |

Each LT XML is shared with the curve-mode sibling, so the weights and
model parameters are guaranteed to be byte-identical across modes.

## Outputs (v4 layout)

Each `job.ini` writes under its own `out/`:

```
<scenario>/out/
├── branch_configs/branch_XXXX.ini         # per-branch virtual INI
├── per_source_models/source_<sid>.xml     # per-source subset
├── branches/branch_XXXX.h5                # (n_sites, n_d0) rates + attrs
├── aggregate/
│   ├── rates_mean.h5                      # (n_sites, n_d0)
│   ├── rates_fractiles.h5                 # (5, n_sites, n_d0),
│   │                                      # order [0.05, 0.16, 0.50, 0.84, 0.95]
│   ├── displacement_map_mean.csv          # one row per site, '# return_period =' header
│   ├── displacement_map_p05.csv
│   ├── displacement_map_p16.csv
│   ├── displacement_map_p50.csv
│   ├── displacement_map_p84.csv
│   └── displacement_map_p95.csv
├── manifest.json
└── validator_report.txt
```

On this demo's 4×4-ish grid the driver resolves **28 sites** (16 active
grid sites + 12 trace sampling sites) and **9 D0 levels**.

## Running the validation

```bash
cd examples/logic_tree_validation/map_mode

# Runs the three scenarios if their outputs are missing; otherwise
# just re-runs the numerical checks over the cached outputs.
python verify_map.py

# Reproduces the 2x2 comparison figure.
python plot_hazard_maps.py
```

`verify_map.py` runs five checks. The interesting ones from a
"framework == post-processing" point of view are:

- **Step 3** — rate-space 50/50 identity on every site and every D0:

      lambda_blend_LT(i, j) == 0.5 * lambda_A(i, j) + 0.5 * lambda_B(i, j)

  Tolerance 1e-12. This is the map-mode counterpart of the rate-space
  check done in the curve-mode `verify.py`.

- **Step 5** — end-to-end post-processing check in displacement space.
  Runs A alone and B alone, takes their `rates_mean.h5` grids,
  computes the analytical blend `0.5·lambda_A + 0.5·lambda_B`, inverts
  it at RP = 2475 yr with
  `openquake.fdha.calc.utils.interpolation.get_map_from_curves`, and
  compares the result to the LT framework's own
  `blend_50_50_map/out/aggregate/displacement_map_mean.csv`.
  Tolerance 1e-10 m.

Because aggregation is strictly on **rates**, the post-processing path
is `rates -> blend -> invert`. Averaging displacements directly would
be physically wrong and is deliberately not tested.

## Reference result (current codebase)

```
STEP 2  Grid / D0 / site-ordering consistency
    reference: single_bilinear_map  shape=(28, 9)  n_sites=28  n_d0=9
    OK: all three runs share grid, D0, site ordering.

STEP 3  Per-site arithmetic (test 7): blend = 0.5 A + 0.5 B on every site
    max|blend - 0.5*(A+B)| per-site = 0.000e+00  (tol 1e-12)

STEP 4  Inversion consistency (test 8)
    scenario: single_bilinear_map       max|manual - csv| = 0.000e+00
    scenario: single_elliptical_map     max|manual - csv| = 0.000e+00
    scenario: blend_50_50_map           max|manual - csv| = 0.000e+00

STEP 5  End-to-end post-processing check
    max|displ_LT_blend - displ_analytical_blend| = 0.000e+00 m  (tol 1e-10)
    OK: LT blend hazard map matches per-model runs + analytical weighted average.
```

The exact `0.000e+00` max differences are by design — the LT framework
and the post-processing share the same rate samples and the same
log-log inversion, so any non-zero gap would indicate a real bug.

`hazard_maps_comparison.png` shows the four panels side-by-side:

1. `(a) single_bilinear_map`
2. `(b) single_elliptical_map`
3. `(c) LT framework blend_50_50.mean`
4. `(d) analytical invert[0.5·lambda_A + 0.5·lambda_B]`

Panels (c) and (d) are pixel-for-pixel identical, which is the visual
restatement of the Step-5 numerical result.

## Related Phase G / v4 tests

The pytest versions of the checks above live under
`openquake/fdha/test/integration/logic_tree/`:

- `test_map_per_site_arithmetic.py` — Step 3 (per-site 50/50).
- `test_map_inversion_consistency.py` — Step 4 (single-scenario inversion).
