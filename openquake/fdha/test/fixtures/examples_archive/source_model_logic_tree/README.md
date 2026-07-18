# Source-Model Logic Tree (SMLT) example

This example demonstrates the OpenQuake-aligned **source-model logic
tree** support in PFDHA.  It exercises the new
`[calculation].source_model_logic_tree_file` field in **both** hazard-curve
and hazard-map mode using a single shared SMLT XML.

## Files

| File | Role |
|---|---|
| `source_model_a.xml` | Base source model A (single SS fault, `aValue=4.2`) |
| `source_model_b.xml` | Base source model B (same fault, `aValue=4.5`) |
| `source_model_logic_tree.xml` | Two `sourceModel` branches **plus** a `bGRRelative` branch set on top |
| `fdha_logic_tree.xml` | Minimal FDHA model logic tree (Pizza2023 PSR + Petersen2011 PFD) |
| `fdha_curve.ini` | Hazard-curve INI (single site) |
| `fdha_map.ini` | Hazard-map INI (4×4 km grid) |
| `run.py` | Runs both modes to `out_curve/` & `out_map/` and prints a manifest summary |
| `verify_sm_lt_curve.py` | Checks `Σ w_i λ_i` ≡ LT `mean`; writes `smlt_validation_hazard_curves.png` |

## What the SMLT does

The SMLT enumerates `2 sourceModel × 2 bGRRelative = 4` source-model
realisations, with weights:

```
sm_a × bg0 (Δb=0.0) → 0.40 × 0.5 = 0.20
sm_a × bg1 (Δb=+0.2) → 0.40 × 0.5 = 0.20
sm_b × bg0 (Δb=0.0) → 0.60 × 0.5 = 0.30
sm_b × bg1 (Δb=+0.2) → 0.60 × 0.5 = 0.30
                                   ----
                                   1.00
```

Each is multiplied by every FDHA-LT branch weight to form the final
PFDHA realisation set: `combined_weight = sm_w × fdha_w`.

The `bGRRelative` modification is applied via OpenQuake's own
`openquake.hazardlib.lt.apply_uncertainty`, so every NRML uncertainty
type the engine recognises (`maxMagGRRelative`, `simpleFaultDipRelative`,
`abGRAbsolute`, geometry overrides, …) works identically in PFDHA.

## Run it

Programmatic runner (recommended: **separate** output directories):

```bash
python examples/source_model_logic_tree/run.py
```

From this directory with the **`fdha` CLI** (or `python -m openquake.fdha.main`).
The logic-tree driver's default output folder is `./out/` next to the INI file.

```bash
cd examples/source_model_logic_tree
fdha fdha_curve.ini --plot lt_curve_builtin.png
fdha fdha_map.ini --plot lt_map_builtin.png
```

If you run **curve and map in a row without changing `out/`**, `manifest.json` will
reflect the **last** calculation while older `aggregate_hazard.csv` /
`hazard_curves/` files may be **left over** from the previous mode.  For a
clean layout, pass explicit `outdir`s (as `run.py` does) or keep one run per
directory, e.g. `out_curve/` and `out_map/`.

### Validation (analytic mean vs framework)

After curve outputs exist under `out_curve/`:

```bash
cd examples/source_model_logic_tree
python verify_sm_lt_curve.py
# or regenerate curve outputs inside out_curve/, then validate:
python verify_sm_lt_curve.py --run
```

## Outputs

### Hazard curves (`out_curve/`)

```
out_curve/
├── manifest.json                                  # full audit trail
├── aggregate_hazard.csv                           # SMLT × FDHA mean & fractiles
├── hazard_curves/branch_NNNN.csv                  # per-realisation curves
└── source_model_branches/<idx>_<branch_id>/...    # per-SMLT subdirs
```

### Hazard maps (`out_map/`)

```
out_map/
├── manifest.json
├── aggregate/
│   ├── rates_mean.h5                       # SMLT-weighted mean rates
│   ├── rates_fractiles.h5                  # 5 fractile cubes
│   ├── displacement_map_mean.csv           # mean displacement at return period
│   └── displacement_map_p{05,16,50,84,95}.csv
└── source_model_branches/<idx>_<branch_id>/
    ├── branches/branch_NNNN.h5             # one HDF5 per FDHA branch
    └── aggregate/rates_mean.h5             # per-SMLT-branch mean
```

## How aggregation works

Hierarchical, mirroring the OpenQuake logic-tree composition:

* **Inside** each source-model realisation, the FDHA logic tree is
  collapsed into a per-realisation mean (and fractile) rate cube using
  the FDHA branch weights - the same per-source-id grouped sum used
  before this feature.
* **Across** source-model realisations, the per-realisation cubes are
  combined with the SMLT weights:

  ```
  total_mean(site, D0) = Σ over smlt branches  sm_w × per_smlt_mean(site, D0)
  ```

  The SMLT-weighted fractile cubes are computed analogously
  (recorded in the manifest as an approximation to the full multi-source
  convolution; for a single dominant SMLT branch per site it is exact).

## Backward compatibility

Drop `source_model_logic_tree_file` from `[calculation]` and the run
falls back to the legacy `source_model_file` path with a single implicit
realisation (`branch_id = "sm0"`, `weight = 1.0`).  All existing INIs
keep producing identical output.
