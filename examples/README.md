# PFDHA examples

This directory intentionally contains only the small, public-facing examples
needed to get started.

## Hazard curve

**What it does:** computes the annual frequency of exceeding a range of
fault-displacement levels at a single site, using a one-branch logic tree
built on the Youngs et al. (2003) primary/secondary rupture and displacement
models. This is the simplest end-to-end PFDHA calculation and a good first run
to verify your installation.

Run a single-site hazard curve:

```bash
mkdir -p examples/outputs
fdha examples/hazard_curve_minimal.ini \
  --output examples/outputs/hazard_curve_minimal_results.json \
  --plot examples/outputs/hazard_curve_minimal.png
```

Input files:

- `hazard_curve_minimal.ini`
- `hazard_curve_minimal_source_model_logic_tree.xml`
- `hazard_curve_minimal_fdha_logic_tree.xml`
- `source_model.xml`

## Hazard map

**What it does:** computes fault-displacement hazard over a small region for a
fixed return period, producing a spatial map of displacement at the target
hazard level. It uses the same Youngs et al. (2003) model chain as the hazard
curve and introduces the `region` and `return_period` configuration options.

Run a small hazard map:

```bash
mkdir -p examples/outputs
fdha examples/hazard_map_minimal.ini \
  --output examples/outputs/hazard_map_minimal_results.json \
  --plot examples/outputs/hazard_map_minimal.png
```

Input files:

- `hazard_map_minimal.ini`
- `hazard_map_minimal_source_model_logic_tree.xml`
- `hazard_map_minimal_fdha_logic_tree.xml`
- `source_model_char.xml`

The `examples/outputs/` directory is for generated local results and should not
be committed.

