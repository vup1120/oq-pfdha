# PFDHA examples

This directory intentionally contains only the small, public-facing examples
needed to get started.

## Hazard curve

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

