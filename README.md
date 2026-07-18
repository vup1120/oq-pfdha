# Probabilistic Fault Displacement Hazard Analysis (PFDHA)


This repository contains tools for performing Probabilistic Fault Displacement Hazard Analysis (PFDHA), a methodology to quantify the hazard of earthquake-induced surface fault rupture and fault displacement. It provides a command-line interface (CLI) to run hazard curves and maps.

## Features

-   **Hazard Curve & Hazard Map Calculation**: Compute hazard curves for specific sites or hazard maps over a region.
-   **Extensible Model Library**: Includes a variety of published models for primary and secondary surface rupture and displacement.
-   **Flexible Configuration**: Use INI files (OpenQuake-style) to define calculation parameters, site locations, and model choices.
-   **Standard-based Inputs**: Leverages NRML (XML) format for seismic source models, compatible with tools like the OpenQuake Engine.

## Quickstart

### 1. Installation

This project requires **Python 3.11+**.

```bash

# Install the package in editable mode
pip install -e .
```

Installing in editable mode (`-e`) makes the `fdha` command-line tool available in your environment.

### 2. Run a Demo

A minimal example is provided to quickly verify your installation. This command computes a hazard curve for a single site using the Youngs et al. (2003) model.

```bash
# Run the minimal hazard curve example
mkdir -p examples/outputs
fdha examples/hazard_curve_minimal.ini \
     --plot examples/outputs/hazard_curve_minimal.png
```

The INI configuration file (`examples/hazard_curve_minimal.ini`) includes:
- `[geometry]`: site coordinates (for hazard curves) or a region definition (for hazard maps)
- `[erf]`: rupture mesh spacing and MFD bin width
- `[calculation]`: source-model logic-tree file, FDHA logic-tree file, and displacement levels
- FDHA model choices in `examples/hazard_curve_minimal_fdha_logic_tree.xml`

The calculation type is detected automatically: `sites` in `[geometry]` runs a
hazard curve, `region` runs a hazard map.

Like the OpenQuake engine, results are written to an output directory rather
than to a single file. The directory defaults to `out/` next to the INI file
(here `examples/out/`). For this hazard-curve example the key files are:

-   `examples/out/aggregate_hazard.csv`: the aggregate (mean and quantile) hazard
    curve - the main result.
-   `examples/out/hazard_curves/branch_0000.csv`: the per-logic-tree-branch curve.
-   `examples/out/manifest.json`: an index of everything the run produced.

The `--plot` argument additionally saves a PNG of the hazard curve to the path
you give (here `examples/outputs/hazard_curve_minimal.png`).

See the [Outputs](./docs/UserManual_Enhanced/08-Outputs.md) chapter for the full
file layout and column definitions.

## Documentation

For the complete User Manual (installation, configuration, models, CLI, workflows), see:

-   **[User Manual Index](./docs/UserManual_Enhanced/index.md)**
-   **[Quick Start](./docs/UserManual_Enhanced/02-QuickStart.md)**
-   **[Command-Line Interface (CLI)](./docs/UserManual_Enhanced/04-CLI.md)**
-   **[Configuration File Guide](./docs/UserManual_Enhanced/05-Configuration.md)**
-   **[Scientific Models](./docs/UserManual_Enhanced/06-Models.md)**

## Citing

If you use this software in your research, please cite it using the metadata in
[`CITATION.cff`](./CITATION.cff).

## License

This project is licensed under the GNU Affero General Public License v3.0.
See [`LICENSE`](./LICENSE) for the full text.
