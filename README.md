# Probabilistic Fault Displacement Hazard Analysis (PFDHA)


This repository contains tools for performing Probabilistic Fault Displacement Hazard Analysis (PFDHA), a methodology to quantify the hazard of earthquake-induced surface fault rupture and fault displacement. It provides a command-line interface (CLI) to run hazard curves and maps.

## Features

-   **Hazard Curve & Hazard Map Calculation**: Compute hazard curves for specific sites or hazard maps over a region.
-   **Extensible Model Library**: Includes a variety of published models for primary and secondary surface rupture and displacement.
-   **Flexible Configuration**: Use INI files (recommended, OpenQuake-style) to define calculation parameters, site locations, and model choices.
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

A minimal example is provided to quickly verify your installation. This command computes a hazard curve for a single site using the Youngs (2003) model.

**Note**: The configuration file now includes all required parameters following OpenQuake Engine conventions.

```bash
# Run the minimal hazard curve example
mkdir -p examples/outputs
fdha examples/hazard_curve_minimal.ini \
     --output examples/outputs/hazard_curve_minimal_results.json \
     --plot examples/outputs/hazard_curve_minimal.png
```

The INI configuration file (`examples/hazard_curve_minimal.ini`) includes:
- `[erf]` section: Rupture mesh spacing and MFD bin width
- `[calculation]` section: source-model logic-tree file, FDHA logic-tree file, and displacement levels
- `[geometry]`: Site coordinates (for hazard curves) or region definition (for hazard maps)
- FDHA model choices in `examples/hazard_curve_minimal_fdha_logic_tree.xml`

**Note**: The calculation type (hazard curve or hazard map) is automatically detected from the configuration file:
- If `[geometry]` contains `region`, it's a hazard map calculation
- If `[geometry]` contains `sites`, it's a hazard curve calculation

After running, you should find the following files:
-   `examples/outputs/hazard_curve_minimal_results.json`: The calculated hazard curve data.
-   `examples/outputs/hazard_curve_minimal.png`: A plot of the hazard curve.

## Documentation

For the complete User Manual (installation, configuration, models, CLI, workflows), see:

-   **[User Manual Index](./docs/UserManual_Enhanced/index.md)**
-   **[Quick Start](./docs/UserManual_Enhanced/02-QuickStart.md)**
-   **[Command-Line Interface (CLI)](./docs/UserManual_Enhanced/04-CLI.md)**
-   **[Configuration File Guide](./docs/UserManual_Enhanced/05-Configuration.md)**
-   **[Scientific Models](./docs/UserManual_Enhanced/06-Models.md)**
-   **[Validation suite](./VALIDATION.md)**




## Citing

If you use this software in your research, please cite it using the metadata in
[`CITATION.cff`](./CITATION.cff).

## License

This project is licensed under the GNU Affero General Public License v3.0.
See [`LICENSE`](./LICENSE) for the full text.


