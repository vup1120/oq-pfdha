# Quick Start Guide

This guide walks through a minimal example to run your first PFDHA calculation. This will verify that your installation is working correctly and introduce you to the basic workflow.

The goal is to compute a fault displacement hazard curve for a single site using a sample configuration and source model provided in the repository.

## Running the Demo

!!! note "Prerequisites"
    Please ensure you have successfully completed all the steps in the [Installation](01-Installation.md) guide before proceeding.

1.  **Navigate to the Repository Root**:
    Make sure you are in the root directory of the `oq-pfdha` project where you cloned it. If you used a virtual environment, ensure it is activated.

2.  **Create an Output Directory**:
    It's good practice to have a dedicated directory for your results.
    ```bash
    mkdir -p examples/outputs
    ```

3.  **Execute the `fdha` Command**:
    Run the following command in your terminal. It computes a hazard curve using a pre-configured INI file. The source model path is specified in the configuration file.

    ```bash
    fdha examples/hazard_curve_minimal.ini \
         --plot examples/outputs/hazard_curve_minimal.png
    ```

## Understanding the Inputs

The command you just ran used two main input files:

1.  **Configuration File (INI format, recommended)**:
    This INI file specifies all the parameters for the calculation. It defines the site location, the displacement levels of interest, the source-model logic tree, and the FDHA model logic tree. For a full breakdown of the options, see the [Configuration](05-Configuration.md) guide.

    ```ini
    [general]
    description = minimal_hazard_curve
    calculation_mode = fdha_classical

    # --- Site Location ---
    [geometry]
    sites = 16.16573727 39.64704451

    [erf]
    rupture_mesh_spacing = 2.0
    width_of_mfd_bin = 0.1

    # --- Calculation Parameters ---
    [calculation]
    investigation_time = 1.0
    displacement_measure_levels = {"FD": [0.0001, 0.001, 0.005, 0.01, 0.015, 0.03, 0.05, 0.075, 0.1, 0.15, 0.3, 0.5, 0.75, 1.0, 3.0, 5.0, 7.5, 10.0]}
    r_threshold_km = 0.1
    fdha_logic_tree_file = hazard_curve_minimal_fdha_logic_tree.xml
    source_model_logic_tree_file = hazard_curve_minimal_source_model_logic_tree.xml

    # --- Output Statistics ---
    [output]
    mean = true
    quantiles = 0.05 0.16 0.5 0.84 0.95
    ```

    !!! note "Configuration File Format"
        The PFDHA toolkit uses INI format (OpenQuake-style). Source models are referenced through `[calculation].source_model_logic_tree_file`, and FDHA scientific models are selected in `[calculation].fdha_logic_tree_file`.

2.  **Source Model File**:
    This is a standard NRML XML file that describes the earthquake source—in this case, a simple fault. It defines the fault's geometry, magnitude-frequency distribution, and other seismological parameters. The path to this file is wrapped by the source-model logic tree referenced in `[calculation].source_model_logic_tree_file`. See the [Inputs](03-Inputs.md) guide for more information.

## Reviewing the Outputs

Like the OpenQuake engine, `fdha` writes its results to an **output directory**,
not to a single file. The directory defaults to `out/` next to the INI file — so
for this example it is `examples/out/`. (The `--plot` argument additionally saves
a PNG of the hazard curve to the path you gave,
`examples/outputs/hazard_curve_minimal.png`.)

For a hazard **curve** job, the key files in the output directory are:

| File | Contents |
|------|----------|
| `aggregate_hazard.csv` | The aggregate (mean/quantile) hazard curve — **the main result** |
| `hazard_curves/branch_0000.csv` | Per-logic-tree-branch hazard curve |
| `manifest.json` | Index of everything produced by the run |

This particular INI also carries a hazard-**map** source-model branch, so the run
additionally writes map outputs under `aggregate/` and
`source_model_branches/.../aggregate/`:

-   `displacement_map_mean.csv` and `displacement_map_quantile-0.05|0.16|0.5|0.84|0.95.csv`
    — displacement maps as CSV.
-   `rates_mean.h5`, `rates_fractiles.h5` — the underlying rate grids in HDF5.

To learn about the output format in detail, refer to the [Outputs](08-Outputs.md) chapter.

## Next Steps

Congratulations, you have successfully run your first PFDHA calculation! You are now ready to explore more advanced topics:

-   Learn about the **[Command-Line Interface](04-CLI.md)** in detail.
-   Dive deep into the **[Configuration](05-Configuration.md)** options to customize your own analyses.
-   Explore the available **[Scientific Models](06-Models.md)**.
