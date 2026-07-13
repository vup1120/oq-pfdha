# Analysis Workflows

This chapter describes the high-level workflows for performing common PFDHA tasks, explaining how the inputs, configuration, and CLI commands work together.

## The PFDHA Calculation Flow

All analyses follow a general sequence, from parsing inputs to generating final results. This process is visualized below.

```mermaid
graph TD
    A[Input: Source-Model Logic Tree XML] --> C{fdha CLI};
    B[Input: FDHA Logic Tree XML] --> C;
    J[Input: Configuration INI] --> C;
    C --> D[1. Parse Inputs, Validate Logic Trees & Setup Sites];
    D --> E[2. Enumerate Source and FDHA End Branches];
    E --> F[3. Unified Hazard Calculation per Branch];
    F --> G[4. Aggregate Weighted Hazard Rates and Fractiles];
    G --> H[5. Generate Final Product];
    H --> I[Output directory: manifest.json, aggregate + per-branch rates<br/>optional PNG plot];

    subgraph "Unified Hazard Calculation (calculate_fdha_hazard)"
        direction LR
        F1[P(Primary Rupture)] --> F2[P(Disp > d | Primary)];
        F2 --> F3[Principal Hazard];
        F4[P(Secondary Rupture)] --> F5[P(Disp > d | Secondary)];
        F5 --> F6[Distributed Hazard];
        F3 --> F7[Total Hazard Rate];
        F6 --> F7;
    end

    style A fill:#f9f,stroke:#333,stroke-width:2px
    style B fill:#f9f,stroke:#333,stroke-width:2px
    style I fill:#ccf,stroke:#333,stroke-width:2px
```

**Key Steps Explained:**

1.  **Parse Inputs**: The tool reads the INI job, the source-model logic tree, and the FDHA model logic tree. It sets up the site location for a hazard curve or a grid of sites for a hazard map. The calculation type is automatically detected from `[geometry]`.
2.  **Enumerate Branches**: The driver expands source-model realizations and FDHA model end branches, validates branch weights and model class names, and materializes internal branch configs.
3.  **Unified Hazard Calculation**: For every site and every branch, the tool uses `calculate_fdha_hazard()` to compute the combined hazard from both principal (on-fault) and distributed (off-fault) displacement using the four models selected in the FDHA logic tree.
4.  **Aggregate Hazard Rates**: Hazard contributions are combined with branch weights to produce mean and fractile hazard curves.
5.  **Generate Final Product**:
    -   For a **hazard curve**, this aggregated curve is the final result.
    -   For a **hazard map**, the tool uses the hazard curve at each site to find the displacement that corresponds to the target `return_period` through interpolation.

---

## Workflow 1: Calculating a Hazard Curve

This is the most fundamental workflow, used to assess the hazard at a specific point or small area.

**Goal**: To calculate the annual frequency of exceedance for a range of displacement levels at **one or more discrete sites**.

!!! note "Multiple Sites"
    A hazard-curve job may evaluate several sites in a single run. List them in `[geometry].sites` as comma-separated `lon lat [depth]` entries, or load them from a CSV with `[geometry].sites_csv`. Per-site curves are written with `site_id`, `lon`, and `lat` columns; see [Outputs](08-Outputs.md).

**Steps**:

1.  **Prepare Inputs**: Create a source model XML file and a configuration file (INI format recommended).
2.  **Configure INI**:
    -   Define a `[geometry]` section with `sites` for one or more point sites (format: `lon1 lat1, lon2 lat2, ...`), or `sites_csv` for sites from a CSV file.
    -   Define the `displacement_measure_levels` in `[calculation]` section (JSON format: `{"FD": [0.001, 0.01, ...]}`).
    -   Set `[calculation].source_model_logic_tree_file` and `[calculation].fdha_logic_tree_file`.
    -   Select the four required model categories in the FDHA logic-tree XML, not in `[models]` INI sections; to model only one side, neutralize the other with the constant `Fixed*SR` models (see [Models — Modeling only principal or only distributed displacement](06-Models.md#modeling-only-principal-or-only-distributed-displacement)).
3.  **Run Calculation**:
    Execute the `fdha` command with your configuration file:
    ```bash
    fdha job_curve.ini --plot curve.png
    ```
    The calculation type is automatically detected from the configuration (presence of `region` indicates hazard map, otherwise hazard curve).
4.  **Interpret Results**: Results are written to the output directory (`out/`
    next to the INI by default). The hazard curve data is written as CSV —
    `aggregate_hazard.csv` (aggregate mean/quantile curve) and
    `hazard_curves/branch_XXXX.csv` (per branch). The `--plot` PNG visualizes the
    mean curve. See [Outputs](08-Outputs.md).

---

## Workflow 2: Generating a Hazard Map

This workflow visualizes the spatial distribution of hazard across a region for a specific probability level.

**Goal**: To calculate the displacement level that is expected to be exceeded at a given return period (e.g., 100,000 years) for every point on a grid.

**Steps**:

1.  **Prepare Inputs**: Prepare your source model(s) and a configuration file (INI format recommended).
2.  **Configure INI**:
    -   Define a `[geometry]` section with `region` to define the rectangular area for the grid (format: `lon1 lat1, lon2 lat2, lon3 lat3, lon4 lat4`).
    -   **Note**: For hazard maps, you must use `region`.
    -   Set `region_grid_spacing` (in degrees) to control the grid resolution.
    -   Set the desired `return_period` in the `[calculation]` section (e.g., `100000` years).
    -   Define `displacement_measure_levels` in `[calculation]` section.
    -   Set `[calculation].source_model_logic_tree_file` and `[calculation].fdha_logic_tree_file`.
    -   Select the four required model categories in the FDHA logic-tree XML, not in `[models]` INI sections; to model only one side, neutralize the other with the constant `Fixed*SR` models (see [Models — Modeling only principal or only distributed displacement](06-Models.md#modeling-only-principal-or-only-distributed-displacement)).
3.  **Run Calculation**:
    ```bash
    fdha job_map.ini --plot map.png
    ```
    The calculation type is automatically detected from the presence of `region` in the `[geometry]` section.
4.  **Interpret Results**: Results are written to the output directory (`out/`
    next to the INI by default). The gridded results are the displacement maps as
    CSV (`aggregate/displacement_map_mean.csv` and one
    `displacement_map_quantile-<q>.csv` per configured quantile) and the
    underlying rate grids as HDF5 (`aggregate/rates_mean.h5`,
    `aggregate/rates_fractiles.h5`). The `--plot` PNG shows the mean map colored by
    displacement. See [Outputs](08-Outputs.md).

---

## References (OQ)

-   OpenQuake Engine User Guide — User Guide Index: https://docs.openquake.org/oq-engine/manual/latest/user-guide/index.html
