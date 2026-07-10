# Command-Line Interface (CLI)

The PFDHA toolkit is operated through the `fdha` command-line interface. The current public interface runs INI jobs and dispatches through the logic-tree driver.

## Synopsis

```bash
fdha <job_ini> [--output <output_file>] [--plot [plot_file]]
```

## Arguments

| Argument | Type | Default | Required? | Description |
| :--- | :--- | :--- | :--- | :--- |
| `job_ini` | string | none | **Yes** | Path to the v5 canonical INI job file. |
| `--output`, `-o` | string | none | No | Save a JSON or CSV summary. In logic-tree mode this is a convenience summary that points to the output directory and aggregate files. |
| `--plot`, `-p` | optional string | none | No | Generate a plot. Use `--plot` with no value to display it, or `--plot curve.png` to save it. |
| `--verbose`, `-v` | flag | `False` | No | Enable verbose logging. |
| `-h`, `--help` | flag | `False` | No | Show help and exit. |

!!! info "Required public INI keys"
    The job file must include `[calculation].fdha_logic_tree_file` and `[calculation].source_model_logic_tree_file`. Legacy public `source_model_file`, `fdha_logic_tree_files`, and `[models.*]` sections are rejected by the loader.

## Automatic Calculation Type Detection

The calculation type is detected from `[geometry]`:

- If `[geometry]` contains `region`, the job is a hazard map.
- Otherwise the job is a hazard curve.

Hazard-curve jobs take one or more point sites from `sites` (or from a
`sites_csv` file); hazard maps must use `region`.

## Examples

### 1. Calculate a Hazard Curve

Configuration file (`job_curve.ini`):

```ini
[general]
description = hazard_curve_example
calculation_mode = fdha_classical

[geometry]
sites = 16.16573727 39.64704451

[erf]
rupture_mesh_spacing = 2.0
width_of_mfd_bin = 0.1

[calculation]
investigation_time = 1.0
displacement_measure_levels = {"FD": [0.001, 0.01, 0.1, 0.5, 1.0]}
fdha_logic_tree_file = hazard_curve_minimal_fdha_logic_tree.xml
source_model_logic_tree_file = hazard_curve_minimal_source_model_logic_tree.xml
```

Model choices are stored in the FDHA logic-tree XML, for example:

```xml
<logicTreeBranchSet branchSetID="bs_1_primary_surf_rup"
                    uncertaintyType="fdhaPrimarySRModel">
  <logicTreeBranch branchID="B1_PRIMARY_SURF_RUP">
    <uncertaintyModel>
      [Youngs2003PrimarySR]
      style = all
    </uncertaintyModel>
    <uncertaintyWeight>1.0</uncertaintyWeight>
  </logicTreeBranch>
</logicTreeBranchSet>
```

Command:

```bash
fdha job_curve.ini --output results.json --plot curve.png
```

### 2. Calculate a Hazard Map

Configuration file (`job_map.ini`):

```ini
[general]
name = hazard_map_example

[geometry]
region = 15.9 39.2, 16.3 39.2, 16.3 39.8, 15.9 39.8
region_grid_spacing = 0.01

[erf]
rupture_mesh_spacing = 2.0
width_of_mfd_bin = 0.1

[calculation]
investigation_time = 1.0
displacement_measure_levels = {"FD": [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]}
return_period = 100000
max_distance_km = 10.0
fdha_logic_tree_file = hazard_map_minimal_fdha_logic_tree.xml
source_model_logic_tree_file = hazard_map_minimal_source_model_logic_tree.xml
```

Command:

```bash
fdha job_map.ini --output map_results.json --plot map.png
```

### 3. Run Without Extra Output

```bash
fdha job.ini
```

### 4. Display a Plot

```bash
fdha job.ini --plot
```
