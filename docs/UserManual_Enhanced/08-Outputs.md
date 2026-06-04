# Understanding Outputs

The current public CLI path runs through the logic-tree driver. The optional `--output` file is a small JSON summary; detailed results are written under the logic-tree output directory.

By default, the output directory is `out/` next to the INI job file unless the logic-tree driver is called directly with a different output directory.

---

## CLI Summary JSON

When you run:

```bash
fdha job.ini --output results.json
```

the JSON file contains paths and shape metadata, not the full hazard arrays.

### Hazard Curve Summary

```json
{
  "logic_tree": true,
  "outdir": "/path/to/out",
  "mode": "hazard_curve",
  "n_sites": 1,
  "n_displ": 18,
  "manifest_json": "/path/to/out/manifest.json",
  "validator_report": "/path/to/out/validator_report.txt",
  "aggregate_hazard_csv": "/path/to/out/aggregate_hazard.csv"
}
```

### Hazard Map Summary

```json
{
  "logic_tree": true,
  "outdir": "/path/to/out",
  "mode": "hazard_map",
  "n_sites": 2500,
  "n_displ": 18,
  "manifest_json": "/path/to/out/manifest.json",
  "validator_report": "/path/to/out/validator_report.txt",
  "rates_mean_h5": "/path/to/out/aggregate/rates_mean.h5",
  "rates_fractiles_h5": "/path/to/out/aggregate/rates_fractiles.h5",
  "displacement_map_mean_csv": "/path/to/out/aggregate/displacement_map_mean.csv",
  "target_return_period": 100000.0
}
```

---

## Shared Logic-Tree Files

Every logic-tree run writes:

- `manifest.json`: branch records, weights, selected model classes, source-model branch metadata, validator file paths, and output layout metadata.
- `validator_report.txt`: FDHA logic-tree validation results. A clean file contains `OK`.
- `branch_configs/branch_XXXX.ini`: internal materialized branch INI files. These may contain `[models.*]` sections because they are generated for one end branch.

---

## Hazard Curve Outputs

Hazard-curve runs write:

- `hazard_curves/branch_XXXX.csv`: one CSV per combined branch.
- `aggregate_hazard.csv`: weighted mean and fractile annual exceedance rates.

For one site, branch CSV columns are:

```text
D0,annual_rate
```

For multiple sites, branch CSV columns are:

```text
site_id,lon,lat,D0,annual_rate
```

The aggregate CSV columns are (with the default quantiles):

```text
D0,mean,quantile-0.05,quantile-0.16,quantile-0.5,quantile-0.84,quantile-0.95
```

The `mean` column is present when `[output].mean` is true (the default), and one
`quantile-<q>` column is written per configured `[output].quantiles` value. For
multiple sites, the aggregate CSV includes `site_id`, `lon`, and `lat` before `D0`.

All rates are annual exceedance rates. `D0` values are displacement thresholds in meters.

---

## Hazard Map Outputs

Hazard-map runs write one HDF5 file per end branch plus aggregate HDF5 and CSV products:

- `branches/branch_XXXX.h5`: branch annual-rate grid.
- `aggregate/rates_mean.h5`: weighted mean annual-rate grid.
- `aggregate/rates_fractiles.h5`: fractile annual-rate grids for the configured quantiles (default 0.05, 0.16, 0.50, 0.84, 0.95; see `[output].quantiles`).
- `aggregate/displacement_map_mean.csv`: mean displacement map at the configured return period.
- `aggregate/displacement_map_quantile-0.05.csv`, `quantile-0.16`, `quantile-0.5`, `quantile-0.84`, `quantile-0.95`: fractile displacement maps at the configured return period (one file per configured quantile, named OpenQuake-style `quantile-<q>`).

Displacement-map CSV columns are:

```text
site_id,lon,lat,is_trace,displ_mean
```

Fractile map files use the matching displacement column label, such as `displ_quantile-0.84`.

The first line records the return period as a comment:

```text
# return_period = 100000.0
```

`is_trace` is `1` for trace sites appended by the map site builder and `0` for active grid sites.

---

## Plot Output

Use `--plot output.png` to save a plot.

- Hazard curves plot the logic-tree mean annual exceedance rate versus displacement for site 0, with a 16-84% band when available.
- Hazard maps plot mean displacement by longitude/latitude.
