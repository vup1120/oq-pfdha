# Understanding Outputs

Like the OpenQuake engine, `fdha` writes its results to an **output directory**
rather than to a single result file. By default this directory is `out/` next to
the INI job file (unless the logic-tree driver is called directly with a different
output directory). Everything described in this chapter lives inside it.

The sections below cover the shared logic-tree files, the hazard-curve outputs,
the hazard-map outputs, and the optional plot.

---

## Shared Logic-Tree Files

Every logic-tree run writes:

- `manifest.json`: branch records, weights, selected model classes, source-model branch metadata, validator file paths, and output layout metadata.
- `source_model_branches/<NN_branchid>/validator_report.txt`: FDHA logic-tree validation results, one per source-model branch. A clean file contains `OK`.
- `source_model_branches/<NN_branchid>/branch_configs/branch_XXXX.ini`: internal materialized branch INI files. These may contain `[models.*]` sections because they are generated for one end branch.

---

## Hazard Curve Outputs

Hazard-curve runs write:

- `hazard_curves/branch_XXXX.csv`: one CSV per combined branch.
- `aggregate_hazard.csv`: weighted mean and fractile annual exceedance rates.

For one site, branch CSV columns are:

```text
D0,annual_rate,annual_rate_principal,annual_rate_distributed
```

For multiple sites, branch CSV columns are:

```text
site_id,lon,lat,D0,annual_rate,annual_rate_principal,annual_rate_distributed
```

`annual_rate` is the total annual exceedance rate;
`annual_rate_principal` and `annual_rate_distributed` are the contributions
from principal (on-fault) and distributed (off-fault) rupture, so users can
inspect both components at the site. They satisfy
`annual_rate = annual_rate_principal + annual_rate_distributed` exactly.

The aggregate CSV columns are (with the default quantiles):

```text
D0,mean,mean_principal,mean_distributed,quantile-0.05,quantile-0.16,quantile-0.5,quantile-0.84,quantile-0.95
```

The `mean`, `mean_principal` and `mean_distributed` columns are present when
`[output].mean` is true (the default), and one `quantile-<q>` column is
written per configured `[output].quantiles` value. The logic-tree mean is
linear in the branch rates, so `mean_principal + mean_distributed = mean`
exactly; quantiles are reported for the total hazard only. For multiple
sites, the aggregate CSV includes `site_id`, `lon`, and `lat` before `D0`.

All rates are annual exceedance rates. `D0` values are displacement thresholds in meters.

---

## Hazard Map Outputs

Hazard-map runs write one HDF5 file per end branch plus aggregate HDF5 and CSV products:

- `branches/branch_XXXX.h5`: branch annual-rate grid. Datasets: `rates`
  (total), `rates_principal`, `rates_distributed` (per-component grids with
  `rates = rates_principal + rates_distributed`), plus the `d0`,
  `site_lons`, `site_lats` axes.
- `aggregate/rates_mean.h5`: weighted mean annual-rate grid. Datasets:
  `rates_mean`, `rates_mean_principal`, `rates_mean_distributed` (the mean
  is linear, so the components sum exactly to the total), plus axes.
- `aggregate/rates_fractiles.h5`: fractile annual-rate grids for the configured quantiles (default 0.05, 0.16, 0.50, 0.84, 0.95; see `[output].quantiles`). Fractiles are computed for the total hazard only.
- `aggregate/displacement_map_mean.csv`: mean displacement map at the configured return period.
- `aggregate/displacement_map_quantile-0.05.csv`, `quantile-0.16`, `quantile-0.5`, `quantile-0.84`, `quantile-0.95`: fractile displacement maps at the configured return period (one file per configured quantile, named OpenQuake-style `quantile-<q>`).

Mean displacement-map CSV columns are:

```text
site_id,lon,lat,is_trace,displ_mean,displ_mean_principal,displ_mean_distributed
```

`displ_mean` inverts the site's mean total-rate curve at the target return
period. `displ_mean_principal` and `displ_mean_distributed` invert each
component's own mean rate curve, i.e. they answer "what displacement has
this return period considering only principal / only distributed
faulting". The inversion is nonlinear, so the two component columns do
**not** sum to `displ_mean`.

Fractile map files contain only the total column, with the matching label such as `displ_quantile-0.84`.

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

---

## Programmatic access

To locate outputs from a script, read `manifest.json` in the output directory -
it records the run mode, branches, weights, and output layout. The result files
themselves have stable names (`aggregate_hazard.csv`, `aggregate/rates_mean.h5`,
etc.) relative to the output directory documented above.
