# Demo: single-site hazard curves with recent FD models

Two ready-to-run jobs on the same site and fault geometry as
`examples/hazard_curve_minimal`, but with **reverse** rake (+90° instead of
−90°), showcasing the newer primary displacement models beyond the stock
Youngs et al. (2003) chain:

- **`job_kuehn2024.ini`** — `Kuehn2024PrimaryFD` with
  `epistemic_uncertainty = true`, so the FD term carries the model's
  posterior-sample ensemble (the job pins `style = normal` via the model
  parameter, which takes precedence over the source's rake-derived style).
  The `[parameters]` section reduces those Monte Carlo samples at the
  **84th percentile**
  (`primary_sr_reduction = {"method": "percentile", "q": 84}`) instead of
  the default median — a conservative single-curve summary of the model's
  internal epistemic spread.
- **`job_moss2024.ini`** — `Moss2024PrimaryFD` (journal `source = EQS`
  coefficients, `version = MD`, all-data completeness) paired with
  `Moss2013PrimarySR` (reverse, Vs30 600 m/s).

Run either from the repository root:

```bash
fdha openquake/fdha/demo/hazard_curve/job_kuehn2024.ini --plot kuehn2024.png
fdha openquake/fdha/demo/hazard_curve/job_moss2024.ini --plot moss2024.png
```

Logic-tree outputs (manifest, per-branch and aggregate curves) are written to
`out/` next to the INI files; `out/` is git-ignored.
