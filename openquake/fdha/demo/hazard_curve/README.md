# Demo: single-site hazard curves with recent FD models

Two ready-to-run jobs on the same fault geometry as
`examples/hazard_curve_minimal`, but with **reverse** rake (+90° instead of
−90°), showcasing the newer primary displacement models beyond the stock
Youngs et al. (2003) chain. The site sits **on the fault trace** (r = 0),
so the principal FD model drives the hazard curve — unlike the minimal
example's off-trace site, where only the distributed models contribute:

- **`job_kuehn2024.ini`** — `Kuehn2024PrimaryFD` with
  `epistemic_uncertainty = true`, so the FD term carries the model's
  posterior-sample ensemble (the job pins `style = normal` via the model
  parameter, which takes precedence over the source's rake-derived style).
  The ensemble is averaged into the hazard curve, so the result is the
  mean hazard with the model's internal epistemic spread fully included —
  no extra configuration is needed.
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
