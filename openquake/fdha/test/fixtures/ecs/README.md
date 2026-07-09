# ECS validation fixtures (R `mgcv` oracle → pure-Python port)

These fixtures are the **ground truth** for validating the pure-Python ECS
reference-line module against the R reference implementation
(`ecs_functions.R` / `CalcEventCoordinateSystem.R`, Lavrentiadis et al. 2024).

**Runtime never invokes R.** R + `mgcv` is run *once* on a developer machine as a
throwaway oracle to emit the CSVs below; thereafter the committed CSVs are the
ground truth and the test suite diffs the Python implementation against them. R
is re-run only if the spline code itself changes.

## Directory layout

```
ecs/
├── README.md                      ← this file (the fixture contract)
├── <event>/                       ← one dir per validation event
│   ├── flatfile_measurements.csv  ← INPUT: displacement points (R FLATFILE_MEASUREMENTS)
│   ├── flatfile_ruptures.csv      ← INPUT: rupture vertices    (R FLATFILE_RUPTURES)
│   ├── ecs_trace.csv              ← ORACLE OUT: final ECS trace (R ecs_trace)
│   ├── lpmatrix.csv               ← ORACLE OUT: mgcv model matrix at final iteration
│   └── penalty_S.csv              ← ORACLE OUT: mgcv penalty matrix S (per smooth)
```

## Validation gates (both must pass before any forward x/L use)

- **Gate 1 — real event:** a displacement+rupture event with R's `ecs_trace`.
  Exercises the full `ecs_main` pipeline (disp + rupture weighting, MRS/PCA start,
  iteration). Confirms the GC2 reuse + iteration control + spline reconstruction
  agree with R end-to-end.
- **Gate 2 — vertices-only:** the same pipeline fed *only* rupture vertices
  (`abs_wt_rup` weights, `start_sol='MRS'`) — the actual forward source-model
  path used for `multiFaultSource`.

## Fixture contract (exact columns)

### Inputs (provided by you, copied verbatim from the R flatfiles)
- `flatfile_measurements.csv`: must include `EQ_ID, longitude_degrees,
  latitude_degrees, rank, <field_disp_wt>` where `<field_disp_wt>` is the
  displacement-weight field used in `CalcEventCoordinateSystem.R`
  (default `recommended_net_preferred_for_analysis_meters`).
- `flatfile_ruptures.csv`: must include `EQ_ID, longitude_degrees,
  latitude_degrees, RUP_ID, NODE_ID, rank`.

### Oracle outputs (from the R run)
- `ecs_trace.csv`: exactly the data frame `ecs_trace` written by
  `CalcEventCoordinateSystem.R` — columns
  `EQ_ID, eq_name, REF_ID, Longitude, Latitude, ecs_u, ecs_t, ecs_curv`.
- `lpmatrix.csv`: the `mgcv` linear-predictor (model) matrix for the **final**
  iteration's GAM fit, i.e. `predict(fit_gam_xy, type = "lpmatrix")`. This is the
  numeric realization of the `tp` thin-plate basis the Python side must reproduce.
- `penalty_S.csv`: the smoothing penalty matrix `S` for the `s(u)` term, i.e.
  `fit_gam_xy$smooth[[1]]$S[[1]]` (and the scalar `sp = lambda_p / fault_len`
  recorded in a header comment), so the Python penalty matches mgcv's.

## Minimal R extraction snippet (dev-only oracle, not shipped at runtime)

Run inside `update_ecs`/`ecs_main` at the final iteration (or re-fit on the
converged `fault_disp`) with the project's existing R functions sourced:

```r
# fit_gam_xy is the mgcv::gam object from update_ecs() on the converged ECS
write.csv(predict(fit_gam_xy, type = "lpmatrix"),
          file.path(out, "lpmatrix.csv"), row.names = FALSE)
write.csv(as.matrix(fit_gam_xy$smooth[[1]]$S[[1]]),
          file.path(out, "penalty_S.csv"), row.names = FALSE)
# header note: sp used = lambda_p / fault_len   (record the numeric value)
write.csv(ecs_trace, file.path(out, "ecs_trace.csv"), row.names = FALSE)
```

## Tolerances (proposed; revisit at Gate 1)

- `ecs_trace` Longitude/Latitude: agreement to within the R iteration's own
  convergence threshold `flt_max_ds = 50 m` (`ecs_main`), reported in km.
- `lpmatrix`, `penalty_S`: `assert_allclose(rtol=1e-6)` after resolving mgcv's
  basis centering/identifiability constraint (the `-1` no-intercept in the R
  formula) — this is the precise check that the Python `tp` reconstruction is
  faithful rather than merely close.
