# r_threshold_km MODE A regression baseline

Golden fixtures freezing the current (pre-epistemic-r_threshold) outputs of
the two public example jobs. They pin MODE A behaviour: a job using the
scalar `[calculation].r_threshold_km` (or relying on the implementation
default when the key is absent) must keep producing these outputs
bit-identically after the `fdhaCalcRThreshold` logic-tree feature lands.

Consumed by
`openquake/fdha/test/integration/logic_tree/test_r_threshold_baseline.py`.

## Provenance

Generated on 2026-07-05 from the code state at commit `f7070be3` by running
`FdhaLogicTree.from_ini(<job>).run(outdir=...)` on:

**`map_default/` regenerated on 2026-07-08** with the same procedure after
the intentional MODE A output change from the original-trace fix: FDHA site
metrics (r, rx sign, x/L) for characteristicFaultSource are now computed
from the exact NRML trace retained at parse time instead of the
mesh-resampled top edge, making them independent of
``rupture_mesh_spacing`` (the map example's characteristic sources have
sinuous traces, so their distances shifted). ``curve_explicit/`` is
unchanged: the curve example uses a simpleFaultSource, whose path is
untouched.

**`curve_explicit/aggregate_hazard.csv`, `curve_explicit/branch_0000.csv`
and `map_default/displacement_map_mean.csv` regenerated on 2026-07-11**
after the principal/distributed component columns were added to the output
files (`mean_principal`/`mean_distributed`,
`annual_rate_principal`/`annual_rate_distributed`,
`displ_mean_principal`/`displ_mean_distributed`). Before re-freezing, the
new files were verified byte-identical to the previous fixtures once the
added columns were stripped, i.e. all previously frozen numbers are
unchanged. The fractile map CSVs and `rates_baseline.npz` did not change.

**`map_default/` regenerated in full on 2026-07-12** after the principal-zone
trace-sampling fix: on-trace (principal) sites are now placed from the exact
``original_trace`` densified to the map grid resolution
(``region_grid_spacing``), instead of the ERF mesh top edge whose node count
collapses when ``rupture_mesh_spacing`` is large. For the minimal map example
this raised the trace-site count 29 -> 98 (denser principal band), so all six
displacement-map CSVs and ``rates_baseline.npz`` (site axis 1594 -> 1663) had
to be re-frozen. Verified surgical before re-freezing: the 1565 distributed
*grid* sites are byte-identical (max |Δ| = 0 on displ_mean/principal/
distributed); only appended trace-site rows changed. ``curve_explicit/`` is a
single-site curve and is unaffected.

**`map_default/` regenerated in full on 2026-07-13** after the principal-zone
trace sampling changed from *densify* to *resample*: the trace is now resampled
at a uniform grid-resolution step instead of retaining every native NRML vertex
(which are often digitised at sub-kilometre spacing, placing far more principal
sites than the grid can resolve — the 48 onshore Taiwan faults produced 3 246
trace sites against a 1 040-site 0.1° grid). For the minimal map example this
lowered the trace-site count 98 -> 68, so all six displacement-map CSVs and
``rates_baseline.npz`` (site axis 1663 -> 1633) were re-frozen. Verified
surgical before re-freezing: the 1565 distributed *grid* sites are
byte-identical (max |Δ| = 0 on displ_mean/principal/distributed); only the
on-trace rows changed. ``curve_explicit/`` is unaffected.

- `curve_explicit/` — `examples/hazard_curve_minimal.ini`
  (sets `r_threshold_km = 0.1` explicitly in `[calculation]`):
  - `aggregate_hazard.csv` — weighted mean + quantile curves (top-level).
  - `branch_0000.csv` — per-branch annual rates
    (`hazard_curves/branch_0000.csv`).
  - `manifest.json` — realization list with weights.
- `map_default/` — `examples/hazard_map_minimal.ini`
  (omits `r_threshold_km`, exercising the implementation default):
  - `displacement_map_mean.csv`, `displacement_map_quantile-*.csv` —
    aggregate displacement maps (top-level `aggregate/`).
  - `manifest.json` — realization list with weights.
  - `rates_baseline.npz` — numeric arrays extracted from the run's HDF5
    outputs: `branch_rates` (from `branches/branch_0000.h5`), `rates_mean`,
    `rates_fractiles` + axes (`d0`, `site_lons`, `site_lats`, `quantiles`)
    and the branch `weight`/`fingerprint`.

Determinism was verified by running each job twice and byte-comparing all
text outputs before freezing.

## Comparison contract

- CSV outputs are deterministic text -> compared **byte-exact**.
- `manifest.json` is compared structurally with absolute filesystem paths
  normalised to basenames (the file embeds the resolved source-model path).
- HDF5-derived rates are compared numerically at **<= 1e-12 relative**
  tolerance.

Do NOT regenerate these fixtures to make a failing test pass unless the
change to MODE A outputs is intentional and reviewed.
