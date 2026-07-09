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
