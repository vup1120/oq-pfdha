# Configuration guide for Visini et al. (2025)

> NOTICE — current scope: In this build, the Visini (2025) secondary models are wired into hazard-curve calculations only. Hazard maps or other calculators are not yet supported for Visini (2025) SR/FD.

The Visini secondary models are activated through the FDHA logic-tree XML referenced by the public INI job. This section lists only the public `uncertaintyModel` parameters accepted by `Visini2025SecondarySR` and `Visini2025SecondaryFD`; rupture context values such as magnitude, distance, and combination are supplied by the calculator.

## Model selection

Select `Visini2025SecondarySR` in an FDHA logic-tree branch with `uncertaintyType="fdhaSecondarySRModel"` and `Visini2025SecondaryFD` in a branch with `uncertaintyType="fdhaSecondaryFDModel"`.

## Rank 1.5 trace definitions
Define every nearby Rank 1.5 splay under `[rank1p5_ruptures]`. Each `trace` entry supplies a name and geometry.

## Global parameters

These parameters are specified in the `[calculation]` section (INI format):

| Name | Type | Units | Default | Allowed | Required? | Description | Dependencies |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `displacement_measure_levels` | JSON string | meters | – | positive numbers | Yes | Displacement thresholds for hazard curves. Used to size probability arrays. Use JSON format: `{"FD": [0.001, 0.01, ...]}`. | Needed for both primary and secondary displacement calls. |
| `case` | string | – | `case1` | `case1`, `case2`, `case3` | Yes | User-mandated combination set: Case 1 → A/B/C, Case 2 → A/B, Case 3 → A only. | Drives which combinations feed the Visini probabilities. |
| `near_far_threshold_km` | float | kilometers | `0.2` | >0 | No | Distance cutoff between "near" and "far" regimes for along‑strike Monte Carlo (SR Rank 2). | Affects `calculate_rank2_total_probability` inputs. |
| `r_threshold_km` | float | kilometers | `0.1` | >0 | No | Distance split between principal (≤ threshold) and distributed (> threshold) branches in the hazard curve. | Only used in hazard‑curve integration. |
| `primary_sr_reduction` / `secondary_sr_reduction` | JSON string | – | median, 50th percentile | method∈{`median`,`mean`,`percentile`} | No | Reduction statistic for Monte Carlo arrays. Use JSON format: `{"method": "median", "q": 50}`. | `secondary_sr_reduction` defaults to `primary_sr_reduction`.

## `Visini2025SecondarySR` parameters
| Name | Type | Units | Default | Allowed | Required? | Description | Dependencies |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `style` | string | – | inferred from source rake when omitted | `normal`, `reverse` | No | Dip‑slip style family. Dip‑slip only. | Must match the earthquake mechanism. |
| `pixel_size` | integer | meters | – | {10, 20, 50, 100, 200, 500} | Yes | **Across-strike width**: Site width perpendicular to the fault. Used for P(across) coefficient lookup (Table 2) and F-ratio lookup (Table 3). | Choose to reflect site/slice footprint. |

## `Visini2025SecondaryFD` parameters
| Name | Type | Units | Default | Allowed | Required? | Description | Dependencies |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `style` | string | – | inferred from source rake when omitted | `normal`, `reverse` | No | Mirrors the rupture mechanism indicator (normal=1, reverse=0). | Keep consistent with rupture style. |
| `scaling_model` | string | – | `WC1994` | `WC1994`, `THINGBAIJAM2017`, `LEONARD2010` | No | Selects the scaling relation used to compute TPFm when not provided explicitly. | Affects throw smoothing window and magnitude–slip link. |
| `tpfm` | float or list[float] | meters | computed | >0 | No | Directly sets mean throw on the principal fault; when set, it overrides the scaler‑derived value. | If set, `scaling_model` is ignored for those sites. |
| `truncation_eps` | float | standard deviations | `3.0` | >0 | No | Half-width of the truncated lognormal residual distribution. | Advanced use only. |

### Units and sign conventions
- Distances `s` (FD), `r` (SR slice minimum distance), and `rx` are in meters in the Visini regressions.
- Internally, the hazard‑curve integrator computes site–fault distances in kilometers, but converts to meters before calling the Visini SR/FD models. The principal/distributed mask uses `r_threshold_km`.
- In both SR and FD, `rx < 0` denotes footwall (FW); `rx ≥ 0` denotes hanging wall (HW). This toggles the FW indicator in the regressions.
- Displacement metric Y is the vertical throw (net vertical component) of Rank 2 distributed ruptures.

### Supported rupture mechanisms
- Dip‑slip only: `style ∈ {normal, reverse}`.
- Not supported: strike‑slip and oblique‑slip styles.

## Rank 1.5 traces & Combination B (≤ 1 km)
Combination selection follows the paper’s intent. In this implementation, users explicitly set `case` in the public INI; do not set `combination` inside the FDHA logic-tree `uncertaintyModel`. When considering Combination B (near‑site Rank 1.5 control):
- You should declare nearby Rank 1.5 splays under `[rank1p5_ruptures]` with `trace.name` and `trace.geometry` (type `Line`, lon/lat pairs).
- Apply Combination B as a candidate only when the site lies within about 1 km of the declared Rank 1.5 trace (per the paper’s recommended applicability ranges, cf. Fig. 7). Beyond ~1 km, prefer Combination A (Rank 1 association).
- Screen faults near the site as potential hosts for Rank 1.5 ruptures that can produce Rank 2 distributed rupture at the site.

### Rank 1.5 definitions (INI format)

For INI format, Rank 1.5 traces can be defined in two ways:

1. **Using XML file** (recommended for INI):
```ini
[calculation]
rank1p5_traces_file = rank1p5_traces.xml
```

2. **Using INI configuration** (for simple cases):
```ini
[rank1p5_ruptures]
trace = [{"name": "Splay_A", "geometry": {"type": "Line", "coords": [[16.338416, 39.640532], [16.351064, 39.653699]]}}]
```

## Case selection (user‑specified)

Set exactly one of:
- `case = case1` → enable combinations A, B, C
- `case = case2` → enable combinations A, B
- `case = case3` → enable combination A only

### Decision aid
- Use `case1` when Rank 1.5 splays at/near the site are plausible and should be considered (A+B+C).
- Use `case2` when splays are near but not mapped directly at the site (A+B only).
- Use `case3` when no Rank 1.5 influence is expected within ~1 km (A only).

## TPFm options (direct vs. scaling) in `Visini2025SecondaryFD.get_prob`

1) Direct TPFm value (overrides scaling):

**FDHA logic-tree uncertainty model block**:
```ini
[Visini2025SecondaryFD]
style = reverse
tpfm = 1.2
```

2) Scaling model estimation (default):

```ini
[Visini2025SecondaryFD]
style = reverse
scaling_model = WC1994
```
Internally, the scaling path computes an along‑strike mean throw profile and smooths it using a half‑window equal to ½·s (with s in km), capped by 0.5 of fault length in normalized coordinates.

## Complete Configuration Examples

### Public INI job plus FDHA logic-tree branches

```ini
[general]
description = Visini_secondary_example
calculation_mode = fdha_classical

[geometry]
sites = 16.16573727 39.64704451

[site_params]
reference_vs30_value = 600.0

[erf]
rupture_mesh_spacing = 1.0
width_of_mfd_bin = 0.1

[calculation]
source_model_logic_tree_file = source_model_logic_tree.xml
fdha_logic_tree_file = fdha_logic_tree.xml
displacement_measure_levels = {"FD": [0.0001, 0.001, 0.005, 0.01, 0.015, 0.03, 0.05, 0.075, 0.1, 0.15, 0.3, 0.5, 0.75, 1.0, 3.0, 5.0, 7.5, 10.0]}
case = case3
near_far_threshold_km = 0.2
r_threshold_km = 0.1
rank1p5_traces_file = rank1p5_traces.xml
```

```xml
<logicTreeBranchSet branchSetID="bs_secondary_sr"
                    uncertaintyType="fdhaSecondarySRModel">
  <logicTreeBranch branchID="B_VISINI2025_SECONDARY_SR">
    <uncertaintyModel><![CDATA[[Visini2025SecondarySR]
style = reverse
pixel_size = 50]]></uncertaintyModel>
    <uncertaintyWeight>1.0</uncertaintyWeight>
  </logicTreeBranch>
</logicTreeBranchSet>

<logicTreeBranchSet branchSetID="bs_secondary_fd"
                    uncertaintyType="fdhaSecondaryFDModel">
  <logicTreeBranch branchID="B_VISINI2025_SECONDARY_FD">
    <uncertaintyModel><![CDATA[[Visini2025SecondaryFD]
style = reverse
scaling_model = LEONARD2010
# tpfm = 1.2]]></uncertaintyModel>
    <uncertaintyWeight>1.0</uncertaintyWeight>
  </logicTreeBranch>
</logicTreeBranchSet>
```

This consolidates the configuration pathways for the Visini (2025) models and clarifies units, case selection, Rank 1.5 usage, pixel‑size binning, site geometry, and TPFm options.

## Site Footprint Configuration

The public hazard-curve configuration supports a point site. For the Visini secondary rupture model, `pixel_size` defines the square site footprint used by the model.

```ini
[Visini2025SecondarySR]
style = normal
pixel_size = 100  # 100m × 100m square site
```
