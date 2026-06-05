# Taiwan FDHA Logic-Tree Hazard Map — unified demo

One INI drives the whole calculation. The logic tree is expressed in three
NRML XML files (one per faulting style) and the source model is a single
file holding one representative fault per style.

## Logic-tree structure (per the user's image)

Each `logicTreeBranchSet` uses equal weights across its branches.

### Reverse (applies to Chelungpu, source id = 17)

| Level | Branches | Equal weight |
|-------|----------|--------------|
| Primary SR  | WC1993, Youngs2003, MossRoss2011, Takao2013, Moss2013, Pizza2023, MammarellaEtAl2024 | 1/7 each |
| Primary FD  | Youngs2003 [AD], Youngs2003 [MD], Takao2013 [AD], Takao2013 [MD], Moss2024 [AD], Moss2024 [MD] | 1/6 each |
| Secondary   | 3 locked (SR, FD) pairs: (Youngs, Youngs), (Takao2014, Youngs), (Visini, Visini) | 1/3 each |

End-branches for this style: **7 × 6 × 3 = 126**.

### Strike-slip (applies to Meishan, source id = 20)

| Level | Branches | Equal weight |
|-------|----------|--------------|
| Primary SR  | WC1993, Youngs2003, Takao2013, Pizza2023, MammarellaEtAl2024 | 1/5 each |
| Primary FD  | Youngs2003 [AD], Youngs2003 [MD], Petersen2011 (bilinear), Takao2013 [AD], Takao2013 [MD] | 1/5 each |
| Secondary   | Full Cartesian 3 SSR × 2 SFD = 6 pairs. SSR = {Youngs2003, Petersen2011, Takao2014}, SFD = {Youngs2003, Petersen2011} | 1/3 × 1/2 |

End-branches for this style: **5 × 5 × 6 = 150**.

### Normal (applies to Shanchiao, source id = 01)

| Level | Branches | Equal weight |
|-------|----------|--------------|
| Primary SR  | WC1993, Youngs2003, Pizza2023, MammarellaEtAl2024 | 1/4 each |
| Primary FD  | Youngs2003 [AD], Youngs2003 [MD] | 1/2 each |
| Secondary   | 2 locked (SR, FD) pairs: (Youngs, Youngs), (Visini, Visini) | 1/2 each |

End-branches for this style: **4 × 2 × 2 = 16**.

**Grand total end-branches: 126 + 150 + 16 = 292.**

Each style's weights sum to 1.0 independently (so weight per end-branch = 1 /
n_selections-within-its-style).

## Files

| File | Purpose |
|------|---------|
| `fdha_map_taiwan.ini` | Single entry point; references the three LT XMLs below. |
| `fdha_logic_tree_reverse.xml` | Reverse tree (126 end-branches). |
| `fdha_logic_tree_strike_slip.xml` | Strike-slip tree (150 end-branches). |
| `fdha_logic_tree_normal.xml` | Normal tree (16 end-branches). |
| `taiwan_representative_3faults.xml` | Shanchiao (normal), Chelungpu (reverse), Meishan (strike-slip) extracted from `F_onshore_1_M_YM_mean_S_mean.xml`. |
| `F_onshore_1_M_YM_mean_S_mean.xml` | Full Taiwan onshore source model (unchanged; kept for reference). |
| `hazard_map_B1-3.ini` | Legacy single-model job (unchanged; reference for parameters). |
| `run_taiwan_hazard_map.py` | Runs the LT and renders the mean displacement map. |

## How the single INI imposes the logic tree

The crux: **each `logicTreeBranchSet` carries `applyToStyle="…"`**, so the
Cartesian enumeration inside the driver only activates a branch set when the
source's rake-derived style matches. A single merged spec therefore handles
all three faulting styles simultaneously.

At run time, the driver:

1. Parses the three LT XMLs, merges their branching levels, and validates
   (FDLT-001 … FDLT-007).
2. Enumerates Cartesian end-branches **per source**:
   16 for Shanchiao (normal), 126 for Chelungpu (reverse), 150 for Meishan
   (strike-slip).
3. For each end-branch, writes a per-source-subset XML
   (`out_taiwan_lt/per_source_models/source_<sid>.xml`) and runs the FDHA
   kernel against that source only. This guarantees that, e.g., reverse
   branch Moss2024+WC1993 is **not** applied to the strike-slip Meishan.
4. Builds a combined SiteCollection (grid inside the region, filtered by
   `max_distance_km`, plus exact fault-trace points).
5. Aggregates rates correctly for independent faults:

       total_mean(site, D0) = Σ_sources  weighted-mean-within-source-LT(site, D0)

   Fractiles are the per-source weighted empirical 5/16/50/84/95 summed
   across sources (exact for single-dominant-source sites, approximate for
   sites where more than one fault contributes).
6. Inverts each site's total-mean rate curve at RP = 2475 yr to get the
   hazard-map displacement.

## Outputs (v4 layout)

```
out_taiwan_lt/
├── branch_configs/branch_XXXX.ini        # per-branch INI (for debugging)
├── per_source_models/source_<sid>.xml    # single-source subsets
├── branches/branch_XXXX.h5               # per end-branch per-site rates
│                                         #   rates (n_sites, n_d0); + weight,
│                                         #   fingerprint, source_id, style
├── aggregate/
│   ├── rates_mean.h5                     # LT-weighted mean rate grid
│   ├── rates_fractiles.h5                # (5, n_sites, n_d0) in order
│   │                                     #   [0.05, 0.16, 0.50, 0.84, 0.95]
│   ├── displacement_map_mean.csv         # mean inverted at RP=2475 yr
│   ├── displacement_map_p05.csv          # one CSV per aggregation level
│   ├── displacement_map_p16.csv
│   ├── displacement_map_p50.csv
│   ├── displacement_map_p84.csv
│   └── displacement_map_p95.csv
├── manifest.json                         # branches + fingerprints + weights
│                                         #   + advisories + layout + quantiles
└── validator_report.txt
```

Per-branch per-site rates are stored as compressed HDF5 (one file per
end-branch); the aggregated mean and fractile grids are HDF5 as well.
Only the derived displacement maps are CSV (one per aggregation level,
with a `# return_period = T` header comment).

## Running

```bash
cd examples/logic_tree_validation_taiwan
python run_taiwan_hazard_map.py
```

With 292 end-branches × 419 sites × 24 D0 levels, this typically takes
~20–40 minutes on a reference box (hardware dependent). Produces
`taiwan_hazard_map_mean.png`.

If you need a quicker turnaround (e.g. for iterating on the tree design),
either drop secondary branches temporarily or coarsen the grid
(`region_grid_spacing = 0.2` cuts site count by 4× without changing the LT).

## Notes / model-specific caveats

- `Youngs2003PrimarySR` only ships with `style = {all, normal}`; for the
  reverse and strike-slip trees we pass `style = "all"` to stay inside its
  documented envelope.
- `Youngs2003PrimaryFD` same convention (`style = "all"` for reverse/SS;
  `style = "normal"` for normal).
- `Visini2025SecondarySR` requires `style ∈ {normal, reverse}` and
  `pixel_size`; we set these explicitly in the reverse tree
  (`style = "reverse", pixel_size = 100`) and the normal tree
  (`style = "normal", pixel_size = 100`).
- `Moss2013PrimarySR` requires `vs30` explicitly (the adapter does not
  forward `[site_params] vs30` to primary-SR models). We pass
  `vs30 = 521.0` to match `[site_params] reference_vs30_value = 521`.
- `Petersen2011PrimaryFD`: the user image shows a single "Petersen2011 AD"
  branch in the strike-slip PFD column. Petersen et al. (2011) is
  parameterised with shape variants (bilinear / elliptical / quadratic)
  rather than AD/MD, so we use the **bilinear** flavour as the
  representative (`Petersen2011PrimaryFD_bilinear`). Switch to
  `_elliptical` or `_quadratic` in `fdha_logic_tree_strike_slip.xml` if
  you want a different interpretation.
- Grid spacing of 0.1° (~10 km) is coarse by design — the demo is about
  *demonstrating that the LT framework works on a real-world three-style
  setup*, not delivering a production map. Drop spacing to 0.05° / 0.02°
  in `fdha_map_taiwan.ini` if you want a finer map (expect ~4× / ~25× run
  time respectively).
