# Demo: scoping the r_sigma_km epistemic branches to specific sources (applyToSources)

Companion to **`demo/r_sigma_epistemic`** (read that one first for what
`r_sigma_km` is and how a `fdhaCalcRSigma` branch set propagates it). This
demo answers the follow-up question: *can the mapping-accuracy uncertainty
be applied to one or two specific faults only, while every other fault
keeps the default σ = 0 boxcar (`r_threshold_km`)?*

**Yes** — a `fdhaCalcRSigma` branch set accepts the standard
`applyToSources` filter:

```xml
<logicTreeBranchSet branchSetID="bs_5_r_sigma"
                    uncertaintyType="fdhaCalcRSigma"
                    applyToSources="src_A">
  <logicTreeBranch branchID="RS_ACCURATE">
    <uncertaintyModel>0.02689</uncertaintyModel>
    <uncertaintyWeight>0.5</uncertaintyWeight>
  </logicTreeBranch>
  <logicTreeBranch branchID="RS_CONCEALED">
    <uncertaintyModel>0.06552</uncertaintyModel>
    <uncertaintyWeight>0.5</uncertaintyWeight>
  </logicTreeBranch>
</logicTreeBranchSet>
```

This is the physically natural use: mapping accuracy is a property of
**each fault's trace** (one fault well surveyed, another concealed under
sediment), not of the whole source model.

## Setup

Two parallel M7.0 characteristic faults ~0.69 km apart
(`source_model_2src.xml`; fault A is the trace of the sibling demo, fault B
the same trace shifted 0.008° east), and two sites, each **40 m from "its"
fault**. The separation exceeds both near-field scales (boxcar h = 50 m,
widest Gaussian support 2σ = 131 m), so each fault's σ treatment can only
matter at its own site.

Three runs of the *same* two-source job, differing only in the FDHA tree:

| Run      | Sigma branch set                                   | Fault A  | Fault B  |
|----------|----------------------------------------------------|----------|----------|
| baseline | none                                               | boxcar   | boxcar   |
| scoped   | Accurate/Concealed, `applyToSources="src_A"`       | Gaussian | boxcar   |
| unscoped | the same set, **no** filter                        | Gaussian | Gaussian |

## Checks asserted by the script

1. **scoped == baseline at site 2, bit-for-bit** (max rel diff `0.0`) —
   fault B genuinely stayed on the boxcar inside the scoped job; fault A's
   far-field contribution there is σ-independent by construction (beyond
   both h and 2σ the two W_p paths give identical contributions).
2. **scoped == unscoped at site 1 to the last ulp** (~2e-16) — fault A gets
   the same branches whether or not the set is scoped. (Not bit-for-bit:
   the two jobs enumerate different realization sets, so the weighted
   aggregation accumulates in a different order.)
3. scoped ≠ baseline at site 1 (max rel diff ~0.42) — the scoping acted.
4. unscoped ≠ baseline at site 2 (~0.42) — without the filter, fault B's
   site would have changed too.
5. per-source manifest weights each sum to 1 (fault A: 0.5 + 0.5;
   fault B: 1.0). Realizations are **per source**; independent sources are
   *summed*, so weights normalise within each source group, not globally.

## Rules to know (validator FDLT-012)

Each source may be covered by **at most one** `fdhaCalcRSigma` branch set:

- several sets with **disjoint** `applyToSources` scopes may coexist
  (e.g. a wide σ for a concealed fault, a narrow σ for a well-mapped one);
- a set **without** `applyToSources` covers every source and therefore
  cannot coexist with any other sigma set;
- overlapping scopes are rejected at validation, loudly.

A source outside every scope simply never gets a sigma slot: it runs with
the default `r_sigma_km = 0`, i.e. the boxcar with the INI's
`r_threshold_km`.

## Run

```
python openquake/fdha/demo/r_sigma_scoped/run_demo.py
```

Outputs go to `out/` next to this file: the three job variants, per-case
results, and `out/r_sigma_scoped_demo.png` — hazard curves at both sites
for the three runs. At site 1 (near fault A) the scoped curve departs from
the boxcar baseline (the unscoped curve coincides with it, hidden
underneath); at site 2 (near fault B) the scoped curve sits exactly on the
baseline and only the unscoped run departs.

All σ values are Petersen et al. (2011, Tables 2–3) two-sided
mapping-accuracy classes, used here as demo inputs, **not**
recommendations.
