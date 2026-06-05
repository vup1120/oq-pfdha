# Reference for Valentini et al. (2025) / Kumamoto Case 2

## Cited reference

International Atomic Energy Agency. (2025). *Benchmarking Current Practices in Probabilistic Fault Displacement Hazard Analysis for Nuclear Installations*. IAEA-TECDOC-2092, Vienna. DOI: 10.61092/iaea.74us-dn4n. The repository directory is retained under `valentini_et_al_2025/` because Ye plans to expand it beyond the current Kumamoto Case 2 contents.

## Source of the numerical values

Published benchmark reference vector `REF_EXCEED` embedded in `test_chiou2025_case2.py` and reused by `plot_kumamoto_logic_tree_vs_paper.py`, corresponding to the Kumamoto principal-fault-displacement Case 2 exercise implemented for the Chiou (2025) primary displacement branch in this repository.

## What this case validates

This benchmark runs the Kumamoto Case 2 configuration for the Chiou (2025) primary fault-displacement model and compares the computed exceedance curve against the bundled paper reference vector for the same displacement grid.

## Tolerance and its justification

The coded tolerance in `test_chiou2025_case2.py` is `rtol=1.2e-1`, `atol=0.0`. The test docstring says the tolerance allows minor geometric and numerical differences across environments and notes an observed maximum relative difference of about 11%; the precision of the published source vector is not otherwise documented.

## How to reproduce

```bash
pytest openquake/fdha/test/benchmark/valentini_et_al_2025/test_chiou2025_case2.py -q
```

## Status on HEAD

BLOCKED on 2026-05-12, commit `fd3cb88`: the test invokes the `fdha` CLI, which failed in this environment with `cannot cache function 'idx_start_stop': no locator available for file '<oq-engine>/openquake/baselib/performance.py'`.
