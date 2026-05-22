# Reference for Norcia sensitivity Youngs2003

## Cited reference

International Atomic Energy Agency. (2025). *Benchmarking Current Practices in Probabilistic Fault Displacement Hazard Analysis for Nuclear Installations*. IAEA-TECDOC-2092, Vienna. DOI: 10.61092/iaea.74us-dn4n.

## Source of the numerical values

The repository hard-codes the Youngs2003 AD 85 reference rate vector in the Norcia comparison scripts, and Ye confirmed this embedded vector as the accepted baseline for this repository. The scientific context is the IAEA TECDOC-2092 Norcia principal fault displacement exercise, with benchmark inputs summarized in Tables 13-14 and the hazard-curve comparison shown in Figure 21.

## What this case validates

This benchmark compares a Norcia source-model sensitivity case using the Youngs2003 AD 85 model configuration against the principal-fault-displacement hazard-curve context documented for the Norcia exercise in IAEA TECDOC-2092.

## Tolerance and its justification

The comparison workflow on the `tecdoc2092` branch records `rtol=1e-2`, `atol=1e-6` for the Norcia vector in `compare_results.py`; the current HEAD plotting scripts compute differences against the same embedded `ref_Youngs_N` vector but do not expose that comparison as a pytest assertion.

## How to reproduce

```bash
python openquake/fdha/test/benchmark/norcia_sensitivity_youngs2003/run_test.py
```

## Status on HEAD

BLOCKED on 2026-05-12, commit `fd3cb88`: the embedded baseline is confirmed by Ye, but the current HEAD workflow remains orchestration-script based rather than a pytest assertion.
