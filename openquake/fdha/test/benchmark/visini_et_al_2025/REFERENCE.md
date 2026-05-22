# Reference for Visini et al. (2025)

## Cited reference

Visini, F., Boncio, P., Valentini, A., Scotti, O., Nurminen, F., Baize, S., and Pace, B. (2025). Empirical regressions for distributed faulting of dip-slip earthquakes. *Earthquake Spectra*, 41(4), 2968-3001. DOI: 10.1177/87552930241308860.

## Source of the numerical values

Digitised CSV reference data under `fig13/reference_data/visini2025_case*.csv` for Visini Figure 13; `test_fig14.py` reproduces Visini Figure 14 in plot form and optionally overlays a CSV if one is present, but no Figure 14 CSV is tracked in the repository.

## What this case validates

This benchmark exercises the Visini et al. (2025) secondary surface-rupture and displacement model components, the decision-tree selection of case combinations, and Figure 13 reproduction workflows comparing logic-tree hazard curves for cases 1-3 against digitised reference curves.

## Tolerance and its justification

`test_visini_integration.py` mostly checks monotonicity, finite positive values, and exact case-set selection, without explicit numeric paper tolerances. `fig13/plot_fig13_cases_vs_reference.py` reports relative error plots but does not enforce a pytest tolerance. `test_fig14.py` is a reproduction script collected by pytest and does not define a test tolerance. Missing enforced tolerances should be decided by Ye.

## How to reproduce

```bash
pytest openquake/fdha/test/benchmark/visini_et_al_2025 -q
```

## Status on HEAD

BLOCKED on 2026-05-12, commit `fd3cb88`: pytest collection/execution is blocked in this environment by the OpenQuake/numba cache error for model imports and by script-style collection in `test_fig14.py` and the standalone Fig. 13 helper script `fig13/run_case1_after_fix.py`.
