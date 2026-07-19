# pfdha vs fdhpy comparison results

**Date:** 2026-07-19
**Reference:** `fdhpy` at commit `0ca38e2` (see [REFERENCE.md](REFERENCE.md)), installed editable
**Result:** 198 / 198 tests PASS

## Summary

All six primary fault-displacement model implementations in pfdha reproduce the
external `fdhpy` reference, both for exceedance probabilities and for the
inspected aleatory parameters (mu, sigma, statistical-distribution internals).

### Exceedance probabilities (`tests/test_simple_reference.py`, 111 cases)

| Model | Cases | Tolerance | Status |
|-------|------:|-----------|--------|
| Youngs et al. (2003) | 12 | rtol 1e-4 | PASS |
| Petersen et al. (2011) | 18 | rtol 1e-6 | PASS |
| Moss et al. (2024) | 24 | rtol 1e-6 | PASS |
| Kuehn et al. (2024) | 27 | rtol 1e-6 | PASS |
| Lavrentiadis & Abrahamson (2023) | 18 | rtol 1e-6 | PASS |
| Chiou et al. (2025) | 12 | rtol 1e-6 | PASS |

Youngs (2003) uses a relaxed tolerance because the two implementations
integrate the normalized-displacement mixture differently; all other models
are strict numerical-equivalence checks.

### Aleatory parameters (87 cases)

| Test file | Cases | Status |
|-----------|------:|--------|
| `test_aleatory_petersen2011.py` | 36 | PASS |
| `test_aleatory_kuehn2024.py` | 30 | PASS |
| `test_aleatory_lavrentiadis2023.py` | 21 | PASS |

Tolerance rtol 1e-6, atol 1e-10.

## Environment note

Under NumPy 2.x, `fdhpy` at the pinned commit raises a `TypeError` in
`chiou_et_al_2025.py` (`float()` on a 1-element record array) that its
`_required` decorator silently converts to `None`, failing the 12 Chiou
cases. A one-line local shim in the fdhpy checkout (`float(...)` ->
`(...).item()`, numerically identical and matching the idiom used elsewhere
in that file) restores correct behaviour. This is a reference-side
compatibility issue, not a pfdha discrepancy.

## Reproduce

```bash
pip install -e /path/to/fdhpy   # commit 0ca38e2
pytest openquake/fdha/test/benchmark/fdhpy_comparison/FDHI_Tests/tests -q
```
