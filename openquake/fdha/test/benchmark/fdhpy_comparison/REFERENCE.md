# Reference for fdhpy comparison

## Cited reference

Sarmiento, A., Lavrentiadis, G., Bozorgnia, Y., Chen, R., Chiou, B., Dawson, T., Kottke, A., Kuehn, N., Madugo, C., Moss, R., Thompson, S., and Zandieh, A. (2025). Comparisons of FDHI fault displacement models for principal and aggregate displacement. *Earthquake Spectra*, 41(4), 2691-2720. DOI: 10.1177/87552930251327894.

## Source of the numerical values

Live outputs from the installed `fdhpy` package during pytest comparisons, pinned to `fdhpy==1.0.3` in the `test` extras of `pyproject.toml` (`pip install -e .[test]`). Earlier runs used the editable checkout at commit `0ca38e269f99937477d9a43f331e4a432cf41b7c` (1.0.2) with a local NumPy-2 shim; 1.0.3 ships that fix. The Lavrentiadis and Abrahamson (2023) outputs are also frozen from 1.0.3 as golden fixtures in [FDHI_Tests/fixtures/la23](FDHI_Tests/fixtures/la23/README.md), so those comparisons do not need `fdhpy` at test time.

## What this case validates

This benchmark compares pfdha primary fault-displacement probability and aleatory-parameter implementations against the external `fdhpy` package for Youngs et al. (2003), Petersen et al. (2011), Moss et al. (2024), Kuehn et al. (2024), Lavrentiadis and Abrahamson (2023), and Chiou et al. (2025) model interfaces.

## Tolerance and its justification

`tests/test_simple_reference.py` uses `rtol=1e-4`, `atol=1e-8` for Youngs2003 exceedance probabilities and `rtol=1e-6`, `atol=1e-10` for most other exceedance comparisons. The aleatory tests inspected use `rtol=1e-6`, `atol=1e-10`. The README describes Youngs2003 as relaxed for integration differences and all others as stricter equivalence checks; the installed `fdhpy` reference is pinned to version 1.0.3.

## How to reproduce

```bash
pytest openquake/fdha/test/benchmark/fdhpy_comparison/FDHI_Tests/tests -q
```

## Status on HEAD

198/198 PASS on 2026-07-19 against the pinned reference checkout `0ca38e2` (with a local NumPy-2 compatibility shim on the reference side; see [RESULTS.md](RESULTS.md)).
