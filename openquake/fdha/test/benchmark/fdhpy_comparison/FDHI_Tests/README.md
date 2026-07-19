# PFDHA Reference Test Suite

Comparison tests between `pfdha` (this implementation) and `fdhpy` (reference
implementation). See [../RESULTS.md](../RESULTS.md) for the latest results and
[../REFERENCE.md](../REFERENCE.md) for the reference citation and pinned commit.

## Structure

```
FDHI_Tests/
├── pytest.ini                     # pytest configuration
├── README.md                      # This file
└── tests/
    ├── test_simple_reference.py   # Exceedance-probability comparisons (111 cases)
    ├── test_aleatory_petersen2011.py
    ├── test_aleatory_kuehn2024.py
    ├── test_aleatory_lavrentiadis2023.py
    ├── conftest.py                # pytest fixtures
    ├── comparison_helpers.py      # Helper utilities
    └── __init__.py
```

## Run

```bash
pip install -e /path/to/fdhpy   # reference, commit pinned in ../REFERENCE.md
python -m pytest tests -q       # from this directory
```

Tests skip cleanly when `fdhpy` is not installed.

## Models tested

| Model | Notes |
|-------|-------|
| Youngs et al. (2003) | D/AD version, style="all" |
| Petersen et al. (2011) | elliptical, quadratic versions |
| Moss et al. (2024) | D/AD, D/MD with GIRS/EQS |
| Kuehn et al. (2024) | normal, reverse, strike-slip |
| Lavrentiadis & Abrahamson (2023) | aggregate metric |
| Chiou et al. (2025) | model7, model8.2 |
