# Benchmark Tests

This directory contains benchmark tests that compare our implementation against:
1. Published results from research papers
2. Reference implementations (e.g., fdhpy)

## Structure

- **`valentini_et_al_2025/`** - Tests against Valentini et al. (2025) paper results
- **`visini_et_al_2025/`** - Tests against Visini et al. (2025) paper results
- **`mammarella_et_al_2024/`** - Tests against Mammarella et al. (2024) paper results
- **`fdhpy_comparison/`** - Comparison tests against fdhpy reference implementation

## Running Benchmark Tests

```bash
# Run all benchmark tests
pytest openquake/fdha/test/benchmark/ -m benchmark

# Run specific benchmark suite
pytest openquake/fdha/test/benchmark/valentini_et_al_2025/
pytest openquake/fdha/test/benchmark/visini_et_al_2025/
pytest openquake/fdha/test/benchmark/mammarella_et_al_2024/
pytest openquake/fdha/test/benchmark/fdhpy_comparison/
```

## Test Data

Reference data and expected results are stored in each benchmark directory's `data/` or `expected/` subdirectories.

## Notes

- Benchmark tests may be slower than unit/integration tests
- Some tests require reference data files to be present
- Results are compared with published tolerances (see individual test files)


