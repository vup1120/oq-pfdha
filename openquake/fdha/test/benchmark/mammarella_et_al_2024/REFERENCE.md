# Reference for Mammarella et al. (2024)

## Cited reference

Mammarella, L., Visini, F., Boncio, P., Baize, S., Scotti, O., Beauval, C., Pace, B., and Thompson, S. (2025). Conditional probability of surface rupture: A numerical approach for principal faulting. *Earthquake Spectra*, 41(4). First published online December 20, 2024. DOI: 10.1177/87552930241293570. Developer-provided CPSR reference code repository: `MammarellaLisa/CPSR`.

## Source of the numerical values

CPSR.m developer-code output stored in `expected/mammarella2024_cpsr.csv`; CPSR public repository `main` resolved on 2026-05-12 to commit `d106b1e1e5f124bff9b784369ad1ca284b733aeb`.

## What this case validates

This benchmark checks `Mammarella2024PrimarySR` against CPSR.m golden probabilities for primary surface-rupture occurrence across magnitude, magnitude-scaling relation, fault style, hypocentral-depth-distribution flag, seismogenic thickness, dip distribution, and truncation parameters.

## Tolerance and its justification

The coded tolerance is `abs=2e-3`, `rel=2e-3` in `test_cpsr_golden.py`. The local file describes this as strict parity against CPSR.m golden outputs; the tolerance is retained exactly as coded for comparison against the archived CSV generated from the external CPSR workflow.

## How to reproduce

```bash
pytest openquake/fdha/test/benchmark/mammarella_et_al_2024/test_cpsr_golden.py -q
```

## Status on HEAD

PASS on 2026-07-19 against the archived CPSR.m golden CSV (9 rows).
