# LA23 reference implementation (authors' code)

`LA23_flt_rup_model.py` is the Lavrentiadis & Abrahamson (2023) implementation published by the authors,
copied unmodified from <https://github.com/NHR3-UCLA/LA23_PFDHA_model> at commit `c3ba1136` (2023-10-04).
It is distributed under the MIT licence in `LICENSE` (© 2023 Natural Hazards Risk and Resiliency Research Center).

`tests/test_lavrentiadis2023_parity.py` imports it directly and compares the pfdha classes with it. It covers the full-rupture
aggregate and principal terms (`LavrentiadisAbrahamson2023SlipProfile`), their percentiles without the zero terms
(`…SlipProfilePrc`), and the average (`…AvgDisp`) and maximum (`…MaxDisp`) displacements. It has no
individual-segment principal model, CDF, mean or zero-displacement mixture. For those, the reference is fdhpy
1.0.3 (`fixtures/la23`) or the paper.
