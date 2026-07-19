This benchmark compares pfdha model outputs and selected aleatory parameters against the external `fdhpy` reference implementation, with the citation and reference checkout documented in [REFERENCE.md](REFERENCE.md). Results are summarized in [RESULTS.md](RESULTS.md).

Run: `pytest openquake/fdha/test/benchmark/fdhpy_comparison/FDHI_Tests/tests -q`

Requires `fdhpy` installed (e.g. `pip install -e /path/to/fdhpy` at the commit pinned in REFERENCE.md); the tests skip if it is missing.

Status: 198/198 PASS on 2026-07-19 (see RESULTS.md).
