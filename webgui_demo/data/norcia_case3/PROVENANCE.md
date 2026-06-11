# Provenance of pre-computed demo data

These files are verbatim copies of outputs produced by the repository's own
Norcia Case 3 benchmark runner:

    python openquake/fdha/test/benchmark/norcia_case3_iaea/run_test.py

run against the tracked inputs in
`openquake/fdha/test/benchmark/norcia_case3_iaea/`
(`job_norcia_case3_iaea_curve.ini`, `job_norcia_case3_iaea_map.ini`,
`source_model_norcia_case3.xml`, and the two logic-tree XMLs).
The benchmark implements the IAEA TECDOC-2092 Case 3 distributed fault
displacement exercise (DOI 10.61092/iaea.74us-dn4n).

| File here | Original path (gitignored output) |
|---|---|
| `aggregate_hazard.csv` | `out_curve/aggregate_hazard.csv` |
| `hazard_curves/branch_000{0,1}.csv` | `out_curve/hazard_curves/` |
| `manifest.json` | `out_curve/manifest.json` |
| `map/displacement_map_mean.csv` | `out_map/aggregate/displacement_map_mean.csv` |
| `map/displacement_map_quantile-*.csv` | `out_map/aggregate/` |

The originals are gitignored (`.gitignore:22,24`); copies are committed here
so the demo GUI is self-contained. To regenerate, run the command above and
re-copy.
