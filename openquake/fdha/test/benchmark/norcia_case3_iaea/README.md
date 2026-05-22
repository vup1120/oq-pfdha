This benchmark implements the IAEA TECDOC-2092 Norcia Case 3 distributed fault displacement exercise for the SRL paper application section.

Run:

```bash
python openquake/fdha/test/benchmark/norcia_case3_iaea/run_test.py
python openquake/fdha/test/benchmark/norcia_case3_iaea/plot_results.py
```

Main outputs:

- `out_curve/aggregate_hazard.csv`
- `out_map/aggregate/displacement_map_mean.csv`
- `Figures/norcia_case3_curve_vs_iaea.pdf`
- `Figures/norcia_case3_map.pdf`

The curve job uses one multi-site INI with the three application sites in this order: PF (`13.278 42.767`), MS (`13.188 42.749`), and SL (`13.212 42.853`).

Important implementation note: this repository exposes `Visini2025SecondarySR` and `Visini2025SecondaryFD`, but not `Visini2025PrimarySR` or `Visini2025PrimaryFD`. The V24 end-branch therefore uses the same Youngs2003 principal SR/FD slots as the Y03 branch and switches only the distributed SR/FD slots to Visini2025. This keeps the branch runnable without inventing unavailable primary Visini classes.

Faulting style is inferred from the source rake (`-90`, normal) during context construction. A `style` entry can still be supplied inside an FDHA logic-tree `uncertaintyModel` as an explicit model-parameter override, but this single-style Norcia benchmark leaves it to the rake-derived classification.

The source XML uses the fault labels and recurrence values from the author-provided Norcia workbook. The source IDs/names are `MVFS` and `NFS`; their incremental MFDs sum to `4.03281e-4 /yr` and `5.407072e-3 /yr`, respectively.
