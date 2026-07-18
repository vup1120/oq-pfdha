# Reference

## Cited references

- Valentini, A., Fukushima, Y., Contri, P., Gulerce, Z. - *The IAEA Exercise
  on Probabilistic Fault Displacement Hazard Assessment*, Earthquake Spectra
  (IAEA PFDHA exercise paper). Figures 4 (principal) and 6 (distributed) are
  the comparison targets.
- International Atomic Energy Agency (2025). *Benchmarking Current Practices
  in Probabilistic Fault Displacement Hazard Analysis for Nuclear
  Installations*. IAEA-TECDOC-2092, Vienna. DOI: 10.61092/iaea.74us-dn4n.
  Sections 5–7 define the cases; the paper's Figs 4/6 correspond to TECDOC
  Figs 14/18/21 and 15/19/22.

## Source of the reference curves

The vectors in `reference/*.csv` are transcribed verbatim (not digitized
from images) from the MATLAB plotting scripts provided by the exercise
coordinators - `Figure4 (1).m` and `Figure6 (2).m` (A. Valentini, IAEA) -
which contain the hazard curves supplied by each hazard-analyst team on an
18-point displacement grid (0.01–1000 cm); the M11 team curve uses its own
28-point grid. `reference/make_reference_csvs.py` embeds the transcription
and regenerates the CSVs.

## Source of the model inputs

Author-provided exercise input packages (site coordinates, fault traces,
geometry, magnitudes, recurrence):

- `PFDHA Benchmarking Study - Kumamoto - Input data.xlsx` (+ Futagawa
  shapefiles) - Uto trace vertices, single-segment M6.5 rate 18.9e-5/yr,
  floating-option rate 64.2e-5/yr, sites.
- `PFDHA Benchmarking Study - Le Teil - Input data.xlsx` (+ Cévennes
  shapefiles) - the `LRF2_Extended5.8km` trace is the 5.8-km principal
  source; sensitivity-case-2/base-case sheets give the principal
  (44.531 N, 4.670 E) and distributed (44.537 N, 4.667 E) sites, M5.5
  (workbook: 5.52), rate 4.6e-5/yr.
- `PFDHA Benchmarking Study - Norcia - Input data.xlsx` - MVFS trace, the
  Gaussian characteristic MFD (M6.4–7.0, Σ = 4.03281e-4/yr) and the
  principal/distributed sites.

These match Tables 9–14 of TECDOC-2092 and Table 1 of the exercise paper.

## Tolerances and their justification

Per-entry tolerances live in `manifest.py` next to the entry they govern.
They wrap the observed agreement (0.4–22%) with margin and reflect the
paper's own findings: team-to-team differences at the curve head are driven
by the surface-rupture-probability model choice, and the deep tail
(AFOE < 1e-8) by each team's unprescribed aleatory-truncation choice -
hence the `assert_dmax_m = 1.0` restriction on the L23/T13 tails. The
documented `post_factor` on the T13 entries stands in for the
not-yet-implemented Takao P2p term. M11 and V24 are compared qualitatively
(see README.md, "Known, documented deviations").

## Status

PASS - 16/16 quantitative entries on 2026-07-10
(`pytest openquake/fdha/test/benchmark/IAEA/ -m benchmark`: 16 passed).
