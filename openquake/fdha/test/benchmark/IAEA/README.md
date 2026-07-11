# IAEA PFDHA exercise benchmark

Reproduces the three case studies of the IAEA probabilistic fault
displacement hazard (PFDHA) exercise — **Kumamoto** (strike-slip), **Le
Teil** (reverse) and **Norcia** (normal) — and compares the oq-pfdha hazard
curves with the curves supplied by the exercise's hazard-analyst teams
(Valentini et al., IAEA exercise paper, Figs 4 and 6; the same curves appear
in IAEA TECDOC-2092 as Figs 14/18/21 and 15/19/22). See
[REFERENCE.md](REFERENCE.md) for provenance.

## Quick start

```bash
# run everything, write computed/ CSVs + comparison_summary.json + summary table
PYTHONPATH=. python openquake/fdha/test/benchmark/IAEA/run_all.py           # all jobs
PYTHONPATH=. python openquake/fdha/test/benchmark/IAEA/run_all.py norcia    # filter

# overlay figures (requires run_all.py output)
PYTHONPATH=. python openquake/fdha/test/benchmark/IAEA/plot_comparison.py

# the same comparisons as assertions
pytest openquake/fdha/test/benchmark/IAEA/ -m benchmark -v
```

Outputs: `computed/<case>_<job>.csv`, `comparison_summary.json`,
`Figures/iaea_fig4_principal_comparison.png`,
`Figures/iaea_fig6_distributed_comparison.png`.

## Layout

- `reference/` — the teams' published hazard curves, transcribed verbatim
  from the exercise coordinators' MATLAB plotting scripts
  (`make_reference_csvs.py` regenerates the CSVs and documents the source).
- `kumamoto/`, `le_teil/`, `norcia/` — one job INI + FDHA logic tree per
  hazard-analyst model chain, plus the source models built from the
  author-provided input workbooks and shapefiles.
- `manifest.py` — job ↔ reference-curve mapping, documented post-factors and
  assertion tolerances (single source of truth for `run_all.py` and the
  pytest test).

## Exercise definition used here

Source parameters come from the author-provided workbooks
(`PFDHA Benchmarking Study - <case> - Input data.xlsx`) and shapefiles:

| Case | Source | Magnitude / rate | Principal site | Distributed site |
|---|---|---|---|---|
| Kumamoto | Uto segment (22 km, dip 60 NW, 14 km) | M6.5, 18.9e-5/yr | on-fault vertex, l/L≈0.23 | r = 5.2 km (HW) |
| Kumamoto (floating) | Futagawa zone (78 km, dip 64.6 NW) | M6.5, 64.2e-5/yr floated | same site | — |
| Le Teil | La Rouvière `LRF2_Extended5.8km` (5.8 km, dip 45 SE, 4 km) | M5.5, 4.6e-5/yr | on-fault, l/L≈0.46 | r = 0.6 km (FW) |
| Norcia | Mount Vettore (34 km, dip 47 SW, 11 km) | Gaussian characteristic MFD M6.4–7.0, Σ = 4.03e-4/yr | near SE end, l/L≈0.05 | r = 7.6 km (HW) |

The floating option uses a `simpleFaultSource` with `PeerMSR` and aspect
ratio 1.3 so the M6.5 rupture spans the full seismogenic width and floats
along strike only, matching the exercise convention (a width-limited scaling
relation would also float ruptures down-dip and dilute the
surface-rupturing rate by ~3x).

## Model chains (one job per exercise team)

| Team | oq-pfdha chain | Notes |
|---|---|---|
| P11 | `WC1993PrimarySR` + `Petersen2011PrimaryFD (bilinear)`; distributed: `Petersen2011SecondarySR/FD (cell 100 m)` | paper: WC93 gives P_sr≈0.70 at M6.5 |
| C24 | `WC1993PrimarySR` + `Chiou2025PrimaryFD (model7)` | same chain as `benchmark/valentini_et_al_2025` |
| K24 | `FixedPrimarySR (1.0)` + `Kuehn2024PrimaryFD`, full-posterior **mean** (the default ensemble reduction) | paper: K24 team assumed P_sr = 1 |
| T13 | `Takao2013PrimarySR` + `Takao2013PrimaryFD (AD, n_sigma = 5)` | published curves include the P2p site-rupture factor (~0.45–0.48), not yet a library model → applied as a documented `post_factor` in `manifest.py` |
| L23 | `FixedPrimarySR (case-specific)` + `Lavrentiadis2023PrimaryFD (disp_prnc_prime, include_zero_slip)` | the team computed a case-specific P_sr by floating the rupture over the dip surface; adopted values 0.85 / 0.65 / 1.0 (Kumamoto / Le Teil / Norcia) back-calculated from their published curves and consistent with the rupture-width / fault-width ratio |
| M11 | `Moss2013PrimarySR (stiff, vs30 760)` + `MossRoss2011PrimaryFD (AD)` | **qualitative only** — see below |
| Y03 | `Youngs2003PrimarySR/FD (normal, AD)`; distributed: `Youngs2003SecondarySR (v3)` + `Youngs2003SecondaryFD (85th)` | |
| V24 | gate `FixedPrimarySR (1.0)` + `Visini2025SecondarySR (100 m)` + `Visini2025SecondaryFD (WC1994)` | **qualitative only** — see below |

## Result summary (2026-07-10)

All quantitative entries pass; max |relative error| against the published
team curves over the asserted range:

| Entry | max rel. err | Entry | max rel. err |
|---|---|---|---|
| Kumamoto P11 | 4.3% | Le Teil K24 | 8.5% |
| Kumamoto C24 | 11.0% | Le Teil T13 (≤1 m) | 10.7% |
| Kumamoto K24 | 3.2% | Le Teil L23 (≤1 m) | 5.3% |
| Kumamoto T13 | 7.2% | Norcia Y03 | 16.8% |
| Kumamoto L23 (≤1 m) | 5.6% | Norcia K24 | 19.7% |
| Kumamoto distributed P11 | **0.4%** | Norcia L23 | 13.3% |
| Kumamoto floating K24/T13/L23 | 17.4 / 5.9 / 22.3% | Norcia distributed Y03 (≤1 m) | 19.1% |

Where an entry is asserted only up to 1 m displacement, the excluded tail
sits at AFOE < 1e-8–1e-12, where the paper itself notes the curves are
controlled by each team's (unprescribed) aleatory-truncation choice.

## Known, documented deviations

- **T13 P2p factor.** The Takao et al. (2013) conditional probability of
  principal rupture at the site (P2p, TECDOC-2092 §3.2.2) is not implemented
  as a library model; the published T13 curves include it. It is applied as
  a constant `post_factor` (0.4801 Kumamoto, 0.4748 Le Teil, back-calculated
  from the published first point; the paper quotes "a factor of 0.45").
- **T13 distributed curves** (Figs 6a/6b) need the Takao et al. (2014/2016)
  distributed displacement model, which is not yet implemented
  (`secondary_surf_displ` has no Takao class). They are plotted from the
  reference only.
- **M11 (Le Teil).** Diagnosed against the Moss & Ross (2011) paper (see
  `le_teil/diagnose_M11.py` and `Figures/iaea_leteil_M11_diagnosis.png`):
  our implementation reproduces the paper's Eq. 5 P_sr to 4 decimals and
  carries Eqs. 7–9 verbatim, but the published M11 curve **cannot be the
  single-M5.5 scenario** — normalized, it is a lognormal with median ≈1.4 m
  (⇔ MD median at **M ≈ 6.6** via the paper's own Eq. 9), whereas M5.5 caps
  the median at AD = 0.37 m / MD = 0.41 m and gives P(D > 10 m) ≈ 5e-7
  against the published 1e-2 of head. Running our chain on the workbook's
  10-km-seismogenic-thickness branch (M6.23/6.48/6.73, w = 0.2, the team's
  own Eq. 5 P_sr) lands within ~1.7x of the published head and crosses the
  curve near 1 m — evidently the team's submission is dominated by the
  large-magnitude thickness alternative (consistent with Moss & Ross's own
  magnitude-integration methodology, cf. their Los Osos example), not the
  single-M5.5 case the other Fig. 4c curves represent. The paper's text
  ("Moss et al. 2013, 10% at M5.5") is itself inconsistent with the
  published head (0.0122 x rate). The entry stays qualitative.
- **V24 (Le Teil / Norcia).** The V24 team's principal-P_sr assumption is
  not stated in the paper; the jobs use P_sr = 1 (the K24-style neutral
  assumption). Gate trials (head of curve, computed/paper): Norcia — Takao
  gate 0.14, P_sr = 1 0.19 (still ~5x low); Le Teil — Takao gate 0.10,
  P_sr = 1 17 (overshoots). No single P_sr convention reconciles both
  cases, and the mid-curve slopes differ as well, consistent with the
  paper's Fig-6 V24 curves having been computed with the earlier (2024,
  under-review) revision of the Visini et al. model. Our
  `Visini2025Secondary*` implementation follows the final published model
  (validated directly in `benchmark/visini_et_al_2025`) and with the Takao
  gate lands on the TECDOC-2092 versions of the V24 curves (Figs 19/22).
- **Norcia Y03** sits 6–17% above the published curve (assertion at 20%).
  The head offset (~6%) corresponds to the difference between our Youngs
  Great-Basin P_sr and the team's effective value; the exercise paper
  itself attributes team-to-team spread of this size to the P_sr choice.
