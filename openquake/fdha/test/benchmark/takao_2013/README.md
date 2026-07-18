# Takao et al. (2013) - paper-figure reproduction benchmark

Validates the Takao et al. (2013) model implementations against the
paper's own worked examples (Figs. 10 and 11, cases (a) and (b)), using
both the anchors quoted in the text and the full curves extracted from
the paper PDF's vector graphics.

Reference: Takao, M., Tsuchiyama, J., Annaka, T., & Kurita, T. (2013).
Application of probabilistic fault displacement hazard analysis in Japan.
*Journal of Japan Association for Earthquake Engineering*, 13(1), 17-36.
https://doi.org/10.5610/jaee.13.17

## Reference data

`extract_reference.py` converts the figure page to SVG (`pdftocairo`),
parses every vector path (solid curves = stroked polylines, dashed/dotted
curves = sequences of small filled rectangles), calibrates each chart's
plot-area rectangle against its axis labels, and writes the data points to
`reference/fig10_points.csv`, `reference/fig11a_points.csv`,
`reference/fig11b_points.csv`. The PDF itself is not committed; the CSVs
are. Curve identity is resolved by nearest-curve assignment in the
comparison scripts (`figcompare.py`).

## Fig. 11 case (a) - principal faulting (`reproduce_fig11a.py`)

Site at the midpoint of a 22-km fault; six curves (R 3,000/30,000 yr ×
Mw 6.6/6.8/7.0) via the paper's Eq. (1), including the full P2p placement
machinery: Takemura (1998) rupture length (`log10 L = 0.5 Mw − 2.07`,
giving the paper's 21.38 km at Mw 6.8), the four segment-length-ratio
lines of Fig. 3, 1-km discretization, and half-open `[s, s+L)` site
coverage. P3p uses the `Takao2013PrimaryFD` building blocks with the
actual segment length as `srl` (AD method, per the paper).

| anchor (paper text) | paper | computed | agreement |
|---|---|---|---|
| P1p(Mw 6.8) | 0.784 | 0.7841 | 0.02% |
| P2p(Mw 6.8) | 0.751 | 0.7505 | 0.06% |
| placement counts (8 values, Mw 6.8) | 1/21, 10/11, 3/3, 1/1 and 1/21, 11/11, 3/3, 1/1 | identical | exact |
| nu(0.01 m \| Mw 6.6, R 30,000) | 1.2e-5 /yr | 1.189e-5 /yr | 0.9% |

Digitized-curve agreement (median / max abs. relative error per curve;
gray ink = R 3,000, black = R 30,000):

| curve | median | max | curve | median | max |
|---|---|---|---|---|---|
| Mw6.6 R3000 | 2.3% | 9.3% | Mw6.6 R30000 | 2.7% | 9.4% |
| Mw6.8 R3000 | 0.4% | 3.2% | Mw6.8 R30000 | 0.6% | 3.5% |
| Mw7.0 R3000 | 1.3% | 4.0% | Mw7.0 R30000 | 1.3% | 3.7% |

## Fig. 11 case (b) - distributed faulting (`reproduce_fig11b.py`)

Four faults (r = 5/5/10/10 km, R = 3,000/30,000/3,000/30,000 yr, Mw 6.8)
via Eq. (13) with `Takao2013PrimarySR` × `Takao2013SecondarySR` ×
`Takao2013SecondaryFD` (AD). `n_sigma = 5` reproduces the paper's
untruncated Eq. 12 integral (converged; the default 3 starves the tails
below AFOE ~1e-9).

| anchor | paper | computed | agreement |
|---|---|---|---|
| four-fault sum at 0.01 m | 7.0e-7 /yr | 6.99e-7 /yr | 0.2% |
| fault-2/fault-3 crossover (Sec. 4(2)) | present | present (~0.6 m) | ✓ |
| recurrence scaling | 10× | 10.000× | exact |

Digitized-curve agreement over the full 8-decade range (1e-10..1e-3):

| curve | median | max | curve | median | max |
|---|---|---|---|---|---|
| fault 1 | 1.0% | 5.2% | fault 4 | 0.9% | 3.2% |
| fault 2 | 0.8% | 4.1% | sum | 1.4% | 5.1% |
| fault 3 | 1.2% | 4.3% | | | |

## Fig. 10 - probability densities (`reproduce_fig10.py`)

The paper plots probability **mass on a 0.01-decade log grid** (the
digitized AD peak height 0.01108 equals the Eq. 10 lognormal log-density
1.108 × 0.01 exactly). The dotted/dashed curves correspond to the
**short-segment gamma of Eq. (8)** (a = 1.53, b = 0.58, the L < 10 km
branch) - not the long-segment Eq. (6) at x/L = 0 that the caption text
suggests: Eq. 6 at x/L = 0 would peak at 0.50 (D/AD) and 0.34 (D),
excluded by the digitized 0.89 and 0.58.

| curve | quantity | agreement (no fitted scale) |
|---|---|---|
| AD (Eq. 10 lognormal) | peak position / height / pointwise | 0.13% / 0.003% / max 0.8% |
| D (Eq. 12 convolution) | peak position / height / pointwise | 0.45% / 0.04% / max 1.0% |
| D/AD (Eq. 8 gamma) | peak position / shape (one fitted scale 0.127) | 1.0% / max 2.3% |

The D curve is the strongest single test: with **no adjustable
parameter** it validates the AD-grid-mass × gamma composition exactly as
`Takao2013PrimaryFD` implements it. The dotted D/AD curve carries an
unexplained uniform display scale in the paper (~1/7.9; its log-mass peak
would otherwise be 0.0108, comparable to AD's) - shape and position are
validated, the scale is fitted and reported.

## Files

- `extract_reference.py` - PDF → `reference/fig*_points.csv` (provenance)
- `figcompare.py` - nearest-curve matching statistics
- `reproduce_fig10.py` / `reproduce_fig11a.py` / `reproduce_fig11b.py` -
  computation, `Figures/*.png` overlays, `fig*_agreement.json`
- `test_paper_reproduction.py` - 23 pytest assertions on all of the above

## Not covered

The P2p machinery (Takemura length, Fig. 3 lines, placement counting)
lives in `reproduce_fig11a.py` as example-specific code - it is validated
here against the paper but is not yet a library model (the same P2p gap
documented in the IAEA benchmark, where the T13 principal chains use
constant post-factors).
