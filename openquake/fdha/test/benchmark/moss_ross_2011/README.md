# Moss & Ross (2011) benchmark - Los Osos example (Figs 7 and 10)

Validates the `MossRoss2011PrimarySR` / `MossRoss2011PrimaryFD`
implementation against the source paper's own worked example: the Los Osos
fault zone PFDHA (Moss & Ross 2011, BSSA; reference curve in their Fig. 7,
surface-rupture-distribution comparison in their Fig. 10).

```bash
PYTHONPATH=. python openquake/fdha/test/benchmark/moss_ross_2011/reproduce_mr2011_fig10.py
pytest openquake/fdha/test/benchmark/moss_ross_2011/ -m benchmark
```

Outputs: `Figures/mossross2011_fig10_reproduction.png`,
`fig10_agreement.json`.

## Setup (paper, Example section)

Slip rate 0.5 mm/yr, area 44 km x 14 km, shear modulus 3.75e11 dyne/cm²,
truncated-exponential magnitudes with b = 0.8 on Mw 5.0–7.0 (moment-balanced
via log M0 = 1.5 Mw + 16.05 → α(M≥5) = 1.15e-2/yr), x/L = 0.25, AD method
with the gamma D/AD distribution (paper Eq. 7). The paper's hazard integral
(Eqs. 2–3) evaluates every event at the fixed site position - no rupture
floating - so the benchmark computes the direct integral with our model
classes rather than an engine job.

Prior coefficient verification (see `IAEA/le_teil/diagnose_M11.py`):
`MossRoss2011PrimarySR` matches Eq. 5 (a = 7.30, b = −1.03) to 4 decimals
and the FD class carries Eqs. 7–9 verbatim.

## Results (2026-07-10)

| Anchor | paper | computed | ratio |
| --- | --- | --- | --- |
| reverse, 2% in 50 yr (1/2475) | 55 cm | 70.2 cm | 1.28 |
| reverse, 1% in 50 yr (1/4975) | 105 cm | 102.5 cm | 0.98 |
| all-slip-types, 2% in 50 yr | 95 cm | 96.0 cm | 1.01 |
| all-slip-types, 1% in 50 yr | 142 cm | 133.3 cm | 0.94 |
| plateau shift (all/reverse) | "nearly 45%" | +48% | - |

Three of four anchors agree within 6% and the plateau shift matches; the
reverse 2%-in-50-yr anchor is the lone outlier (+28%). Three observations
establish that the paper's *text value* (55 cm) is the inconsistent
quantity, not our curve:

1. **The anchors are over-determined and 55 cm is the one that cannot
   coexist with the other three.** Both curves share the same activity rate
   α, so any α revision rescales them together. Forcing our reverse curve
   through (55 cm, 1/2475 yr) requires α × 0.70 - but that same factor
   drags the other three anchors to ratios 0.80–0.82, destroying the
   otherwise 1–6% agreement. No common activity rate satisfies all four
   published values simultaneously.
2. The paper's own Fig. 8 shows the percent difference between the
   all-slip-types and reverse curves at 0.5–1 m displacement to be roughly
   30–45%. Our curves give +37% at that hazard level - consistent - whereas
   the paper's text anchors imply 95/55 = +73%, inconsistent with its own
   Fig. 8.
3. The discrepancy is geometrically amplified where it sits: the 2%-in-50-yr
   level crosses the reverse curve on its shoulder (local log-log slope
   ≈ −1.6, vs −2.1 at the 1%-in-50-yr level), so the 28% horizontal offset
   corresponds to only a ~44% vertical one - about one contour of the
   paper's unseeded ~5000-sample Monte Carlo noise plus curve-reading
   error on a log axis.

An all-slip-types-only variant of the figure (no reverse curve/anchors) is
written alongside as `Figures/mossross2011_fig10_reproduction_allslip.png`;
both its paper anchors fall on the computed curve.

Sensitivity checks (documented in the session that built this benchmark):
the result is robust to the aleatory truncation (n_sigma 3 vs 5) and
degrades under a characteristic (Youngs & Coppersmith) MFD or a 16.1
Hanks–Kanamori constant - confirming the truncated-exponential / 16.05
reconstruction.

Tolerances in `test_fig10_reproduction.py` encode exactly this picture:
15% on the three consistent anchors, 35% on the documented outlier, and
[1.30, 1.60] on the plateau ratio.
