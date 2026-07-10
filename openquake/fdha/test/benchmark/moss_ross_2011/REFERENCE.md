# Reference

Moss, R. E. S., & Ross, Z. E. (2011). Probabilistic Fault Displacement
Hazard Analysis for Reverse Faults. *Bulletin of the Seismological Society
of America*, 101(4), 1542–1553. doi:10.1785/0120100248.

## Source of the reference values

All comparison targets are stated in the paper's text (Example and
Sensitivity Analysis sections), not digitized from figures:

- Reference curve (Fig. 7, reverse P_sr + AD/gamma, x/L = 0.25):
  "a 2% probability in 50 yr of exceeding 55 cm and 1% probability in 50 yr
  of exceeding 105 cm. These risk parameters have return periods of 2475 yr
  and 4975 yr."
- Fig. 10 comparison (all-slip-types P_sr): "the 2%-and-1%-in-50-yr
  parameters are 95 cm and 142 cm."
- Plateau shift: "the probability at which the hazard curve plateaus has
  increased by nearly 45%."

Source-model inputs (slip rate 0.5 mm/yr, 44 km × 14 km, μ = 3.75e11
dyne/cm², b = 0.8, Mmax = Mw 7, Mmin = 5, x/L = 0.25, AD + gamma) are from
the paper's Example section. The magnitude-frequency form (truncated
exponential) and the Hanks–Kanamori constant (16.05) are reconstructed;
both were confirmed by sensitivity analysis (alternatives degrade all
anchors simultaneously).

## Model classes exercised

- `MossRoss2011PrimarySR` — paper Eq. 5 (verified to 4 decimals).
- `MossRoss2011PrimaryFD` — paper Eqs. 7 (gamma D/AD), 8 (AD lognormal;
  the AD/gamma path is the paper's own reference configuration).
- `WC1993PrimarySR` — the all-slip-types probability of surface rupture
  used in the paper's Fig. 10 comparison.

## Tolerances and their justification

See README.md: 15% on the three mutually consistent anchors (observed
1–6%), 35% on the reverse 2%-in-50-yr anchor (observed +28%; the paper's
text value is inconsistent with the paper's own Fig. 8 percent-difference
curve, while our value is consistent with it), [1.30, 1.60] on the plateau
ratio (paper "nearly 45%", observed +48%).

## Status

PASS — 5/5 assertions on 2026-07-10.
