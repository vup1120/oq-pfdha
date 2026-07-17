# Petersen et al. (2011) benchmark — Fig. 9c rupture-location weight

Validates the tool's rupture-location weight kernel
(`openquake.fdha.calc.location_weight.location_weight`, σ > 0 path) against
Petersen et al. (2011, BSSA 101(2), 805–825), Fig. 9c (paper p. 820): the
across-strike displacement-hazard profile for the four mapping-accuracy
classes (Accurate / Approximate / Concealed / Inferred).

```bash
cd /path/to/oq-pfdha
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=$PWD python -m pytest openquake/fdha/test/benchmark/petersen_2011 -v
```

Runs in ~3 s. Pure numpy/scipy — no full hazard-job runs, no digitization
code executed at test time (the digitized arrays are frozen data, see below).

## What this tests

The worked example from the paper (p. 819): characteristic M 7.0 every
140 yr, `P(sr != 0 | m)` from eq. 5 (`WC1993PrimarySR`), bilinear principal
displacement at the mid-rupture branch l/L = 0.5 ≥ 0.3
(`Petersen2011PrimaryFD_bilinear`, eq. 8), distributed occurrence + eq. 18
displacement at 25-m cells (`Petersen2011SecondarySR`,
`Petersen2011SecondaryFD`, near-field floored per the tool's D7 constant
`NEAR_FIELD_FLOOR_KM`), and the rupture-location weight `W_p(r)` from the
tool's real kernel (`location_weight`, σ > 0 path — pinned, pure Gaussian,
truncated beyond ±2σ). The four class curves are built by solving, for each
across-strike distance r, the displacement `D0` at which

```
lambda_total(D0, r) = alpha * P_sr * [W_p(r) * P(D_principal > D0)
                                      + P(d != 0 | r) * P(D_distributed > D0 | r)]
```

equals λ* = −ln(0.9)/50 (10% in 50 yr), per the additive kernel spec (D1,
σ > 0 → summed, not complementary) in
`docs/design/rupture_location_uncertainty.md`. All model calls go through
the tool's real classes; nothing is re-implemented locally except the
bisection solver and the RMS bookkeeping.

## Provenance: the digitized figure

`fig9c_digitized.npz` (copied from `docs/design/figures/fig9c_digitized.npz`
so this benchmark is self-contained; that source tree is untracked) holds
the four printed curves digitized at 300 dpi from the paper's p. 820 panel
(script of origin: `digitize_fig9c.py`, session scratchpad, not ported here —
it is a one-off image-processing tool, not part of the runtime benchmark).
Each array is `(2, N)`: row 0 = across-strike distance r (m, signed), row 1 =
displacement D₀ (cm). A verification overlay
(`digitize_check.png`, not shipped) confirmed the extracted dots sit on the
printed curves before freezing the npz.

## The documented text-vs-figure inconsistency

Inverting the digitized curves through the paper's own published equations
(scripts `invert_wp.py`, `roundtrip.py`, session scratchpad) shows that the
**printed** Fig. 9c curves encode an implied rupture-location weight

```
W_p_figure(r) ≈ 0.90 * exp(-r^2 / (2 * (1.65 * sigma_table)^2))
```

— i.e. an on-trace pin of ~0.90 (not 1) and an effective Gaussian width
~1.65× the Table 2/3 two-sided σ. **Neither the 0.90 pin nor the 1.65 width
factor is derivable from anything stated in the paper.** This is the third
documented text-vs-figure inconsistency found in this codebase's Petersen
(2011) work (after the on-trace pinning ambiguity and the Fig. 10a absolute
levels — see `docs/design/rupture_location_uncertainty.md`).

The 0.90 factor also closes a long-standing peak-value gap: the exact
closed-form peak of the **stated** method (α = 1/140, eq. 5 `P_sr` = 0.8654,
λ* = −ln(0.9)/50, eq. 8 μ = 4.4644 / σ_ln = 0.9624 at X/L = 0.5) is 128.9 cm
(`peak_audit.py`); 0.90 × that ≈ 124 cm ≈ the printed "123 cm" figure text.

## Policy

The tool implements the paper's **stated** method — the only citable
anchor: a pure Gaussian at the Table 2/3 σ, pinned to 1 on the trace,
truncated at ±2σ. This benchmark compares that stated method against the
digitized figure **on an explicitly documented-offset basis** — it does not
tune anything to hit the printed curves. Concretely:

- the stated method reproduces the digitized curves' *shape* (Gaussian
  widening with σ) but sits high/narrow relative to the print (rms ≈
  37–49 cm across the four classes);
- applying the reverse-engineered figure transform (σ → 1.65σ, pin → 0.90)
  reproduces the print far better (rms ≈ 7.6–11.8 cm) — closing the loop on
  the inversion, without being part of the tool's default behaviour.

**If you want the tool to reproduce the printed Fig. 9c rather than the
paper's stated method, feed `σ × 1.65` as your `r_sigma_km` branch value**
(the tails will then match the print; the on-trace peak will still sit ~10%
*above* the printed value, because `W_p(0) = 1` for every σ — the 0.90 pin
is not exposed as a user-facing option). See `docs/UserManual_Enhanced/09-Uncertainty.md`
(`fdhaCalcRSigma` section) for the full provenance table and the
strike-slip-provenance caveat.

## Assertions (`test_petersen2011_fig9c.py`)

1. **Pinned equal peaks** — the four class curves peak at the same `D0` at
   r = 0 (relative spread < 1e-9: `W_p(0) = 1` for every σ), and that shared
   peak (134.0 cm observed) is within ±10% of the analytic 128.9 cm peak.
2. **Tails widen monotonically with σ**; each class's `W_p(r)` is exactly 0
   strictly beyond 2σ (kernel contract — the interval `|r| ≤ 2σ` is closed).
3. **Documented-offset lock** — per class, `rms(transformed vs digitized) <
   rms(stated vs digitized)`, and the transformed rms stays below 15 cm
   (observed 7.6–11.8 cm vs the stated method's 36.7–49.2 cm).
4. **Peak ratio** digitized/stated ≈ 0.90, asserted in `[0.85, 0.95]`
   (observed 0.905).

## Measured numbers (2026-07-17)

| class | σ (m) | stated peak (cm) | rms stated (cm) | rms transformed (cm) |
| --- | --- | --- | --- | --- |
| Accurate | 26.89 | 134.0 | 36.7 | 7.6 |
| Approximate | 43.82 | 134.0 | 40.2 | 9.2 |
| Concealed | 65.52 | 134.0 | 49.2 | 11.8 |
| Inferred | 72.69 | 134.0 | 43.6 | 10.8 |

Peak ratio (digitized mean / stated mean) = 121.2 / 134.0 = 0.905.
