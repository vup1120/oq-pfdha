# oq-pfdha validation report - IAEA PFDHA benchmarking exercise

**Date:** 2026-07-11 · **Scope:** the three IAEA exercise case studies
(Kumamoto strike-slip, Le Teil reverse, Norcia normal) and the
model-level validation of the FDHA chains they use.

The validation is organized in three tiers: (1) each displacement /
rupture-probability model is validated against its **source
publication**; (2) full engine model chains are validated against the
**published team curves** of the IAEA exercise (TECDOC-2092 / Valentini
et al. Figs 4 and 6); (3) the **complete author-workbook definitions**
- epistemic logic trees and all sensitivity cases - are implemented and
regression-tested, so the benchmark covers everything the exercise
defines, not only the configurations with published curves.

---

## Tier 1 - model-level validation against source publications

| Suite | What is validated | Result |
|---|---|---|
| `benchmark/takao_2013` | Takao et al. (2013) chain against the paper's own Figs 10/11 (curves extracted from the PDF vector graphics) and every text anchor | Fig 11(b): all 5 curves ~1% median / ≤5.2% max over 8 decades; Fig 11(a) incl. full P2p machinery: P2p 0.7505 vs 0.751, all placement counts exact, 6 curves 0.4–2.7% median; Fig 10: AD and the Eq. 12 convolution match with **no free parameter** (peak 0.13–0.45% / height 0.003–0.04%). 23/23 tests |
| `benchmark/moss_ross_2011` | Moss & Ross (2011) reverse chain vs their Los Osos example (Figs 7/10) | anchors within 1–6% (one documented outlier where the paper's text contradicts its own figure). 5/5 tests |
| `benchmark/visini_et_al_2025` | Visini et al. (2025) distributed models vs the publication | direct validation suite |
| `benchmark/valentini_et_al_2025` | Chiou (2025) chain vs the published Kumamoto application | ≤12% |
| `benchmark/mammarella_et_al_2024`, `benchmark/norcia_*`, `benchmark/fdhpy_comparison` | further per-publication and cross-implementation checks | see suite READMEs |

Unit level: 400+ tests including exact anchors (e.g.
`Takao2013SecondaryFD` places exactly 10% exceedance on the paper's
Eq. 15/16 regression levels).

## Tier 2 - IAEA exercise: published-curve comparisons (18 entries)

Each entry = one team model chain, run through the engine as a logic-tree
job and asserted against the team's published curve (`manifest.py`,
`run_all.py`, `test_iaea_benchmark.py`). Max |relative error| over the
asserted range:

| Entry | err | Entry | err |
|---|---|---|---|
| Kumamoto P11 | 4.3% | Le Teil K24 | 8.5% |
| Kumamoto C24 | 11.0% | Le Teil T13 (≤1 m) | 10.7% |
| Kumamoto K24 | 3.2% | Le Teil L23 (≤1 m) | 5.3% |
| Kumamoto T13 | 7.2% | Le Teil distributed T13 (≤1 m) | 36.7% (40% limit)¹ |
| Kumamoto L23 (≤1 m) | 5.6% | Norcia Y03 | 16.8% |
| Kumamoto distributed P11 | **0.4%** | Norcia K24 | 19.7% |
| Kumamoto distributed T13 | 32.5% (40% limit)¹ | Norcia L23 | 13.3% |
| Kumamoto floating K24/T13/L23 | 17.4 / 5.9 / 22.3% | Norcia distributed Y03 (≤1 m) | 19.1% |

¹ The T13 distributed heads match at 0.1% / 9%; the mid-range shape
difference is attributable to the team (the library model reproduces the
Takao 2013 paper itself at ~1%, Tier 1), and the team's aleatory
integration is undocumented in TECDOC-2092.

Qualitative entries (documented reasons, README "Known deviations"):
Le Teil **M11** (our stiff `Moss2013PrimarySR` = 9.9% reproduces the
paper's stated "10% at Mw5.5" and the expected plateau
`4.6e-5 x 0.099 = 4.6e-6`; the paper's *plotted* plateau implies an
effective P_sr of only 1.2%, inconsistent with its own text - see
`le_teil/diagnose_M11_psr.py`), Le Teil / Norcia **V24** (paper curves
computed with the earlier, under-review model revision), and the
base-case-tree entries below.

## Tier 3 - complete workbook coverage

**Epistemic base-case trees** (all published curves correspond to the
single-scenario configurations; the trees quantify the full model):

| Case | Tree | Implementation | Tree vs single-scenario curves |
|---|---|---|---|
| Kumamoto | 4 coexisting rupture sources × independent mag/rate branches (w .2/.6/.2) | `kumamoto/make_basecase.py` (mean-collapse: 3-bin MFDs, exact for the mean) + P11/T13 jobs | head ×1.8–2.1, 10-m tail ×28–32 (Mw 7.0–7.4 multi-segment branches) |
| Le Teil | 36 end branches: thickness × rupture length × magnitude × slip rate | `le_teil/make_basecase.py` (12 explicit smlt branches, slip level collapsed) + M11 job | M11 tail moves toward the published curve (3/5 m ratio 0.18/0.02 → 0.37/0.14); head ×14 - confirms the M11 diagnosis |
| Norcia | none (single MVFS source, Gaussian MFD) | already the standard jobs | - |

**Sensitivity cases** (no published curves → demonstration jobs with
regression snapshots + physical-consistency tests,
`make_sensitivity_cases.py` / `sensitivity.py` /
`test_sensitivity_cases.py`, 11 tests):

| Job | Configuration | Head (AFOE) | Consistency checks |
|---|---|---|---|
| Kumamoto sens1 (T13) | Suizenji Mw 5.8, r = 0.6 km | 1.75e-8 | Petersen chain correctly rejects M5.8 (outside M6–8 calibration) → Takao chain used |
| Kumamoto sens3 (T13) | Suizenji, on-fault | 6.14e-6 | head = rate × P1p(5.8) to 5% |
| Kumamoto sens4 (P11) | 4 sources, r = 10 km | 2.14e-7 | < base-case head at r = 5.22 km |
| Le Teil sens1 (T13) | 4 fault-source combos (w .25), r = 0.6 km | 6.67e-9 | 5-branch smlt from the Cévennes shapefile |
| Le Teil sens4 (T13) | LRF_2 + MRF_3 + PCF_2, r = 0.6 km | 1.16e-9 | 3-source sum |
| Norcia sens2 (Y03) | MVFS, r = 2.4 km | 9.45e-6 | > base-case site (r = 7.6 km) |
| Norcia sens3 (Y03) | MVFS + NFS (GR MFD 5.5–6.8), 3 km from NFS | 1.05e-5 | > MVFS alone |

Overlay: `Figures/iaea_sensitivity_cases.png` (regenerate with
`plot_sensitivity.py` after `run_sensitivity.py`).

## Test inventory

```
pytest openquake/fdha/test/benchmark/IAEA/        # 18 published-curve + 11 sensitivity
pytest openquake/fdha/test/benchmark/takao_2013/  # 23 paper-reproduction
pytest openquake/fdha/test/unit                   # 414 unit tests
```

All pass as of this report. Regenerating inputs from the author data:
`make_basecase.py` (per case) and `make_sensitivity_cases.py` (requires
the workbooks + shapefiles, not committed); regenerating snapshots:
`run_sensitivity.py --update-snapshots`.

## Known limitations

1. **Takao P2p** is validated as example code (`benchmark/takao_2013`,
   0.06% vs the paper) but is not yet a library model; the IAEA T13
   principal entries carry documented constant post-factors.
2. **M11** and **V24** published curves cannot be reproduced from the
   exercise definitions (team-specific configurations / earlier model
   vintage); both are diagnosed and documented rather than asserted.
3. Sensitivity-case jobs are regression-tested, not
   reference-validated - the exercise published no curves for them.
