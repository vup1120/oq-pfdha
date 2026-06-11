# Configuration guide for Nishizaka et al. (2026) — Far-Field Surface Rupture

Reference: Nishizaka, N., Onishi, K., Ikeda, M., Si, H., Yamamoto, K., & Tsuji, T. (2026). Characteristics of far-field surface ruptures caused by two recent strike-slip earthquakes: Insights into fault displacement prediction. *Seismological Research Letters*. https://doi.org/10.1785/0220250293

The Nishizaka et al. (2026) models revise the conventional Japanese far-field (distributed) equations for **strike-slip** faults using the densely mapped surface ruptures of the 2016 Kumamoto and 2019 Ridgecrest earthquakes. Two features distinguish them from the other models in this library:

1. Surface ruptures are **not classified as principal or distributed** — the occurrence probability `P2s` covers all surface ruptures in a unit cell.
2. A second distance predictor is introduced: `r2`, the shortest surface distance from **mapped, pre-existing active faults**, capturing the strong concentration of far-field rupture near reactivated structures.

Both models are magnitude-independent (magnitude enters the displacement model only through the source-fault average displacement, EAD).

## Model classes

| Class | Category | Equation |
| :--- | :--- | :--- |
| `Nishizaka2026SecondarySR` | Secondary surface rupture | `P2s = logistic(c1 + c2·ln(r1+c3) + c4·ln(r2+c5))` (Eq. 2, conventional terms; Table 1, 250 m cells) |
| `Nishizaka2026SecondaryFD` | Secondary surface displacement | `D/EAD ~ Gamma(shape a, scale (c6/c8)·e^(c7·r1))` (Eqs. 3–4; Table 2) |

## `Nishizaka2026SecondarySR` parameters

| Name | Type | Units | Default | Allowed | Required? | Description |
| --- | --- | --- | --- | --- | --- | --- |
| `r2` | float | km | – | ≥ 0 | Yes | Shortest surface distance from mapped, pre-existing active faults to the site. |
| `region` | string | – | `kumamoto` | `kumamoto`, `ridgecrest` | No | Event-specific regression (epistemic alternatives). |
| `side` | string | – | `plain` | `plain`, `mountain` | No | Plain or mountain side of the earthquake source fault. |
| `cell_size` | int | m | `250` | `250` | No | Unit cell size. Only the main-text Table 1 (250 m) coefficients are bundled; the supplemental Table S1 sizes are not included. |

The distance `r` (named `r1` in the paper) is supplied by the calculator. Note the paper measures `r1` as the 3D distance from the **subsurface earthquake source fault**; the calculator's rupture-trace distance is used as an approximation.

## `Nishizaka2026SecondaryFD` parameters

| Name | Type | Units | Default | Allowed | Required? | Description |
| --- | --- | --- | --- | --- | --- | --- |
| `ead` | float | m | – | > 0 | Yes | Average displacement of the earthquake source faults (e.g. from a fault model or empirical scaling; 1.9 m for the 2016 Kumamoto earthquake after Asano & Iwata, 2021). |
| `dataset` | string | – | – | `proximal`, `nonproximal`, `all` | One of `dataset`/`r2` | Table 2 regression subset: proximal (`r2` ≤ 1 km) or non-proximal (`r2` > 1 km) to mapped, pre-existing active faults. |
| `r2` | float | km | – | ≥ 0 | One of `dataset`/`r2` | If given instead of `dataset`, the subset is selected automatically (threshold 1 km). |

Displacement coefficients exist only for the **2016 Kumamoto** dataset; the Ridgecrest far-field displacement data were too sparse to regress (see the paper).

## Example logic-tree branches

```xml
<logicTreeBranch branchID="SS_SSR_NISHIZAKA2026">
  <uncertaintyModel>
    [Nishizaka2026SecondarySR]
    r2 = 0.5
    region = "kumamoto"
    side = "plain"
  </uncertaintyModel>
  <uncertaintyWeight>1.0</uncertaintyWeight>
</logicTreeBranch>

<logicTreeBranch branchID="SS_SFD_NISHIZAKA2026">
  <uncertaintyModel>
    [Nishizaka2026SecondaryFD]
    ead = 1.9
    r2 = 0.5
  </uncertaintyModel>
  <uncertaintyWeight>1.0</uncertaintyWeight>
</logicTreeBranch>
```

## Notes and cautions

- **Strike-slip only**, and calibrated on two events — the authors stress that the coefficients are event/tectonic-setting specific and recommend regional regressions for application elsewhere.
- **Do not combine** `Nishizaka2026SecondarySR` with a separate principal-rupture model near the source: `P2s` already includes all ruptures, so pairing it with a principal model can double-count near-fault hazard. It is intended for **far-field** sites (≳ 3 km).
- The implementation is validated against the values stated in the paper's Results section (P2s at r1 = 10 km for both events/sides, and the 90th-percentile normalized displacements at r1 = 5 and 10 km); see `openquake/fdha/test/unit/test_nishizaka2026.py`.
- The revised equations in the paper (their full Eq. 2) also include the `c4·ln(r2+c5)` term fitted per side; the implemented form is exactly the published Table 1 parameterization.

## References

- Nishizaka, N., Onishi, K., Ikeda, M., Si, H., Yamamoto, K., & Tsuji, T. (2026). Characteristics of far-field surface ruptures caused by two recent strike-slip earthquakes: Insights into fault displacement prediction. *Seismological Research Letters*. https://doi.org/10.1785/0220250293
- Takao, M., Ueta, K., Annaka, T., Kurita, T., Nakase, H., Kyoya, T., & Kato, J. (2014). Reliability improvement of probabilistic fault displacement hazard analysis. *Journal of Japan Association for Earthquake Engineering*, 14(2), 16–36. https://doi.org/10.5610/jaee.14.2_16
- Takao, M., Tani, T., Oshima, T., Annaka, T., & Kurita, T. (2016). Maximum likelihood estimation of the parameters regarding displacement evaluation of distributed fault in PFDHA. *Journal of Japan Association for Earthquake Engineering*, 16, 96–101. https://doi.org/10.5610/jaee.16.2_96
