# Configuration guide for Mammarella et al. (2024) — Principal Surface Rupture Probability

Reference: Mammarella, L., Visini, F., Boncio, P., Baize, S., Scotti, O., Beauval, C., Pace, B., & Thompson, S. (2024). Conditional probability of surface rupture: A numerical approach for principal faulting. *Earthquake Spectra*. https://doi.org/10.1177/87552930241293570

The Mammarella et al. (2024) model computes the conditional probability of principal surface rupture (CPSR) numerically, by integrating over probability distributions of rupture width (from a magnitude scaling relation), fault dip, hypocentral depth ratio, and seismogenic thickness, instead of using an empirical magnitude-only logistic regression.

## Model selection

Select `MammarellaEtAl2024PrimarySR` (or its alias `Mammarella2024PrimarySR`) in an FDHA logic-tree branch with `uncertaintyType="fdhaPrimarySRModel"`.

## Parameters (set inside the `uncertaintyModel` block)

| Name | Type | Units | Default | Allowed | Required? | Description |
| --- | --- | --- | --- | --- | --- | --- |
| `MSR` | int | – | – | `0`, `1`, `2` | Yes | Magnitude scaling relation for the rupture-width distribution: `0` = Leonard (2014) interplate, `1` = Leonard (2014) stable continental region, `2` = Thingbaijam et al. (2017). |
| `HDD_str` | string | – | – | `ITA_N`, `GB_N`, `AGG_N`, `ITA_R`, `TAI_R`, `JAP_R`, `AGG_R`, `CA_S`, `NZ_S`, `JAP_S`, `AGG_S` | Yes | Hypocentral depth distribution label (Table 2 of the paper); regional/style-specific mean and sigma of the hypocentral depth ratio. |
| `dip_sigma` | float | degrees | – | > 0 | Yes | Standard deviation of the fault dip distribution. |
| `t_d` | float | σ units | – | > 0 | Yes | Truncation of the dip distribution, in multiples of `dip_sigma`. |
| `Zs_sigma` | float | km | – | > 0 | Yes | Standard deviation of the seismogenic thickness distribution. |
| `t_z` | float | σ units | – | > 0 | Yes | Truncation of the seismogenic-thickness distribution, in multiples of `Zs_sigma`. |
| `style` | string | – | inferred from rake | `normal`, `reverse`, `strike-slip` | No | Style of faulting; if omitted, it is inferred from the rupture rake. |
| `seismothickness` | float | km | `15.0` | > 0 | No | Mean seismogenic thickness (used as `Zs_mu`). Set this explicitly if 15 km is not appropriate for your region. |

`dip_mu` (the mean dip) is supplied automatically by the calculator from the
fault-source geometry; it is not set in the configuration.

## Example logic-tree branch

```xml
<logicTreeBranch branchID="PSR_MAMMARELLA2024">
  <uncertaintyModel>
    [MammarellaEtAl2024PrimarySR]
    style = "normal"
    MSR = 1
    HDD_str = "AGG_N"
    seismothickness = 12.0
    dip_sigma = 10.0
    t_d = 2.0
    Zs_sigma = 2.0
    t_z = 2.0
  </uncertaintyModel>
  <uncertaintyWeight>0.25</uncertaintyWeight>
</logicTreeBranch>
```

## Notes and cautions

- Invalid `MSR` or `HDD_str` values raise `ValueError`.
- Probability is clipped to [0, 1].
- The default `seismothickness` of 15 km is a generic fallback; for site-specific studies, provide a regionally constrained value.
- Discretizations (number of log-width and dip integration points, hypocentral-depth-ratio grid, seismogenic-thickness step) are internal implementation details and not configurable via the configuration file.
