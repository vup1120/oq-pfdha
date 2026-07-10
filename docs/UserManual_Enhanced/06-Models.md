# Scientific Models Reference

This chapter lists the FDHA model classes registered in the runtime model library. In public jobs, these class names are selected inside FDHA logic-tree XML `<uncertaintyModel>` blocks, not in public `[models.*]` INI sections.

!!! note "Terminology: primary/secondary vs. principal/distributed"
    This toolkit names its four model categories `primary_*` and `secondary_*`,
    following the IAEA convention in which **primary = principal** (rupture on
    the main seismogenic fault) and **secondary = distributed** (off-fault
    rupture on splays, shears, and nearby structures). In the broader PFDHA
    literature — e.g. Valentini et al. (2025), *Reviews of Geophysics* — the
    preferred terms are **principal** and **distributed**, because
    "primary/secondary" is also used at a higher level to separate tectonic
    ("primary") from non-tectonic ("secondary") earthquake effects such as
    landsliding and liquefaction. In this toolkit, `secondary_surf_rup` and
    `secondary_surf_displ` always refer to **distributed tectonic** rupture and
    displacement — never to non-tectonic ground failure.

## Model Categories

For each calculation, select one model in each category; together they provide the conditional pieces of the PFDHA integral.

### Primary Surface Rupture
Models that provide the conditional probability that the principal (seismogenic) rupture reaches the ground surface on the primary fault—commonly denoted as the conditional probability of surface rupture (CPSR or P(Slip|M)). Approaches include empirical logistic regressions (e.g., magnitude‑dependent P(Slip|M)) and numerical, fault‑geometry–aware formulations that incorporate seismogenic thickness, dip, and rupture width distributions; the latter highlight the strong control of seismogenic depth and local geometry on CPSR.

### Primary Surface Displacement
Models that describe the probability distribution of displacement along the primary fault trace at a site, conditioned on surface rupture intersecting the site (and, implicitly, on the event magnitude and along‑strike position). Modern formulations distinguish displacement metrics (e.g., average, maximum, raw D, or normalized D/AD, D/MD) and what is being summed (single‑strand “single principal,” sum‑of‑principal across multiple strands, or aggregate definitions when applicable). Clear specification of vector components (lateral, vertical, net, dip‑slip) is essential because models predict different components.

### Secondary (Distributed) Surface Rupture
Models that give the probability of non‑zero, off‑fault rupture occurring in the distributed zone away from the principal trace, typically expressed as a function of magnitude and fault‑normal distance r to the principal rupture. In practice, this term is a major source of epistemic spread in PFDHA; inter‑team comparisons show order‑of‑magnitude differences in distributed hazard driven by different assumptions about the conditional probability of secondary rupturing.

### Secondary (Distributed) Surface Displacement
Models that describe the distribution of displacement amplitudes for secondary ruptures, conditioned on secondary rupture occurring at the site and parameterized by distance from the principal fault (and sometimes by along‑strike position or style of faulting). This term closes the distributed‑hazard integral (probability of any secondary rupture × probability of non‑zero displacement at r × exceedance distribution P[D>d0|⋅]) and is the piece that, together with the previous item, controls the shape of distributed‑displacement hazard curves.

#### Context notes for users
- The primary pair (items 1–2) applies on the principal trace; the secondary pair (items 3–4) applies off‑trace. This principal vs distributed distinction is standard in IAEA guidance and current exercises.
- For long return periods, hazard is often dominated by the surface‑rupture probability model choices (item 1), especially for moderate magnitudes—so multiple CPSR branches with justified weights are recommended.

### Displacement Metrics and Definitions for PFDHA

## 1) Displacement metrics (what models predict)

Recent reviews report that FD displacement models predict one or more of: average displacement (AD), maximum displacement (MD), normalized displacement (D/AD, D/MD), and the displacement amplitude D itself.

- AD, MD are event‑level quantities with magnitude scaling in several models (e.g., recent numerical and empirical models).
- D/AD, D/MD are normalized metrics used by multiple models (e.g., Youngs et al. 2003; Moss‑family; Mammarella 2024) typically with along‑strike position x/L as a predictor.
- Some models predict D directly (often after a transformation).

## 2) Displacement definitions (how displacement is defined)

Displacement is defined by two axes: (i) vector component — lateral, fault‑normal, vertical; net is vector sum; dip‑slip is vertical + fault‑normal; and (ii) participating ruptures — single principal, sum‑of‑principal, or aggregate (= principal + distributed).

!!! tip "Use models with their calibrated definition"
    Apply each model with the same displacement component and participation definition it was calibrated for; models predict different components and may not be interchangeable without conversion.

## 3) Principal (primary) vs distributed (secondary) surface ruptures

Principal ruptures occur on the primary fault that generated the earthquake; distributed ruptures occur off‑trace (splays, shears). They are modeled separately in data sets and in the PFDHA workflow (primary vs secondary branches).

## 4) Conditional probability of principal surface rupture (CPSR)

A numerical CPSR model for PFDHA computes the conditional probability of surface rupture on the principal fault as a function of magnitude, using down‑dip geometry and probabilistic inputs: hypocenter depth distribution (HDD), rupture width W(M) scaling, and hypocenter‑to‑rupture positioning (HDR).
---

## Primary Surface Rupture Models

These models answer the question: *Given an earthquake of a certain magnitude and style on a fault, what is the probability it will break the ground surface?*

| Model Class Name | Reference | Faulting style | Mw range¹ | Description |
| :--- | :--- | :--- | :--- | :--- |
| `Youngs2003PrimarySR` | Youngs et al. (2003) | Normal | 4.5–7.6 | A logistic regression model on magnitude. Its `style` parameter selects a dataset, not a generic faulting style: `all` uses the Wells & Coppersmith (1993) worldwide regression, `normal` the 32 Great Basin earthquakes (Pezzopane & Dawson 1996 data) from the Youngs et al. (2003) Appendix. |
| `MammarellaEtAl2024PrimarySR` | Mammarella et al. (2024) | All (numerical) | 5.0–8.0 | A numerical integration model considering rupture width, dip, and other physical parameters. |
| `Mammarella2024PrimarySR` | Mammarella et al. (2024) | All (numerical) | 5.0–8.0 | Alias/implementation class for the Mammarella et al. (2024) primary surface rupture model. |
| `MossRoss2011PrimarySR` | Moss & Ross (2011) | Reverse | 5.5–8.0 | A logistic regression model specifically for reverse faults. |
| `Moss2013PrimarySR` | Moss et al. (2013) | Reverse, strike-slip | 4.2–8.7 | A model for reverse and strike-slip faults that accounts for site stiffness (Vs30). |
| `Pizza2023PrimarySR` | Pizza et al. (2023) | All (N / R / SS subsets) | 5.5–7.9 | A logistic regression model from an updated global database of surface ruptures. |
| `Takao2013PrimarySR` | Takao et al. (2013) | Reverse, strike-slip | 5.5–7.4 | A model developed for reverse and strike-slip faults in Japan. |
| `WC1993PrimarySR` | Wells & Coppersmith (1993) | All | 5.0–8.2 | An early logistic regression model for surface rupture probability. |
| `Yang2021PrimarySR` | Yang et al. (2021) | Reverse | 4.7–6.6 | A logistic regression model for reverse faults developed using Australian earthquake data. |
| `FixedPrimarySR` | Fixed value | — | — | Constant primary surface rupture probability model (not data-derived). |

¹ Applicable moment-magnitude range of the underlying empirical/numerical model,
from Valentini et al. (2025), *Reviews of Geophysics*, Table 4. Using a model
outside this range extrapolates beyond its calibration data.

---

## Primary Surface Displacement Models

These models answer the question: *Given that a surface rupture has occurred, what is the probability that the displacement at a point will exceed a certain amount?*

| Model Class Name | Reference | Faulting style | Slip component | Participation | Mw range¹ | Description |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `Youngs2003PrimaryFD` | Youngs et al. (2003) | Normal | Vertical | Single principal | via AD/MD scaling | Models the displacement profile based on Average or Maximum Displacement normalization. |
| `Chiou2025PrimaryFD` | Chiou et al. (2025) | Strike-slip | Net | Sum-of-principal | 6.0–8.3 | A model for sum-of-principal displacement on strike-slip faults (NGA-Displacement). |
| `Kuehn2024PrimaryFD` | Kuehn et al. (2024) | All | Net | Aggregate | 5.0–8.0 (R); 6.0–8.0 (N, SS) | A comprehensive model with options to include epistemic uncertainty via posterior sampling. |
| `Lavrentiadis2023PrimaryFD`| Lavrentiadis & Abrahamson (2023) | All | Net | Aggregate / sum-of-principal | 5.0–8.5 | A model that can account for zero-slip probability and rupture gaps. |
| `MossRoss2011PrimaryFD` | Moss & Ross (2011) | Reverse | Vertical | Single principal | via AD/MD scaling | Restored legacy normalized displacement model with AD/MD scaling. |
| `Moss2022PrimaryFD` | Moss et al. (2022) | Reverse | Vertical | Single principal | 4.7–8.0 | GIRS-2022-05 report formulation with `gamma_mode`, incomplete MD subset, and `sigma_type` selector. |
| `Moss2024PrimaryFD` | Moss et al. (2024) | Reverse | Vertical | Single principal | 4.7–8.0 | Peer-reviewed *Earthquake Spectra* implementation; `source="EQS"` uses journal Table 2 alpha/beta files, while `source="GIRS"` uses GIRS gamma regressions. AD/MD scaling follows Table 3. |
| `Petersen2011PrimaryFD` | Petersen et al. (2011) | Strike-slip | Lateral | Single principal | 6.0–8.0 | Base Petersen primary displacement model. |
| `Petersen2011PrimaryFD_bilinear` | Petersen et al. (2011) | Strike-slip | Lateral | Single principal | 6.0–8.0 | Bilinear strike-slip displacement profile. |
| `Petersen2011PrimaryFD_elliptical` | Petersen et al. (2011) | Strike-slip | Lateral | Single principal | 6.0–8.0 | Elliptical strike-slip displacement profile. |
| `Petersen2011PrimaryFD_quadratic` | Petersen et al. (2011) | Strike-slip | Lateral | Single principal | 6.0–8.0 | Quadratic strike-slip displacement profile. |
| `Takao2013PrimaryFD` | Takao et al. (2013) | Reverse, strike-slip | Net | Single principal | via AD/MD scaling | A normalized displacement model with coefficients dependent on surface rupture length. |

¹ Applicable Mw range, from Valentini et al. (2025), *Reviews of Geophysics*,
Table 4. "via AD/MD scaling" marks normalized models in which magnitude enters
through a separate average/maximum-displacement scaling relation rather than as
a direct model input. *Slip component* and *Participation* (single-principal,
sum-of-principal, or aggregate) are the calibrated definitions — apply each
model only with its own definition (see the tip in "Displacement Metrics").

---

## Secondary Surface Rupture Models

These models answer the question: *Given an earthquake, what is the probability of surface rupture at a distance away from the primary fault?*

| Model Class Name | Reference | Faulting style | Mw range¹ | r range¹ | Description |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `Youngs2003SecondarySR` | Youngs et al. (2003) | Normal | 5.5–7.4 | 0–15 km | A logistic regression model based on magnitude, distance, and hanging wall location. |
| `Petersen2011SecondarySR` | Petersen et al. (2011) | Strike-slip | 6.5–7.5 | 0–2.5 km | A model for strike-slip faults where probability depends on distance and grid cell size. |
| `Takao2014SecondarySR` | Takao et al. (2014) | Reverse, strike-slip | 5.8–7.4 | 0–25 km | A model for reverse and strike-slip faults based on distance and pixel size. |
| `Visini2025SecondarySR` | Visini et al. (2025) | Normal, reverse | 5.5–7.9 (N); 4.9–7.9 (R) | 0–10 km (HW); 0–8 km (FW) | A logistic regression for normal/reverse faults depending on magnitude, distance, and pixel size. |
| `FerrarioLivio2021SecondarySR` | Ferrario & Livio (2021) | Normal | 6.0–7.5 | 0–15.5 km (HW); 0–12.5 km (FW) | A model for distributed surface rupture probability. |
| `Rodriguez2023SecondarySR` | Rodriguez Padilla & Oskin (2023) | Strike-slip | — | 0–3 km | Probability per unit area; recommended for near-field, immature strike-slip faults. |
| `Takao2013SecondarySR` | Takao et al. (2013) | Reverse, strike-slip | 5.8–7.4 | 0–25 km | A model for distributed surface rupture probability. |
| `Petersen2011SecondarySR_default` | Petersen et al. (2011) | Strike-slip | 6.5–7.5 | 0–2.5 km | Default Petersen secondary surface rupture variant exposed by the library. |
| `Moss2022SecondarySR` | Moss et al. (2022) | Reverse | Report-specific | Report-specific | Distributed surface-rupture probability model from GIRS-2022-05 Section 5.2.1; simple mode uses Eq. 5.5 / Table 5.3 and biexponential mode uses Eqs. 5.6-5.7 / Tables 5.4-5.5. |
| `FixedSecondarySR` | Fixed value | — | — | — | Constant secondary surface rupture probability model (not data-derived). |

¹ Applicable Mw and fault-normal-distance (r) ranges, from Valentini et al.
(2025), *Reviews of Geophysics*, Table 4. HW = hanging wall, FW = footwall.
Using a model outside these ranges extrapolates beyond its calibration data.
Moss et al. (2022) secondary models are documented directly in GIRS-2022-05
Section 5, rather than in the Valentini et al. summary table.

---

## Secondary Surface Displacement Models

These models answer the question: *Given that secondary rupture has occurred, what is the probability the displacement will exceed a certain amount?*

| Model Class Name | Reference | Faulting style | Slip component | Mw range¹ | r range¹ | Description |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| `Youngs2003SecondaryFD` | Youngs et al. (2003) | Normal | Vertical | 5.5–7.4 | 0–15 km | Models distributed displacement as a fraction of the Maximum Displacement on the principal fault. |
| `Petersen2011SecondaryFD` | Petersen et al. (2011) | Strike-slip | Lateral | 6.5–7.5 | 0–2.5 km | Provides exceedance probability for distributed displacement on strike-slip faults. |
| `Moss2022SecondaryFD` | Moss et al. (2022) | Reverse | Vertical distributed displacement normalized by MD/AD | Report-specific | Report-specific | Distributed displacement model from GIRS-2022-05 Section 5.2.3; envelope mode uses Eq. 5.8 / Tables 5.7-5.8 and gamma mode combines the report's global gamma distribution with the distance envelope. |
| `Visini2025SecondaryFD` | Visini et al. (2025) | Normal, reverse | Vertical | 5.5–7.9 (N); 4.9–7.9 (R) | 0–10 km (HW); 0–8 km (FW) | A regression model for normal/reverse faults predicting median displacement from magnitude, distance, and mean throw. |

¹ Applicable Mw and fault-normal-distance (r) ranges, from Valentini et al.
(2025), *Reviews of Geophysics*, Table 4. HW = hanging wall, FW = footwall.
Moss et al. (2022) secondary models are documented directly in GIRS-2022-05
Section 5, rather than in the Valentini et al. summary table.

!!! note "Validation metadata"
    The logic-tree validator has citation-traced displacement-definition metadata only for a subset of primary displacement classes (`Youngs2003PrimaryFD`, Petersen bilinear/elliptical/quadratic variants, and `Chiou2025PrimaryFD`). Classes absent from that metadata are still registered if they are imported by the runtime library; they are simply skipped by the mixed-displacement-definition advisory check.

!!! note "Detailed Model Documentation"
    The specific parameters for each model are crucial for correct implementation. Detailed configuration guides are available for the following models:

    **Primary Surface Rupture Models:**
    - [Mammarella et al. (2024)](models/primary/MammarellaEtAl2024.md) — Primary surface rupture probability (numerical)
    - [Yang et al. (2021)](models/primary/Yang2021.md) — Reverse faults, logistic regression (Australian data)

    **Primary Surface Displacement Models:**
    - [Chiou et al. (2025)](models/primary/Chiou2025.md) — Strike-slip faults, sum-of-principal displacement
    - [Kuehn et al. (2024)](models/primary/Kuehn2024.md) — All fault styles, aggregate displacement with epistemic uncertainty
    - [Lavrentiadis & Abrahamson (2023)](models/primary/Lavrentiadis2023.md) — All fault styles, aggregate/principal displacement
    - [Moss et al. (2024)](models/primary/Moss2024.md) — Reverse faults, normalized displacement
    - [Petersen et al. (2011)](models/primary/Petersen2011.md) — Strike-slip faults, multiple functional forms
    - [Takao et al. (2013)](models/primary/Takao2013.md) — Reverse and strike-slip faults, normalized displacement
    - [Youngs et al. (2003)](models/primary/Youngs2003.md) — Normal faults, normalized displacement

    **Secondary Surface Displacement Models:**
    - [Petersen et al. (2011)](models/secondary/Petersen2011.md) — Strike-slip faults, distributed displacement
    - [Visini et al. (2025)](models/secondary/VisiniEtAl2025.md) — Normal and reverse faults, distributed displacement (Note: Model class name is `Visini2025SecondarySR` and `Visini2025SecondaryFD`)
    - [Youngs et al. (2003)](models/secondary/Youngs2003.md) — Normal faults, distributed displacement

    Each guide includes complete parameter tables, INI configuration examples, implementation details, and usage notes.
