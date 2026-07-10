# Aleatory and Epistemic Uncertainty in the Model Library

This chapter documents how each FDHA model registered in this toolkit treats
**aleatory variability** (the irreducible scatter in the conditional
distribution, controlled by σ, τ, φ, or a discrete probability mass) and
**epistemic uncertainty** (the choice between alternative coefficients,
formulations, or datasets, normally captured as logic-tree branches). The
descriptions below are extracted directly from the source code; the file/line
references point to the lines that compute the standard deviation, evaluate the
distribution, and combine the discrete and continuous components.

All four model categories are covered:

- §1 Primary surface rupture (PSR) — `openquake/fdha/primary_surf_rup/`
- §2 Primary surface displacement (PSD) — `openquake/fdha/primary_surf_displ/`
- §3 Secondary (distributed) surface rupture (SSR) — `openquake/fdha/secondary_surf_rup/`
- §4 Secondary (distributed) surface displacement (SSD) — `openquake/fdha/secondary_surf_displ/`
- §5 Cross-cutting patterns and recommendations for logic-tree construction

For each model, the entry lists: the **distribution family** used to evaluate
exceedance probability; the **σ (and τ/φ) parameters** with their functional
dependence on magnitude and along-strike position; the **truncation**
convention; any **discrete zero-displacement** component; and the **epistemic
branches** exposed as user options.

---

## 1. Primary Surface Rupture (PSR) Models

PSR models predict the conditional probability `P(SR|M)` that a rupture of
given magnitude reaches the surface. With one exception they all use a
**logistic regression** of the form

```
P(SR|M) = 1 / (1 + exp(a − b·M))     (or equivalent exp(fx)/(1+exp(fx)))
```

so the output is intrinsically bounded to [0, 1] by the link function and no
explicit aleatory σ is carried on `P` itself. Epistemic uncertainty is exposed
as alternative coefficient sets (style, regional dataset, site stiffness).

### `FixedPrimarySR` — [fixed.py](../../openquake/fdha/primary_surf_rup/fixed.py)
- Deterministic constant probability supplied by the user (default 1.0).
- No aleatory or epistemic structure.

### `WC1993PrimarySR` — [wells_coppersmith1993.py:51](../../openquake/fdha/primary_surf_rup/wells_coppersmith1993.py#L51)
- Logistic, magnitude only: `a = −12.51`, `b = 2.053`.
- Single deterministic model; no branches.

### `Youngs2003PrimarySR` — [youngs2003.py:44-46](../../openquake/fdha/primary_surf_rup/youngs2003.py#L44-L46)
- Logistic, magnitude only.
- **Epistemic branches:** `style="all"` (`a=−12.51`, `b=2.053`) vs `style="normal"` (`a=−16.02`, `b=2.685`).

### `MossRoss2011PrimarySR` — [moss_ross2011.py:50](../../openquake/fdha/primary_surf_rup/moss_ross2011.py#L50)
- Logistic, magnitude only: `a = 7.3`, `b = 1.03`. Reverse-faulting calibration.

### `Moss2013PrimarySR` — [moss2013.py:86-107](../../openquake/fdha/primary_surf_rup/moss2013.py#L86-L107)
- Logistic, magnitude only.
- **Epistemic branches:** four coefficient sets covering style (reverse / strike-slip) × site stiffness (Vs30 > 600 m/s or ≤ 600 m/s).

### `Takao2013PrimarySR` — [takao2013.py:54](../../openquake/fdha/primary_surf_rup/takao2013.py#L54)
- Logistic, magnitude only: `a = −32.03`, `b = 4.9`.

### `Yang2021PrimarySR` — [yang2021.py:59](../../openquake/fdha/primary_surf_rup/yang2021.py#L59)
- Logistic, magnitude only: `a = −24.59`, `b = 4.0`.

### `Pizza2023PrimarySR` — [pizza2023.py:74-83](../../openquake/fdha/primary_surf_rup/pizza2023.py#L74-L83)
- Logistic, magnitude only.
- **Epistemic branches:** four style options — `all`, `normal`, `reverse`, `strike-slip`.

### `Mammarella2024PrimarySR` — [mammarella2024.py](../../openquake/fdha/primary_surf_rup/mammarella2024.py)
- **Departs from the logistic pattern.** Computes `P(SR|M)` by a discrete grid
  integration over geometric and physical inputs:
  - log₁₀ W (rupture width) drawn from a **truncated lognormal**, half-width `T_W·σ = 1.0σ` ([mammarella2024.py:118](../../openquake/fdha/primary_surf_rup/mammarella2024.py#L118));
  - dip from a **truncated normal**, half-width `t_d·σ`;
  - seismogenic thickness `Zs` from a **truncated normal**, half-width `t_z·σ`;
  - hypocentre-depth ratio (HDR) from **uniform** templates; hypocentre-depth distribution (HDD) ratio from discretized **normal** templates (`TAB2`).
- **No aleatory σ on P.** All spread is geometric/epistemic. Magnitude–width
  scaling has three MSR options (codes 0/1/2, [lines 85–95](../../openquake/fdha/primary_surf_rup/mammarella2024.py#L85-L95)) and 11 HDD templates ([`TAB2`, lines 99–104](../../openquake/fdha/primary_surf_rup/mammarella2024.py#L99-L104): 3 normal, 4 reverse, 4 strike-slip) — together a large epistemic branch set.

---

## 2. Primary Surface Displacement (PSD) Models

PSD models predict the conditional exceedance probability `P(D > d | M, x/L)`.
Distributional choices diverge widely, but two structural families dominate:

- **Magnitude-scaling + normalized-displacement convolution** (Youngs 2003,
  Moss & Ross 2011, Takao 2013, Moss 2022/2024): magnitude enters through a
  log₁₀-normal scaling of AD or MD, then the ratio `D/AD` follows a Gamma
  distribution and `D/MD` a Beta distribution, with shape/scale parameters
  parameterised in x/L.
- **Direct lognormal / transformed-normal on `D`** (Petersen 2011, Lavrentiadis
  2023, Kuehn 2024, Chiou 2025): a (possibly power-transformed) Normal CDF is
  evaluated on `ln(D)` or `D^λ`, with magnitude and x/L entering the mean and
  σ.

### `Youngs2003PrimaryFD` — [youngs2003.py](../../openquake/fdha/primary_surf_displ/youngs2003.py)
- **Distribution:** Gamma on `D/AD` via `gamma.sf` ([line 146](../../openquake/fdha/primary_surf_displ/youngs2003.py#L146)); Beta on `D/MD` via `beta.cdf` renormalised at D=1 ([lines 150–151](../../openquake/fdha/primary_surf_displ/youngs2003.py#L150-L151)); magnitude scaling log₁₀-normal. (The standalone `get_prob_D_AD`/`get_prob_D_MD` helpers use `gamma.sf`/`beta.sf` at lines 178/197.)
- **σ:** Fixed in log₁₀ space from Wells & Coppersmith (1994): `σ_AD = 0.36`, `σ_MD = 0.42` ([lines 34–35](../../openquake/fdha/primary_surf_displ/youngs2003.py#L34-L35)). Gamma α, β are x/L-dependent exponentials ([lines 114–119](../../openquake/fdha/primary_surf_displ/youngs2003.py#L114-L119)).
- **τ/φ:** Not separated.
- **Truncation:** ±`n_sigma`σ ε-space integration with step 0.1; `n_sigma` is a constructor parameter (default 6) settable from the logic tree via `[Youngs2003PrimaryFD] n_sigma = <value>`.
- **Zero-displacement:** Implicit through integration bounds.
- **Epistemic branches:** `style="all"` vs `style="normal"`.

### `MossRoss2011PrimaryFD` — [moss_ross2011.py](../../openquake/fdha/primary_surf_displ/moss_ross2011.py)
- **Distribution:** Gamma on `D/AD` ([lines 104–109](../../openquake/fdha/primary_surf_displ/moss_ross2011.py#L104-L109)); Beta on `D/MD` ([lines 134–136](../../openquake/fdha/primary_surf_displ/moss_ross2011.py#L134-L136)).
- **σ:** `σ_AD = 0.17`, `σ_MD = 0.31` in log₁₀ ([lines 125, 143](../../openquake/fdha/primary_surf_displ/moss_ross2011.py#L125)). Gamma `a,b` polynomial in x/L; Beta `α,β` linear in x/L.
- **Truncation:** ±`n_sigma`σ — the conditional AD/MD log₁₀-normal is normalised over a log-spaced grid spanning `mean ± n_sigma·σ` (per displacement type), matching Takao 2013. `n_sigma` is a constructor parameter (default 3) settable from the logic tree via `[MossRoss2011PrimaryFD] n_sigma = <value>`. (Previously normalised over an arbitrary fixed 0.001–10 m grid.)
- **τ/φ:** Not separated. No style branches, though `get_prob_D_AD` exposes a gamma-vs-Weibull `variant` option.

### `Petersen2011PrimaryFD` (and `_bilinear`, `_quadratic`, `_elliptical`) — [petersen2011.py](../../openquake/fdha/primary_surf_displ/petersen2011.py)
- **Distribution:** Lognormal in *natural log* `ln(D in cm)`, evaluated via `norm.cdf` ([line 80](../../openquake/fdha/primary_surf_displ/petersen2011.py#L80)).
- **σ (fixed in ln-cm):**
  - Bilinear: `σ₁ = 1.2906` (low x/L), `σ₂ = 0.9624` (high x/L) ([line 97](../../openquake/fdha/primary_surf_displ/petersen2011.py#L97));
  - Quadratic: `σ = 1.1346` ([line 149](../../openquake/fdha/primary_surf_displ/petersen2011.py#L149));
  - Elliptical: `σ = 1.1348` ([line 128](../../openquake/fdha/primary_surf_displ/petersen2011.py#L128)).
- **Mean:** Magnitude- and x/L-dependent (piecewise / polynomial / elliptical).
- **τ/φ:** Not separated. **Truncation:** none — `norm.cdf` applied directly.
- **Epistemic branches:** the three regression variants are exposed as
  separate registered classes.

### `Takao2013PrimaryFD` — [takao2013.py](../../openquake/fdha/primary_surf_displ/takao2013.py)
- **Distribution:** Gamma on `D/AD` (`1 − gamma.cdf`, [line 139](../../openquake/fdha/primary_surf_displ/takao2013.py#L139)); Beta on `D/MD` (`1 − beta.cdf`, [line 191](../../openquake/fdha/primary_surf_displ/takao2013.py#L191)); log₁₀-normal magnitude scaling.
- **σ:** Wells & Coppersmith — `σ_AD = 0.36`, `σ_MD = 0.42` in log₁₀ ([lines 81, 86](../../openquake/fdha/primary_surf_displ/takao2013.py#L81)). Gamma `α, β` (and Beta `α, β` for MD) switch on SRL < 10 km (fixed) vs ≥ 10 km (x/L-dependent) ([lines 132–137, 184–189](../../openquake/fdha/primary_surf_displ/takao2013.py#L132-L137)).
- **Truncation:** ±`n_sigma`σ; logspace integration over 1000 points. `n_sigma` is a constructor parameter (default 3) settable from the logic tree via `[Takao2013PrimaryFD] n_sigma = <value>`.
- **τ/φ:** Not separated.

### `Moss2022PrimaryFD` — [moss2022.py](../../openquake/fdha/primary_surf_displ/moss2022.py) and `Moss2024PrimaryFD` — [moss2024.py](../../openquake/fdha/primary_surf_displ/moss2024.py)
- **Distribution:** Gamma on `D/XD` via `gamma.cdf` ([moss2024.py:103](../../openquake/fdha/primary_surf_displ/moss2024.py#L103)); magnitude scaling log₁₀-normal via numerical integration.
- **σ (alternatives):**
  - "recommended" σ — 0.20 (MD-complete), 0.20 (AD-complete), 0.25 (AD-"all");
  - "regression" σ — 0.148 (MD-complete) / 0.133 (AD-complete).
  - (These scaling tables live in [moss2022.py:19–28](../../openquake/fdha/primary_surf_displ/moss2022.py#L19-L28).)
- **Gamma α, β source (Moss 2024):**
  - `source="EQS"` (default) — interpolated from *Earthquake Spectra* Table 2 ([lines 133–135](../../openquake/fdha/primary_surf_displ/moss2024.py#L133-L135));
  - `source="GIRS"` — x/L regression from GIRS-2022-05 Figures 4.3–4.4 ([lines 126–130](../../openquake/fdha/primary_surf_displ/moss2024.py#L126-L130)).
- **Truncation:** ±6σ ε-space, step 0.1.
- **τ/φ:** Not separated.
- **Epistemic branches:** completeness (complete / all), σ source (recommended / regression), gamma source (EQS / GIRS).

### `Lavrentiadis2023PrimaryFD` — [lavrentiadis2023.py](../../openquake/fdha/primary_surf_displ/lavrentiadis2023.py)
- **Distribution:** Normal on the **power-transformed displacement** `D^0.3`, via `norm.sf` ([lines 94–96](../../openquake/fdha/primary_surf_displ/lavrentiadis2023.py#L94-L96)). The transformation is applied as `D = μ_prime^(1/0.3)` ([line 345](../../openquake/fdha/primary_surf_displ/lavrentiadis2023.py#L345)), so probability statements live in the "prime" (Box–Cox-like) space.
- **σ — full τ/φ decomposition, both magnitude-dependent** ([lines 329–340](../../openquake/fdha/primary_surf_displ/lavrentiadis2023.py#L329-L340)):
  - between-event `τ_agg = clip(0.115 + 0.060·(M−6), [0.115, 0.205])`;
  - within-event `φ_agg = clip(0.120 + 0.150·(M−6), [0.120, 0.270])`;
  - principal φ adds a component-correlation term with `ρ = −0.15`;
  - additional segmentation variance `φ_add = c₁₈ + c₁₉·M + c₂₀·(M−6.7)²`;
  - total `σ_total = √(τ² + φ² + φ_add²)`.
- **Discrete zero-displacement components — unique in the library:**
  - `P_gap` (segment-gap probability, [lines 293–310](../../openquake/fdha/primary_surf_displ/lavrentiadis2023.py#L293-L310));
  - `P_zero_slip` (logistic, [line 317](../../openquake/fdha/primary_surf_displ/lavrentiadis2023.py#L317));
  - combined as `ccdf_prnc · (1 − P_zero_slip) · (1 − P_gap)` ([line 99](../../openquake/fdha/primary_surf_displ/lavrentiadis2023.py#L99)).
- **Epistemic branches:** style-dependent coefficient sets; three displacement
  metrics (`disp_agg_prime`, `disp_prnc_prime`, `disp_agg_seg`).

### `Kuehn2024PrimaryFD` — [kuehn2024/kuehn2024.py](../../openquake/fdha/primary_surf_displ/kuehn2024/kuehn2024.py)
- **Distribution:** **Box–Cox-transformed Normal** — `norm.cdf` on `(D^λ − 1)/λ` (or `ln D` if λ = 0), [lines 95–104, 161–169](../../openquake/fdha/primary_surf_displ/kuehn2024/kuehn2024.py#L95-L104).
- **σ — explicit between/within decomposition** combined as `σ_total = √(σ_mode² + σ_within²)` ([lines 197, 211, 224](../../openquake/fdha/primary_surf_displ/kuehn2024/kuehn2024.py#L197)):
  - between-event `σ_mode` is magnitude-dependent — bilinear hinge at M = 7.0 for strike-slip ([lines 256–259](../../openquake/fdha/primary_surf_displ/kuehn2024/kuehn2024.py#L256-L259)), sigmoid for normal ([lines 261–263](../../openquake/fdha/primary_surf_displ/kuehn2024/kuehn2024.py#L261-L263)), fixed `s_m,r` for reverse ([line 222](../../openquake/fdha/primary_surf_displ/kuehn2024/kuehn2024.py#L222));
  - within-event `σ_within` is x/L-dependent quadratic for SS/RV ([lines 265–283](../../openquake/fdha/primary_surf_displ/kuehn2024/kuehn2024.py#L265-L283)); constant for normal ([lines 208–209](../../openquake/fdha/primary_surf_displ/kuehn2024/kuehn2024.py#L208-L209)).
- **Epistemic — uniquely propagated as a full ensemble:** with
  `epistemic_uncertainty=True` the model returns one probability curve per
  posterior coefficient sample ([lines 76–120](../../openquake/fdha/primary_surf_displ/kuehn2024/kuehn2024.py#L76-L120)); with `False` it returns a single mean-coefficient curve ([lines 122–180](../../openquake/fdha/primary_surf_displ/kuehn2024/kuehn2024.py#L122-L180)).
- **Truncation:** None explicit.

### `Chiou2025PrimaryFD` — [chiou2025.py](../../openquake/fdha/primary_surf_displ/chiou2025.py)
- **Distribution:** **Negative exponentially-modified Gaussian (nEMG)** via `stats.exponnorm.cdf` ([line 176](../../openquake/fdha/primary_surf_displ/chiou2025.py#L176)) — an asymmetric Gaussian + exponential combination.
- **σ — explicit magnitude/position decomposition** ([lines 158–164](../../openquake/fdha/primary_surf_displ/chiou2025.py#L158-L164)):
  - magnitude component `σ_mag = max(0.4, cv₁·exp(cv₂·max(0, M − 6.1)))` (i.e. a **lower floor** of 0.4);
  - position component `σ_xl = cv₃·exp(cv₄·max(0, x_fold − ccap))`;
  - combined `σ' = √(σ_mag² + σ_xl²)`, with mixing shape `K = cv₅ / σ'`.
- **τ/φ:** Not labelled as between/within event, but the σ_mag / σ_xl split is
  structurally similar.
- **Epistemic branches:** four model variants — `model7`, `model8.1`, `model8.2`, `model8.3`.

---

## 3. Secondary (Distributed) Surface Rupture (SSR) Models

SSR models predict `P(distributed rupture | M, r)`. As with PSR, they carry no
explicit aleatory σ on the probability itself; the output is bounded to [0, 1]
by either a logistic link or by explicit clipping. The functional form (logistic
in `M, r`; power law in `r`; bi-exponential; etc.) is itself the main epistemic
choice.

### `FixedSecondarySR` — [fixed.py](../../openquake/fdha/secondary_surf_rup/fixed.py)
- Deterministic constant probability.

### `Youngs2003SecondarySR` — [youngs2003.py](../../openquake/fdha/secondary_surf_rup/youngs2003.py)
- Logistic in distance with hanging-wall indicator h: `fx = const + (c₁ + c₂·h)·log(r + offset)` (h at [line 59](../../openquake/fdha/secondary_surf_rup/youngs2003.py#L59); v1/v2 formulas at [lines 62, 66](../../openquake/fdha/secondary_surf_rup/youngs2003.py#L62-L66)).
- **Epistemic branches:** three versions — v1, v2, and v3 = average of v1/v2.

### `Petersen2011SecondarySR` — [petersen2011.py](../../openquake/fdha/secondary_surf_rup/petersen2011.py)
- Far-field power law `ln(P) = a·ln(r) + b` ([line 112](../../openquake/fdha/secondary_surf_rup/petersen2011.py#L112)); near-field linear interpolation between `p₀, p₁, p₂` ([lines 121–133](../../openquake/fdha/secondary_surf_rup/petersen2011.py#L121-L133)).
- **Epistemic branches:** cell-size-dependent coefficient tables (25, 50, 100, 150, 200 m, [lines 35–41](../../openquake/fdha/secondary_surf_rup/petersen2011.py#L35-L41)); default-vs-near-field method.

### `Takao2013SecondarySR` — [takao2013.py:87](../../openquake/fdha/secondary_surf_rup/takao2013.py#L87) and `Takao2014SecondarySR` — [takao2014.py:54](../../openquake/fdha/secondary_surf_rup/takao2014.py#L54)
- Takao 2013: logistic with `fx = C₁ + (C₂ + C₃·M)·log(r + C₄)`.
- Takao 2014: logistic — `fx = C₁ + C₂·ln(r + C₃)` then `P = exp(fx)/(1 + exp(fx))` ([takao2014.py:54, 57](../../openquake/fdha/secondary_surf_rup/takao2014.py#L54)); pixel-size-dependent coefficients (50–500 m).

### `Moss2022SecondarySR` — [moss2022.py](../../openquake/fdha/secondary_surf_rup/moss2022.py)
- **Two methods:** simple exponential `P = min(exp(−a·r + b), 1)` ([lines 93–100](../../openquake/fdha/secondary_surf_rup/moss2022.py#L93-L100)); bi-exponential CDF `F(x) = a·exp(b·x) + c·exp(d·x)` with `P = clip(1 − F, 0, 1)` ([lines 118–120](../../openquake/fdha/secondary_surf_rup/moss2022.py#L118-L120)).
- Magnitude-binned (M ≥ 7, 6–7, 5–6) and HW/FW-asymmetric coefficient tables.

### `FerrarioLivio2021SecondarySR` — [ferrario2021.py:114](../../openquake/fdha/secondary_surf_rup/ferrario2021.py#L114)
- Logistic in `ln(r)`. **Epistemic branches:** "regular" vs "conservative" × HW/FW.

### `Rodriguez2023SecondarySR` — [rodriguez2023.py:76](../../openquake/fdha/secondary_surf_rup/rodriguez2023.py#L76)
- Power law `P = a·((r + b)/b)^c`, fixed coefficients calibrated at 1 m pixel.

### `VisiniEtAl2025SecondarySR` — [visini2025.py](../../openquake/fdha/secondary_surf_rup/visini2025.py)
- **Compound model:** logistic `P_slice = 1/(1 + exp(a + b·M + c·r + d·fw))` ([lines 207, 212–213](../../openquake/fdha/secondary_surf_rup/visini2025.py#L207-L213)) **combined with** a Monte-Carlo integration over truncated-lognormal distributed-rupture segment lengths (style/HW/FW-dependent μ, σ at [lines 114–123](../../openquake/fdha/secondary_surf_rup/visini2025.py#L114-L123); MC draw at [lines 381–384](../../openquake/fdha/secondary_surf_rup/visini2025.py#L381-L384)).
- **Truncation:** segment-length sampling default `"truncated"` (clipped to fixed per-mechanism min/max bounds, `_drlengths_min_max`, [lines 120–123](../../openquake/fdha/secondary_surf_rup/visini2025.py#L120-L123)), with a `"legacy"` alternative (raw lognormal clipped to `[10, fault_length]`; validated at [line 334](../../openquake/fdha/secondary_surf_rup/visini2025.py#L334)).
- **Epistemic branches:** combinations A, B, C × style × pixel size × segment-sampling mode.

---

## 4. Secondary (Distributed) Surface Displacement (SSD) Models

### `Youngs2003SecondaryFD` — [youngs2003.py](../../openquake/fdha/secondary_surf_displ/youngs2003.py)
- **Distribution:** Gamma on `D/MD` via `gamma.sf` ([line 149](../../openquake/fdha/secondary_surf_displ/youngs2003.py#L149)); MD log₁₀-normal with `σ = 0.38` (normal-fault calibration, [line 41](../../openquake/fdha/secondary_surf_displ/youngs2003.py#L41)).
- **Scale:** Gamma `b = x / scale_factor` with `x` and decay constants HW/FW- and `r`-dependent ([lines 139, 143](../../openquake/fdha/secondary_surf_displ/youngs2003.py#L139)); shape `a = 2.5` ([line 46](../../openquake/fdha/secondary_surf_displ/youngs2003.py#L46)).
- **Truncation:** ±3σ in log₁₀ ([line 43](../../openquake/fdha/secondary_surf_displ/youngs2003.py#L43)).
- **Epistemic branches:** percentile option (`"85"` / `"95"`, [lines 45–54](../../openquake/fdha/secondary_surf_displ/youngs2003.py#L45-L54)); HW/FW split.

### `Petersen2011SecondaryFD` — [petersen2011.py](../../openquake/fdha/secondary_surf_displ/petersen2011.py)
- **Distribution:** Lognormal in `ln(D in cm)` via `norm.cdf` ([line 133](../../openquake/fdha/secondary_surf_displ/petersen2011.py#L133)).
- **σ:** `σ_dist = 1.1193` in `ln(cm)` ([line 123](../../openquake/fdha/secondary_surf_displ/petersen2011.py#L123)).
- **Mean:** `μ = 1.4016·M − 0.1671·ln(r) − 6.7991` ([line 120](../../openquake/fdha/secondary_surf_displ/petersen2011.py#L120)).
- **τ/φ:** Not separated. **Truncation:** none explicit. Single deterministic
  formulation.

### `Moss2022SecondaryFD` — [moss2022.py](../../openquake/fdha/secondary_surf_displ/moss2022.py)
- **Distribution alternatives:** Gamma method ([line 192](../../openquake/fdha/secondary_surf_displ/moss2022.py#L192)) **or** envelope (lognormal on MD) method ([line 216](../../openquake/fdha/secondary_surf_displ/moss2022.py#L216)).
- **σ:** Gamma uses global α/β (Table 3) with the scale rescaled per site by a distance envelope ([lines 179–181](../../openquake/fdha/secondary_surf_displ/moss2022.py#L179-L181)); envelope method uses log₁₀ σ from Table 4.4 ([lines 160–162](../../openquake/fdha/secondary_surf_displ/moss2022.py#L160-L162)).
- **Truncation:** Gamma ±6σ in ε-space ([lines 184–186](../../openquake/fdha/secondary_surf_displ/moss2022.py#L184-L186)).
- **Epistemic branches:** AD vs MD; completeness (complete / incomplete / all); σ source (recommended / regression); method (gamma / envelope).

### `VisiniEtAl2025SecondaryFD` — [visini2025.py](../../openquake/fdha/secondary_surf_displ/visini2025.py)
- **Distribution:** Lognormal in `ln(Y)` via `norm.cdf` ([lines 188, 191](../../openquake/fdha/secondary_surf_displ/visini2025.py#L188)).
- **σ:** `σ = 1.0271` in `ln(Y)` ([line 78](../../openquake/fdha/secondary_surf_displ/visini2025.py#L78)).
- **Mean:** `ln(Y_med) = a + b·ln(s) + c·ln(TPFm) + d·M + e·I_style + f·I_fw + g_offset` ([lines 232–240](../../openquake/fdha/secondary_surf_displ/visini2025.py#L232-L240)).
- **Truncation:** Symmetric truncated normal at ±`n_sigma · σ` (constructor parameter, default 3, settable via `[Visini2025SecondaryFD] n_sigma = <value>`; the legacy name `truncation_eps` is still accepted); lower clamp at `1e-16` to avoid `log(0)`.
- **Epistemic branches:** combinations A / B / C as additive offsets `g` ([lines 73–76](../../openquake/fdha/secondary_surf_displ/visini2025.py#L73-L76)); style; HW/FW.

---

## 5. Cross-Cutting Patterns

The table below summarises the variability treatment across the displacement
models (rupture models are uniformly logistic with no aleatory σ on `P`).

| Pattern | Models that use it |
| :--- | :--- |
| **Lognormal in log₁₀ via AD/MD scaling, convolved with Gamma(D/AD) and/or Beta(D/MD)** | Youngs 2003, Moss & Ross 2011, Takao 2013, Moss 2022, Moss 2024 (PSD); Youngs 2003, Moss 2022 (SSD) |
| **Lognormal in natural-log `ln D` directly** | Petersen 2011 (PSD and SSD); Visini 2025 (SSD) |
| **Normal on power-transformed `D` (Box–Cox-like)** | Lavrentiadis 2023 (`D^0.3`); Kuehn 2024 (`(D^λ − 1)/λ`) |
| **Asymmetric (nEMG)** | Chiou 2025 |
| **Explicit between-event / within-event (τ/φ) decomposition** | Lavrentiadis 2023, Kuehn 2024 (and Chiou 2025 with magnitude/position split) |
| **Explicit discrete zero-displacement probability inside the model** | Lavrentiadis 2023 only (P_gap, P_zero_slip) |
| **Ensemble-of-coefficients epistemic propagation inside `get_prob`** | Kuehn 2024 (with `epistemic_uncertainty=True`) |
| **Epistemic exposed as named branches (style / dataset / version / σ source / method)** | All other models |
| **Truncation at ±`n_sigma`σ (user-configurable; default 3)** | Moss & Ross 2011 (PSD); Takao 2013 (PSD); Youngs 2003 (SSD); Visini 2025 (SSD) |
| **Truncation at ±`n_sigma`σ ε-space integration (user-configurable; default 6)** | Youngs 2003 (PSD) |
| **Truncation at ±6σ ε-space integration (fixed)** | Moss 2022 / 2024 (PSD and Gamma-method SSD) |
| **No truncation** | Petersen 2011 (PSD and SSD); Chiou 2025; Lavrentiadis 2023; Kuehn 2024 |

### Practical guidance for users

- **Truncation level is user-configurable for the truncated models.** Moss & Ross
  2011, Takao 2013, Youngs 2003 (PSD), and Visini 2025 (SSD) expose the truncation
  half-width as an `n_sigma` constructor parameter (defaults: 3 for Moss & Ross /
  Takao / Visini, 6 for Youngs PSD). Set it from a logic-tree branch, e.g.:
  ```xml
  <uncertaintyModel><![CDATA[[MossRoss2011PrimaryFD]
  n_sigma = 4]]></uncertaintyModel>
  ```
  (For backward compatibility, Visini 2025 SSD also accepts the former name
  `truncation_eps`.)
- **σ is mostly fixed in log space.** Apart from Lavrentiadis 2023, Kuehn 2024,
  and Chiou 2025, the aleatory σ in this library is a scalar regression value
  (commonly the Wells & Coppersmith 1994 values 0.36 / 0.42 in log₁₀) without
  magnitude or distance dependence. Sensitivity studies that vary σ should
  recognise this when comparing models.
- **τ/φ-aware logic trees only make sense for three models.** If your hazard
  framework treats between-event and within-event variability separately (e.g.
  to support partial-correlation aggregation across sites), only Lavrentiadis
  2023, Kuehn 2024, and — with a magnitude/position rather than between/within
  reading — Chiou 2025 expose the necessary components. The other PSD models
  must be used with a single total σ.
- **Zero-displacement / no-slip handling.** Lavrentiadis 2023 is the only
  model that resolves "no displacement at this site even though the rupture
  occurred" inside the displacement model itself. With the other PSD models,
  this branch is the responsibility of the upstream PSR model (`P(SR|M)`) and
  the truncation of the conditional distribution.
- **Epistemic branching is mostly user-driven.** With the exception of
  Kuehn 2024's posterior-sample ensemble, every model presents one curve per
  parameter setting. Logic trees should enumerate the named alternatives
  documented above (style, dataset version, completeness, σ source, method,
  cell/pixel size, HW/FW), each with an explicit weight.
- **The principal/distributed distance threshold can itself be an epistemic
  branch.** The hard-step split at `r_threshold_km` is a simplification of the
  rupture-location term fr(r) of Petersen et al. (2011, p. 810); the toolkit
  exposes the threshold choice as weighted `fdhaCalcRThreshold` logic-tree
  branches (each `<uncertaintyModel>` a bare value in km), following
  IAEA-TECDOC-2092 (2025, §3.3). This is mutually exclusive with the scalar
  `[calculation].r_threshold_km` INI key — see the conflict rule in the
  [Configuration](05-Configuration.md) chapter.
- **Rupture-probability models carry no quantified σ.** When PSR or SSR
  dominate the hazard (typically at moderate magnitudes and long return
  periods, as noted in the context notes of the [Models](06-Models.md) chapter), epistemic
  spread is the only available knob — populate multiple PSR/SSR branches
  with justified weights rather than expecting an aleatory σ on `P` to
  capture model uncertainty.
