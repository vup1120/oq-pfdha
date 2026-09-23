#!/usr/bin/env python
"""Why the published Le Teil M11 curve cannot be the single-M5.5 scenario.

Diagnostic established against Moss & Ross (2011, BSSA 101(4), 1542-1553):

1. Our implementation matches the paper exactly: ``MossRoss2011PrimarySR``
   reproduces Eq. 5 (a=7.30, b=-1.03) to 4 decimals, and
   ``MossRoss2011PrimaryFD`` carries Eqs. 7-9 verbatim (D/AD gamma
   coefficients; log AD = 0.3244 M - 2.2192, sigma 0.17; log MD =
   0.5102 M - 3.1971, sigma 0.31).

2. The published M11 curve (reference/fig4c_leteil_principal_M11.csv),
   normalized by its head, is an almost perfect lognormal with median
   ~1.4 m and sigma_ln ~0.89. Through the paper's own Eq. 9 a median MD of
   1.4 m corresponds to M ~ 6.6 -- not 5.5, where the paper's scaling caps
   the median at AD = 0.37 m / MD = 0.41 m. At M5.5 the model itself gives
   P(D > 10 m) ~ 5e-7, while the published curve carries 1e-2 of its head
   at 10 m: no parameterization (AD/MD, gamma/Weibull, any truncation, any
   P_sr) can produce that tail at M5.5.

3. The Le Teil author workbook's base-case logic tree contains a
   10-km-seismogenic-thickness alternative (weight 0.2, rupture length
   21 km) with magnitudes 6.23/6.48/6.73. Running our Moss & Ross chain on
   those branches alone (their own Eq. 5 P_sr, AD method, x/L = 0.46)
   lands within a factor ~1.7 of the published head and crosses the
   published curve near 1 m -- the right ballpark, unlike the 8x-high /
   20 000x-thin-tailed single-M5.5 chain.

Conclusion: the M11 team's submission is evidently dominated by the
large-magnitude thickness alternative of the base-case tree (consistent
with Moss & Ross's own methodology of integrating over a magnitude
distribution, cf. their Los Osos example), not the single-M5.5 scenario
shown for the other teams in the paper's Fig. 4c. The paper's text ("Moss
et al. 2013, 10% at M5.5") is itself inconsistent with the published
head (0.0122 x rate). Absent the team's actual run configuration the
entry stays qualitative.

Writes ``../Figures/iaea_leteil_M11_diagnosis.png``.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("NUMBA_DISABLE_JIT", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from openquake.pfd.primary_surf_rup import (Moss2013PrimarySR,
                                             MossRoss2011PrimarySR)
from openquake.pfd.primary_surf_displ import MossRoss2011PrimaryFD

HERE = Path(__file__).resolve().parent
RATE_M55 = 4.6e-5
XL = 0.46

# Thickness-10 km alternative of the Le Teil base-case tree (author workbook,
# "Moment balancing" sheet): branch weight 0.2, rupture length 21 km,
# magnitudes weighted .3/.4/.3, rates averaged over the slip-rate branches.
TH10 = [(6.73, 0.3, 5.11e-6), (6.48, 0.4, 1.21e-5), (6.23, 0.3, 2.87e-5)]
W_TH10 = 0.2


def main() -> None:
    ref_rows = list(csv.reader(
        (HERE.parent / "reference" / "fig4c_leteil_principal_M11.csv").open()))
    d_m = np.array([float(r[0]) for r in ref_rows[1:]]) / 100.0
    ref = np.array([float(r[1]) for r in ref_rows[1:]])

    fd = MossRoss2011PrimaryFD(n_sigma=5.0)
    m13 = Moss2013PrimarySR()
    mr11 = MossRoss2011PrimarySR()

    # (a) the nominal single-M5.5 exercise scenario (as in job_principal_M11)
    p55 = np.asarray(fd.get_prob(d=d_m, X_L_ratio=XL, mag=5.5,
                                 norm_disp_type="AD")).reshape(len(d_m), -1)[:, 0]
    single = RATE_M55 * m13.get_prob(mag=5.5, style="reverse", vs30=760.0) * p55

    # (b) reconstruction: thickness-10 branches, team's own Eq. 5 P_sr, AD
    recon = np.zeros_like(d_m)
    for mag, w_m, rate in TH10:
        p = np.asarray(fd.get_prob(d=d_m, X_L_ratio=XL, mag=mag,
                                   norm_disp_type="AD")).reshape(len(d_m), -1)[:, 0]
        recon += W_TH10 * w_m * rate * float(mr11.get_prob(mag=mag)) * p

    fig, ax = plt.subplots(figsize=(8.5, 6))
    ax.loglog(d_m * 100, ref, "o-", color="#17becf", lw=2.5, ms=5,
              label="M11 published (paper Fig. 4c)")
    ax.loglog(d_m * 100, single, "--", color="#d62728", lw=2,
              label="single M5.5 scenario (Moss13 P$_{sr}$=0.10 + MossRoss AD)")
    ax.loglog(d_m * 100, recon, "--", color="#2ca02c", lw=2,
              label="reconstruction: 10-km-thickness branch\n"
                    "(M6.2-6.7, w=0.2, Eq. 5 P$_{sr}$, MossRoss AD)")
    ax.set_xlabel("Principal displacement (cm)")
    ax.set_ylabel("AFOE (yr$^{-1}$)")
    ax.set_ylim(1e-10, 1e-5)
    ax.set_title("Le Teil M11 diagnosis: the published curve implies M$\\approx$6.6\n"
                 "(MD median 1.4 m), not the single-M5.5 scenario")
    ax.grid(True, which="both", ls=":", alpha=0.35)
    ax.legend(fontsize=8.5, loc="lower left")
    fig.tight_layout()
    out = HERE.parent / "Figures" / "iaea_leteil_M11_diagnosis.png"
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=150)
    print("saved", out)


if __name__ == "__main__":
    main()
