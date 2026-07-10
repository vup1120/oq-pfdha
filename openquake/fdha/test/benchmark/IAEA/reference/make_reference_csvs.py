#!/usr/bin/env python
"""Write the digitized IAEA-paper hazard curves as tidy CSV files.

The vectors below are transcribed verbatim from the MATLAB plotting scripts
provided by the exercise coordinators (``Figure4 (1).m`` and
``Figure6 (2).m``, A. Valentini, IAEA), which contain the team-supplied
hazard curves shown in Figs 4 and 6 of Valentini et al. (IAEA PFDHA
exercise paper; the same curves appear as Figs 14/18/21 and 15/19/22 of
IAEA TECDOC-2092).

Displacements are in cm, annual frequencies of exceedance in 1/yr.

Model abbreviations (TECDOC-2092 Table 8):
  P11 Petersen et al. (2011)      T13 Takao et al. (2013/2014/2016)
  C24 Chiou et al. (2025)         K24 Kuehn et al. (2024)
  L23 Lavrentiadis & Abrahamson   M11 Moss & Ross (2011)
  Y03 Youngs et al. (2003)        V24 Visini et al. (2025)

Run once to (re)generate the CSVs next to this script:

    python make_reference_csvs.py
"""

from pathlib import Path

HERE = Path(__file__).resolve().parent

# 18-point displacement grid shared by all curves except M11 (cm)
DISP_CM = [0.01, 0.1, 0.5, 1.0, 1.5, 3.0, 5.0, 7.5, 10.0, 15.0,
           30.0, 50.0, 75.0, 100.0, 300.0, 500.0, 750.0, 1000.0]

# ----------------------------------------------------------- Figure 4 (principal)
FIG4A_KUMAMOTO = {  # single-segment (Uto) case
    "P11": [1.32e-04, 1.32e-04, 1.32e-04, 1.31e-04, 1.30e-04, 1.25e-04,
            1.18e-04, 1.10e-04, 1.03e-04, 8.95e-05, 6.31e-05, 4.35e-05,
            3.00e-05, 2.20e-05, 4.64e-06, 1.82e-06, 7.85e-07, 4.09e-07],
    "T13": [4.13e-05, 4.13e-05, 4.13e-05, 4.11e-05, 4.09e-05, 4.00e-05,
            3.85e-05, 3.63e-05, 3.42e-05, 3.01e-05, 2.08e-05, 1.36e-05,
            8.56e-06, 5.75e-06, 6.80e-07, 1.74e-07, 4.94e-08, 1.83e-08],
    "L23": [1.46e-04, 1.44e-04, 1.38e-04, 1.33e-04, 1.29e-04, 1.20e-04,
            1.09e-04, 9.81e-05, 8.88e-05, 7.36e-05, 4.42e-05, 2.40e-05,
            1.19e-05, 6.20e-06, 7.92e-08, 2.16e-09, 4.08e-11, 9.90e-13],
    "K24": [1.88932e-04, 1.88604e-04, 1.87097e-04, 1.85005e-04, 1.82779e-04,
            1.75756e-04, 1.66281e-04, 1.54888e-04, 1.44249e-04, 1.25382e-04,
            8.46014e-05, 5.31141e-05, 3.18888e-05, 2.03010e-05, 1.60788e-06,
            3.10767e-07, 6.98637e-08, 2.23468e-08],
    "C24": [1.31775e-04, 1.31572e-04, 1.30400e-04, 1.28706e-04, 1.26877e-04,
            1.20965e-04, 1.12737e-04, 1.02679e-04, 9.33119e-05, 7.71261e-05,
            4.54213e-05, 2.48653e-05, 1.32235e-05, 7.71793e-06, 4.62204e-07,
            8.05555e-08, 1.63559e-08, 4.71018e-09],
}

FIG4B_KUMAMOTO_FLOATING = {
    "T13": [5.03e-05, 5.03e-05, 5.02e-05, 5.01e-05, 4.98e-05, 4.89e-05,
            4.72e-05, 4.48e-05, 4.24e-05, 3.78e-05, 2.72e-05, 1.84e-05,
            1.21e-05, 8.43e-06, 1.17e-06, 3.30e-07, 1.01e-07, 3.97e-08],
    "L23": [2.20e-04, 2.15e-04, 2.06e-04, 1.99e-04, 1.93e-04, 1.78e-04,
            1.62e-04, 1.45e-04, 1.32e-04, 1.09e-04, 6.62e-05, 3.64e-05,
            1.84e-05, 9.72e-06, 1.38e-07, 4.07e-09, 8.35e-11, 2.28e-12],
    "K24": [2.06799e-04, 2.05803e-04, 2.02627e-04, 1.99162e-04, 1.95924e-04,
            1.86983e-04, 1.76248e-04, 1.64199e-04, 1.53356e-04, 1.34536e-04,
            9.38583e-05, 6.12689e-05, 3.81821e-05, 2.49982e-05, 2.17364e-06,
            4.30250e-07, 9.73976e-08, 3.11132e-08],
}

FIG4C_LETEIL = {
    "T13": [1.35e-07, 1.35e-07, 1.31e-07, 1.25e-07, 1.18e-07, 1.00e-07,
            8.11e-08, 6.35e-08, 5.08e-08, 3.43e-08, 1.35e-08, 5.38e-09,
            2.21e-09, 1.07e-09, 3.23e-11, 4.15e-12, 6.64e-13, 1.62e-13],
    "K24": [4.59997e-05, 4.59930e-05, 4.58948e-05, 4.56509e-05, 4.53041e-05,
            4.38557e-05, 4.14256e-05, 3.81448e-05, 3.49149e-05, 2.90803e-05,
            1.70429e-05, 8.99485e-06, 4.47966e-06, 2.42664e-06, 7.97272e-08,
            8.44388e-09, 1.03261e-09, 1.95757e-10],
    "L23": [2.75e-05, 2.74e-05, 2.72e-05, 2.70e-05, 2.67e-05, 2.58e-05,
            2.45e-05, 2.28e-05, 2.11e-05, 1.78e-05, 1.03e-05, 4.74e-06,
            1.80e-06, 6.93e-07, 6.72e-10, 1.55e-12, 0.0, 0.0],
}

# Moss & Ross (2011) team curve is tabulated on its own displacement grid
DISP_MOSS_CM = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80,
                90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
FIG4C_LETEIL_M11 = [5.613e-07, 5.611e-07, 5.607e-07, 5.603e-07, 5.598e-07,
                    5.592e-07, 5.586e-07, 5.579e-07, 5.571e-07, 5.563e-07,
                    5.454e-07, 5.308e-07, 5.136e-07, 4.949e-07, 4.751e-07,
                    4.548e-07, 4.344e-07, 4.142e-07, 3.943e-07, 2.322e-07,
                    1.354e-07, 8.05e-08, 4.90e-08, 3.06e-08, 1.96e-08,
                    1.28e-08, 8.5e-09, 5.8e-09]

FIG4D_NORCIA = {
    "L23": [3.69e-04, 3.64e-04, 3.54e-04, 3.45e-04, 3.38e-04, 3.21e-04,
            3.02e-04, 2.82e-04, 2.65e-04, 2.36e-04, 1.74e-04, 1.23e-04,
            8.33e-05, 5.86e-05, 6.26e-06, 1.05e-06, 1.50e-07, 2.58e-08],
    "K24": [4.02941e-04, 4.00644e-04, 3.92324e-04, 3.83067e-04, 3.74595e-04,
            3.52434e-04, 3.28041e-04, 3.02942e-04, 2.81947e-04, 2.48150e-04,
            1.81946e-04, 1.31022e-04, 9.30131e-05, 6.89562e-05, 1.20247e-05,
            3.52947e-06, 1.09417e-06, 4.34675e-07],
    "Y03": [3.29e-04, 3.25e-04, 3.16e-04, 3.08e-04, 3.01e-04, 2.85e-04,
            2.67e-04, 2.49e-04, 2.34e-04, 2.08e-04, 1.56e-04, 1.14e-04,
            8.16e-05, 6.15e-05, 1.31e-05, 4.72e-06, 1.79e-06, 8.09e-07],
}

# --------------------------------------------------------- Figure 6 (distributed)
FIG6A_KUMAMOTO_DIST = {
    "T13": [1.67e-08, 1.67e-08, 1.64e-08, 1.59e-08, 1.54e-08, 1.40e-08,
            1.24e-08, 1.06e-08, 9.22e-09, 7.12e-09, 3.75e-09, 1.92e-09,
            9.90e-10, 5.70e-10, 3.55e-11, 6.57e-12, 1.43e-12, 4.34e-13],
    "P11": [2.95e-07, 2.95e-07, 2.72e-07, 2.32e-07, 1.96e-07, 1.25e-07,
            7.61e-08, 4.59e-08, 3.01e-08, 1.52e-08, 3.60e-09, 1.00e-09,
            3.16e-10, 1.30e-10, 2.42e-12, 2.79e-13, 4.34e-14, 1.07e-14],
}

FIG6B_LETEIL_DIST = {
    "T13": [9.31e-10, 9.26e-10, 8.90e-10, 8.40e-10, 7.91e-10, 6.62e-10,
            5.31e-10, 4.14e-10, 3.31e-10, 2.23e-10, 8.90e-11, 3.59e-11,
            1.50e-11, 7.35e-12, 2.36e-13, 3.13e-14, 5.16e-15, 1.28e-15],
    "V24": [2.57e-09, 2.57e-09, 2.57e-09, 2.57e-09, 2.57e-09, 2.56e-09,
            2.44e-09, 2.05e-09, 1.57e-09, 1.04e-09, 1.98e-10, 1.30e-11,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
}

FIG6C_NORCIA_DIST = {
    "Y03": [1.21e-06, 1.21e-06, 1.19e-06, 1.16e-06, 1.12e-06, 9.92e-07,
            8.40e-07, 6.91e-07, 5.77e-07, 4.19e-07, 1.96e-07, 9.27e-08,
            4.48e-08, 2.47e-08, 1.17e-09, 1.48e-10, 1.73e-11, 2.51e-12],
    "V24": [2.98e-07, 2.98e-07, 2.98e-07, 2.98e-07, 2.98e-07, 2.98e-07,
            2.95e-07, 2.79e-07, 2.52e-07, 2.10e-07, 8.78e-08, 2.61e-08,
            5.69e-09, 9.56e-10, 1.1046719e-20, 0.0, 0.0, 0.0],
}


def write_csv(name: str, disp_cm, curves: dict) -> None:
    path = HERE / name
    cols = list(curves)
    lines = ["disp_cm," + ",".join(cols)]
    for i, d in enumerate(disp_cm):
        lines.append(f"{d:g}," + ",".join(f"{curves[c][i]:.6e}" for c in cols))
    path.write_text("\n".join(lines) + "\n")
    print(f"wrote {path.name}: {len(disp_cm)} rows x {len(cols)} models")


def main() -> None:
    for curves in (FIG4A_KUMAMOTO, FIG4B_KUMAMOTO_FLOATING, FIG4C_LETEIL,
                   FIG4D_NORCIA, FIG6A_KUMAMOTO_DIST, FIG6B_LETEIL_DIST,
                   FIG6C_NORCIA_DIST):
        for name, vals in curves.items():
            assert len(vals) == 18, (name, len(vals))
    assert len(FIG4C_LETEIL_M11) == len(DISP_MOSS_CM) == 28

    write_csv("fig4a_kumamoto_principal.csv", DISP_CM, FIG4A_KUMAMOTO)
    write_csv("fig4b_kumamoto_principal_floating.csv", DISP_CM,
              FIG4B_KUMAMOTO_FLOATING)
    write_csv("fig4c_leteil_principal.csv", DISP_CM, FIG4C_LETEIL)
    write_csv("fig4c_leteil_principal_M11.csv", DISP_MOSS_CM,
              {"M11": FIG4C_LETEIL_M11})
    write_csv("fig4d_norcia_principal.csv", DISP_CM, FIG4D_NORCIA)
    write_csv("fig6a_kumamoto_distributed.csv", DISP_CM, FIG6A_KUMAMOTO_DIST)
    write_csv("fig6b_leteil_distributed.csv", DISP_CM, FIG6B_LETEIL_DIST)
    write_csv("fig6c_norcia_distributed.csv", DISP_CM, FIG6C_NORCIA_DIST)


if __name__ == "__main__":
    main()
