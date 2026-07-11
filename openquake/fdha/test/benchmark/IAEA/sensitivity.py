# -*- coding: utf-8 -*-
"""Registry of the IAEA-exercise SENSITIVITY-CASE demonstration jobs.

These implement the author workbooks' sensitivity cases, for which the
exercise published no reference hazard curves. They are validated by

1. structural checks (positive head, non-increasing hazard curve), and
2. regression snapshots committed in ``reference_snapshots/`` —
   regenerate with ``python run_sensitivity.py --update-snapshots`` after
   an intentional model change, and review the diff.

Inputs are generated from the author workbooks + shapefiles by
``make_sensitivity_cases.py``; see that module's docstring for the case
definitions.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SensitivityJob:
    case: str
    job: str
    description: str


SENSITIVITY_JOBS = [
    # Suizenji is a Mw 5.8 source: the Petersen (2011) models hard-fail
    # below M6 (outside their calibration), so the Kumamoto sens1/sens3
    # demonstrations use the Takao (2013) chains (Japanese strike-slip
    # data, calibrated from Mw 5.8 up).
    SensitivityJob("kumamoto", "distributed_T13_sens1",
                   "Suizenji fault (M5.8), site r = 0.6 km"),
    SensitivityJob("kumamoto", "principal_T13_sens3",
                   "Suizenji fault, on-fault site"),
    SensitivityJob("kumamoto", "distributed_P11_sens4",
                   "four rupture sources (middle branches), site r = 10 km"),
    SensitivityJob("le_teil", "distributed_T13_sens1",
                   "four fault-source combinations (w 0.25), site r = 0.6 km"),
    SensitivityJob("le_teil", "distributed_T13_sens4",
                   "LRF_2 + MRF_3 + PCF_2 summed, site r = 0.6 km"),
    SensitivityJob("norcia", "distributed_Y03_sens2",
                   "MVFS at the r = 2.4 km site"),
    SensitivityJob("norcia", "distributed_Y03_sens3",
                   "MVFS + NFS (GR-like MFD), site 7.6 km HW / 3 km FW"),
]
