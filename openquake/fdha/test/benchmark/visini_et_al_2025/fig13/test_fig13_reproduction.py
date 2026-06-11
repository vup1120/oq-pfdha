"""
Benchmark: reproduction of Visini et al. (2025) Figure 13.

Runs the three decision-tree example cases of the paper end-to-end through
the ``fdha`` CLI and compares the aggregate conditional-probability curves
against the curves digitized from the published Figure 13
(``reference_data/visini2025_case{1,2,3}.csv``).

Scenario (see fig13/README.md for the full derivation):
- Mw 7.0 normal-faulting characteristic rupture, 40 km long, dip 36 N;
- 100 m x 100 m site on the hanging wall, 2.0 km from the trace, x/L = 0.5;
- TPFm = 1.42 m (value stated in the paper);
- n_sigma = 3 with truncation-range normalization (FDHLab eps = 3);
- Case 1 = Combinations A+B+C, Case 2 = A+B, Case 3 = A.

Tolerances account for (a) digitization error of the published log-log
curves and (b) the Monte Carlo noise of the along-strike occurrence factor
(about +/-5 percent per run, not seeded).
"""

import csv
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

pytestmark = [pytest.mark.benchmark, pytest.mark.visini2025, pytest.mark.slow]

BASE = Path(__file__).parent

# Compare only inside the digitized range, above the digitization floor
# (1e-5) and below the +3-sigma truncation cliff of Combination A/B (~2.5 m).
D_MAX = 2.0
POINT_RATIO_BOUNDS = (0.60, 1.50)
MEDIAN_RATIO_BOUNDS = (0.80, 1.25)


def _run_case(case):
    out_json = BASE / f"_pytest_case{case}.json"
    res = subprocess.run(
        [sys.executable, "-m", "openquake.fdha.main",
         str(BASE / f"job_case{case}.ini"), "--output", str(out_json)],
        capture_output=True, text=True, cwd=str(BASE), timeout=900,
    )
    assert res.returncode == 0, f"case{case} run failed:\n{res.stderr[-2000:]}"
    result = json.loads(out_json.read_text())
    out_json.unlink()
    agg = Path(result["aggregate_hazard_csv"])
    rows = np.genfromtxt(agg, delimiter=",", skip_header=1)
    return rows[:, 0], rows[:, 1]  # displacement (m), mean probability


def _load_reference(case):
    rows = []
    with open(BASE / "reference_data" / f"visini2025_case{case}.csv") as f:
        for row in csv.reader(f):
            try:
                rows.append((float(row[0]), float(row[1])))
            except (ValueError, IndexError):
                continue  # header or blank line
    rows.sort()
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1]


@pytest.mark.parametrize("case", [1, 2, 3])
def test_fig13_case_reproduces_published_curve(case):
    d, p = _run_case(case)
    dr, pr = _load_reference(case)

    mask = (d >= dr.min()) & (d <= min(D_MAX, dr.max())) & (p > 0)
    assert mask.sum() >= 6, "too few comparable displacement levels"

    ref_interp = np.exp(
        np.interp(np.log(d[mask]), np.log(dr), np.log(np.maximum(pr, 1e-12)))
    )
    ratio = p[mask] / ref_interp

    lo, hi = POINT_RATIO_BOUNDS
    assert np.all((ratio >= lo) & (ratio <= hi)), (
        f"case{case}: point ratios outside [{lo}, {hi}]: "
        f"{[(float(x), float(rt)) for x, rt in zip(d[mask], ratio)]}"
    )
    mlo, mhi = MEDIAN_RATIO_BOUNDS
    med = float(np.median(ratio))
    assert mlo <= med <= mhi, f"case{case}: median ratio {med:.3f} outside [{mlo}, {mhi}]"
