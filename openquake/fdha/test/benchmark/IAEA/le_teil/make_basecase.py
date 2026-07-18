#!/usr/bin/env python
"""Generate the Le Teil BASE-CASE epistemic tree from the author workbook.

The exercise workbook (``PFDHA Benchmarking Study - Le Teil - Input
data.xlsx``, sheets "Base Case" + "Moment balancing") defines the Le Teil
base case as a logic tree over the 21-km La Rouviere source:

    seismogenic thickness   2 km (w 0.4) -> rupture length 2.8 km (w 0.2)
                                            rupture length 4.6 km (w 0.8)
                            4 km (w 0.4) -> rupture length 5.8 km (w 1.0)
                           10 km (w 0.2) -> rupture length 21 km (w 1.0)
    magnitude               three per geometry, weights 0.3 / 0.4 / 0.3
    slip rate               0.015 / 0.01 / 0.005 mm/yr, weights 0.1/0.4/0.5

The "Moment balancing" sheet gives, for each of the 36 end branches, the
moment-balanced event rate times the probability that the rupture covers
the site ("Rate, rups (1/yr)" - the same site-covering-rate convention as
the single-branch jobs, whose 4.6e-5 is its slip-0.01 / M5.52 row).

This script collapses the slip-rate level (exact for the mean hazard
curve, which is linear in the branch rates) and emits

- ``source_model_basecase_<i>.xml``  (12 = 4 geometries x 3 magnitudes)
- ``smlt_basecase.xml``              (weights = w_geometry x w_magnitude)

Each branch trace runs along the LRF alignment of ``source_model_lrf.xml``
with the principal site kept at the same relative along-strike position
(l/L = 0.46) as in the single-branch job, and lower seismogenic depth
equal to the branch thickness (dip 45).

The generated files are committed; rerun only to regenerate them
(requires the author workbook)::

    PYTHONPATH=. python make_basecase.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORKBOOK = ("/home/ychen/Documents/PFDHA/IAEA_benchmarking/"
            "Input for Modeling-20260319T133033Z-3-001/Input for Modeling/"
            "Le Teil Case - Reverse/"
            "PFDHA Benchmarking Study - Le Teil - Input data.xlsx")

# LRF trace of the single-branch job (SW -> NE) and the principal site
P0 = np.array([4.64198, 44.51204])
P1 = np.array([4.69428, 44.54868])
SITE = np.array([4.67, 44.531])

SLIP_WEIGHTS = {0.015: 0.1, 0.01: 0.4, 0.005: 0.5}
GEOMETRY_WEIGHTS = {  # (rupture length km, thickness km) -> weight
    (2.83, 2.0): 0.4 * 0.2,
    (4.6, 2.0): 0.4 * 0.8,
    (5.8, 4.0): 0.4 * 1.0,
    (21.0, 10.0): 0.2 * 1.0,
}
MAG_WEIGHTS = (0.3, 0.4, 0.3)  # per geometry, highest magnitude first

SM_TEMPLATE = """<?xml version="1.0" encoding="utf-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.5"
      xmlns:gml="http://www.opengis.net/gml">
  <sourceModel name="IAEA Le Teil base case: {name}">
    <sourceGroup name="IAEA Le Teil base case" tectonicRegion="Active Shallow Crust">
      <characteristicFaultSource id="LRF_{tag}" name="La Rouviere Fault {name}" tectonicRegion="Active Shallow Crust">
        <incrementalMFD binWidth="0.1" minMag="{mag}">
          <occurRates>{rate:.6e}</occurRates>
        </incrementalMFD>
        <rake>90.0</rake>
        <surface>
          <simpleFaultGeometry>
            <gml:LineString>
              <gml:posList>
                {x0:.5f} {y0:.5f}
                {x1:.5f} {y1:.5f}
              </gml:posList>
            </gml:LineString>
            <dip>45</dip>
            <upperSeismoDepth>0.0</upperSeismoDepth>
            <lowerSeismoDepth>{thickness}</lowerSeismoDepth>
          </simpleFaultGeometry>
        </surface>
      </characteristicFaultSource>
    </sourceGroup>
  </sourceModel>
</nrml>
"""


def read_branches():
    """36 moment-balancing rows -> 12 (L_rup, thickness, mag) branches with
    slip-collapsed site-covering rates."""
    import openpyxl

    ws = openpyxl.load_workbook(WORKBOOK, data_only=True)["Moment balancing"]
    rows = []
    for row in ws.iter_rows(min_row=1, values_only=True):
        vals = [v for v in row if v is not None]
        if len(vals) >= 14 and isinstance(vals[0], (int, float)):
            (_, _, width, _, slip, _, _, mag, _, _, l_rup, _, rate_site,
             _) = vals[:14]
            thickness = round(float(width) * np.sin(np.radians(45)), 0)
            rows.append((round(float(l_rup), 2), thickness, float(mag),
                         float(slip), float(rate_site)))
    assert len(rows) == 36, f"expected 36 end branches, got {len(rows)}"

    branches = {}
    for l_rup, thickness, mag, slip, rate in rows:
        key = (l_rup, thickness, mag)
        branches.setdefault(key, 0.0)
        branches[key] += SLIP_WEIGHTS[slip] * rate
    assert len(branches) == 12, sorted(branches)
    return branches


def make_trace(l_rup_km):
    """Trace of length l_rup along the LRF alignment, site at l/L = 0.46."""
    seg = P1 - P0
    # local metric: km per degree along the segment direction
    km_per_deg = np.hypot(seg[0] * np.cos(np.radians(SITE[1])) * 111.195,
                          seg[1] * 111.195) / np.linalg.norm(seg)
    u = seg / np.linalg.norm(seg)  # unit vector in degree space
    # site's projection onto the P0->P1 line, in km from P0
    s_site = np.dot(SITE - P0, seg) / np.linalg.norm(seg) * km_per_deg
    frac = s_site / (np.linalg.norm(seg) * km_per_deg)  # ~0.46
    start = P0 + u * (s_site - frac * l_rup_km) / km_per_deg
    end = P0 + u * (s_site + (1 - frac) * l_rup_km) / km_per_deg
    return start, end


def main():
    branches = read_branches()
    entries = []
    for i, ((l_rup, thickness, mag), rate) in enumerate(
            sorted(branches.items())):
        geom_w = GEOMETRY_WEIGHTS[(l_rup, thickness)]
        mags_in_geom = sorted(
            (m for (l, t, m) in branches if (l, t) == (l_rup, thickness)),
            reverse=True)
        mag_w = MAG_WEIGHTS[mags_in_geom.index(mag)]
        (x0, y0), (x1, y1) = make_trace(l_rup)
        tag = f"L{l_rup:g}_M{mag:g}".replace(".", "p")
        name = f"L={l_rup:g} km, T={thickness:g} km, Mw{mag:g}"
        fname = f"source_model_basecase_{tag}.xml"
        (HERE / fname).write_text(SM_TEMPLATE.format(
            name=name, tag=tag, mag=mag, rate=rate, thickness=thickness,
            x0=x0, y0=y0, x1=x1, y1=y1))
        entries.append((fname, geom_w * mag_w, name))
        print(f"{fname}: w={geom_w * mag_w:.3f} rate={rate:.4e}")

    total = sum(w for _, w, _ in entries)
    assert abs(total - 1.0) < 1e-9, total

    lines = ['<?xml version="1.0" encoding="UTF-8"?>',
             '<nrml xmlns:gml="http://www.opengis.net/gml"',
             '      xmlns="http://openquake.org/xmlns/nrml/0.4">',
             '  <logicTree logicTreeID="lt_sm_le_teil_basecase">',
             '    <logicTreeBranchSet uncertaintyType="sourceModel"'
             ' branchSetID="bs_sm_le_teil_basecase">']
    for i, (fname, w, name) in enumerate(entries):
        lines += [f'      <logicTreeBranch branchID="b_basecase_{i:02d}">',
                  f'        <!-- {name} -->',
                  f'        <uncertaintyModel>./{fname}</uncertaintyModel>',
                  f'        <uncertaintyWeight>{w:.3f}</uncertaintyWeight>',
                  '      </logicTreeBranch>']
    lines += ['    </logicTreeBranchSet>', '  </logicTree>', '</nrml>', '']
    (HERE / "smlt_basecase.xml").write_text("\n".join(lines))
    print(f"smlt_basecase.xml: {len(entries)} branches, total weight {total}")


if __name__ == "__main__":
    main()
