#!/usr/bin/env python
"""Generate the Kumamoto BASE-CASE source model from the author workbook.

The exercise workbook (``PFDHA Benchmarking Study - Kumamoto - Input
data.xlsx``, sheet "Base Case") defines the Kumamoto base case at the
r = 5.22 km site as FOUR COEXISTING rupture sources on the
Futagawa-Uto-Uto-Hanto-North system (a rupture-scenario budget: the
single-segment rupture is frequent, multi-segment ruptures rare), each
with independent epistemic magnitude and rate branches (weights
0.2/0.6/0.2 both):

    Uto (2)                       L 22 km, dip 60 NW,  Mw 6.3/6.5/6.7,
                                  rate 6.3/18.9/56.7 e-5
    Futagawa (1) + Uto (2)        L 46 km, dip 68.6 NW, Mw 6.7/6.9/7.1,
                                  rate 0.427/1.28/3.84 e-5
    Uto (2) + Uto-Hanto-North (3) L 54 km, dip 60 NW,  Mw 6.8/7.0/7.2,
                                  rate 1.177/3.53/10.59 e-5
    Futagawa+Uto+Uto-Hanto-North  L 78 km, dip 64.6 NW, Mw 7.0/7.2/7.4,
                                  rate 0.427/1.28/3.84 e-5

Because the mean hazard curve is linear in the branch rates, the rate
level collapses exactly into the weighted-mean rate, and the three
magnitude branches of a source are represented exactly (for the mean) as
a three-bin incremental MFD with bin rates = w_mag x mean rate. The
result is ONE source model with four characteristicFaultSources
(``source_model_basecase.xml``) and a single-branch
``smlt_basecase.xml``.

Traces are read from the workbook "Coordinates" sheet (the combined
traces are listed explicitly); vertex order is NE -> SW so the surfaces
dip NW, matching ``source_model_uto.xml``. Rupture thickness 14 km ->
lowerSeismoDepth 14, rake 0 (strike-slip), as in the single-branch job.

Rerun only to regenerate the committed files (requires the workbook)::

    PYTHONPATH=. python make_basecase.py
"""

from __future__ import annotations

from pathlib import Path

HERE = Path(__file__).resolve().parent
WORKBOOK = ("/home/ychen/Documents/PFDHA/IAEA_benchmarking/"
            "Input for Modeling-20260319T133033Z-3-001/Input for Modeling/"
            "Kumamoto Case - Strike Slip/"
            "PFDHA Benchmarking Study - Kumamoto - Input data.xlsx")

MAG_WEIGHTS = (0.2, 0.6, 0.2)
RATE_WEIGHTS = (0.2, 0.6, 0.2)

# Epistemic table of the "Base Case" sheet (magnitudes and mean rates in
# descending order; both weighted 0.2/0.6/0.2). Rates in 1e-5 /yr.
SOURCES = [
    # (coordinates-sheet key, id, dip, mags desc, rates desc)
    ("UTO (2)", "UTO", 60.0, (6.7, 6.5, 6.3), (56.7, 18.9, 6.3)),
    ("FUTAGAWA (1) + UTO (2)", "FUT_UTO", 68.6,
     (7.1, 6.9, 6.7), (3.84, 1.28, 0.4266666667)),
    ("UTO (2) + UTO HANTO NORTH (3)", "UTO_UHN", 60.0,
     (7.2, 7.0, 6.8), (10.59, 3.53, 1.176666667)),
    ("FUTAGAWA (1) + UTO (2) + UTO HANTO NORTH (3)", "FUT_UTO_UHN", 64.6,
     (7.4, 7.2, 7.0), (3.84, 1.28, 0.4266666667)),
]

SOURCE_TEMPLATE = """      <characteristicFaultSource id="{sid}" name="{name}" tectonicRegion="Active Shallow Crust">
        <incrementalMFD binWidth="0.2" minMag="{min_mag}">
          <occurRates>{rates}</occurRates>
        </incrementalMFD>
        <rake>0.0</rake>
        <surface>
          <simpleFaultGeometry>
            <gml:LineString>
              <gml:posList>
{poslist}
              </gml:posList>
            </gml:LineString>
            <dip>{dip}</dip>
            <upperSeismoDepth>0.0</upperSeismoDepth>
            <lowerSeismoDepth>14.0</lowerSeismoDepth>
          </simpleFaultGeometry>
        </surface>
      </characteristicFaultSource>
"""


def read_traces():
    """Coordinates sheet -> {segment name: [(lon, lat), ...]} (NE -> SW)."""
    import openpyxl

    ws = openpyxl.load_workbook(WORKBOOK, data_only=True)["Coordinates"]
    traces, current = {}, None
    for row in ws.iter_rows(min_row=2, values_only=True):
        name, lat, lon = row[0], row[1], row[2]
        if isinstance(name, str) and name.strip():
            current = name.strip().upper()
            traces[current] = []
        if current and isinstance(lat, (int, float)) and isinstance(lon, (int, float)):
            traces[current].append((float(lon), float(lat)))
    return traces


def main():
    traces = read_traces()
    blocks = []
    for key, sid, dip, mags, rates in SOURCES:
        trace = traces[key.upper()]
        exp_pts = {"UTO": 3, "FUT_UTO": 5, "UTO_UHN": 5, "FUT_UTO_UHN": 7}
        assert len(trace) == exp_pts[sid], (sid, trace)
        mean_rate = sum(w * r for w, r in zip(RATE_WEIGHTS, rates)) * 1e-5
        # ascending magnitude bins; bin rate = magnitude weight x mean rate
        mags_asc = sorted(mags)
        assert all(abs(b - a - 0.2) < 1e-9
                   for a, b in zip(mags_asc, mags_asc[1:])), mags
        bin_rates = [MAG_WEIGHTS[list(mags).index(m)] * mean_rate
                     for m in mags_asc]
        poslist = "\n".join(f"                {lon:.5f} {lat:.5f}"
                            for lon, lat in trace)
        blocks.append(SOURCE_TEMPLATE.format(
            sid=sid, name=key.title(), min_mag=mags_asc[0],
            rates=" ".join(f"{r:.6e}" for r in bin_rates),
            poslist=poslist, dip=dip))
        print(f"{sid}: mean_rate={mean_rate:.4e} bins={bin_rates}")

    sm = ('<?xml version="1.0" encoding="utf-8"?>\n'
          '<nrml xmlns="http://openquake.org/xmlns/nrml/0.5"\n'
          '      xmlns:gml="http://www.opengis.net/gml">\n'
          '  <sourceModel name="IAEA Kumamoto base case (4 rupture sources,'
          ' mean-collapsed epistemic branches)">\n'
          '    <sourceGroup name="IAEA Kumamoto base case"'
          ' tectonicRegion="Active Shallow Crust">\n'
          + "".join(blocks) +
          '    </sourceGroup>\n  </sourceModel>\n</nrml>\n')
    (HERE / "source_model_basecase.xml").write_text(sm)

    smlt = ('<?xml version="1.0" encoding="UTF-8"?>\n'
            '<nrml xmlns:gml="http://www.opengis.net/gml"\n'
            '      xmlns="http://openquake.org/xmlns/nrml/0.4">\n'
            '  <logicTree logicTreeID="lt_sm_kumamoto_basecase">\n'
            '    <logicTreeBranchSet uncertaintyType="sourceModel"'
            ' branchSetID="bs_sm_kumamoto_basecase">\n'
            '      <logicTreeBranch branchID="b_kumamoto_basecase">\n'
            '        <uncertaintyModel>./source_model_basecase.xml'
            '</uncertaintyModel>\n'
            '        <uncertaintyWeight>1.0</uncertaintyWeight>\n'
            '      </logicTreeBranch>\n'
            '    </logicTreeBranchSet>\n  </logicTree>\n</nrml>\n')
    (HERE / "smlt_basecase.xml").write_text(smlt)
    print("wrote source_model_basecase.xml + smlt_basecase.xml")


if __name__ == "__main__":
    main()
