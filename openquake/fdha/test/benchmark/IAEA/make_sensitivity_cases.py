#!/usr/bin/env python
"""Generate the IAEA exercise SENSITIVITY-CASE inputs from the author data.

The three author workbooks define, besides the base cases, sensitivity
cases with no published reference curves. They are implemented here as
demonstration jobs (validated by regression snapshots and structural
checks in ``test_sensitivity_cases.py``; see ``sensitivity.py`` for the
job list):

Kumamoto (workbook "Sensitive Case" sheets + Coordinates):
    sens1  Suizenji fault (5.4 km, T 5 km, dip 60 SW, Mw 5.8,
           23.3e-5/yr), site r = 0.6 km            -> distributed P11
    sens3  same Suizenji source, on-fault site      -> principal P11
    sens4  the four base-case rupture sources at their middle
           magnitude/rate branches, site r = 10 km  -> distributed P11

Le Teil (workbook "Sensitivity Case 1"/"Sensitive Case 4" +
``Cevenne_Fault_Segments`` shapefile):
    sens1  four fault-source combinations (w 0.25 each; LRF_2 itself
           with two rupture-length branches w 0.4/0.6), site r = 0.6 km
                                                     -> distributed T13
    sens4  three separate faults (LRF_2, MRF_3, PCF_2) summed,
           site r = 0.6 km                           -> distributed T13

Norcia (workbook "Sensitivity Case 2"/"Sensitivity Case 3" + Coordinates):
    sens2  MVFS at the r = 2.4 km site               -> distributed Y03
    sens3  MVFS + NFS (29 km, T 12 km, dip 50 SW, GR-like MFD
           Mw 5.5-6.8), site 7.6 km HW / 3 km FW     -> distributed Y03

Regenerating (requires the author workbook + shapefile)::

    PYTHONPATH=. python openquake/fdha/test/benchmark/IAEA/make_sensitivity_cases.py
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
INPUT_ROOT = Path("/home/ychen/Documents/PFDHA/IAEA_benchmarking/"
                  "Input for Modeling-20260319T133033Z-3-001/"
                  "Input for Modeling")

SM_HEAD = ('<?xml version="1.0" encoding="utf-8"?>\n'
           '<nrml xmlns="http://openquake.org/xmlns/nrml/0.5"\n'
           '      xmlns:gml="http://www.opengis.net/gml">\n'
           '  <sourceModel name="{name}">\n'
           '    <sourceGroup name="{name}"'
           ' tectonicRegion="Active Shallow Crust">\n')
SM_TAIL = '    </sourceGroup>\n  </sourceModel>\n</nrml>\n'

SOURCE_TEMPLATE = """      <characteristicFaultSource id="{sid}" name="{sname}" tectonicRegion="Active Shallow Crust">
        <incrementalMFD binWidth="{bin_width}" minMag="{min_mag}">
          <occurRates>{rates}</occurRates>
        </incrementalMFD>
        <rake>{rake}</rake>
        <surface>
          <simpleFaultGeometry>
            <gml:LineString>
              <gml:posList>
{poslist}
              </gml:posList>
            </gml:LineString>
            <dip>{dip}</dip>
            <upperSeismoDepth>0.0</upperSeismoDepth>
            <lowerSeismoDepth>{lsd}</lowerSeismoDepth>
          </simpleFaultGeometry>
        </surface>
      </characteristicFaultSource>
"""


def source_block(sid, sname, trace, dip, lsd, rake, min_mag, rates,
                 bin_width=0.1):
    poslist = "\n".join(f"                {lon:.5f} {lat:.5f}"
                        for lon, lat in trace)
    return SOURCE_TEMPLATE.format(
        sid=sid, sname=sname, bin_width=bin_width, min_mag=min_mag,
        rates=" ".join(f"{r:.6e}" for r in rates), rake=rake,
        poslist=poslist, dip=dip, lsd=lsd)


def write_model(path, name, blocks):
    path.write_text(SM_HEAD.format(name=name) + "".join(blocks) + SM_TAIL)
    print(f"wrote {path.name}")


def write_smlt(path, lt_id, branches):
    """branches: [(source model filename, weight, comment)]"""
    lines = ['<?xml version="1.0" encoding="UTF-8"?>',
             '<nrml xmlns:gml="http://www.opengis.net/gml"',
             '      xmlns="http://openquake.org/xmlns/nrml/0.4">',
             f'  <logicTree logicTreeID="{lt_id}">',
             f'    <logicTreeBranchSet uncertaintyType="sourceModel"'
             f' branchSetID="bs_{lt_id}">']
    for i, (fname, w, comment) in enumerate(branches):
        lines += [f'      <logicTreeBranch branchID="b_{lt_id}_{i:02d}">',
                  f'        <!-- {comment} -->',
                  f'        <uncertaintyModel>./{fname}</uncertaintyModel>',
                  f'        <uncertaintyWeight>{w:g}</uncertaintyWeight>',
                  '      </logicTreeBranch>']
    lines += ['    </logicTreeBranchSet>', '  </logicTree>', '</nrml>', '']
    path.write_text("\n".join(lines))
    print(f"wrote {path.name}")


# --------------------------------------------------------------- Kumamoto

def kumamoto():
    d = HERE / "kumamoto"
    sys.path.insert(0, str(d))
    from make_basecase import read_traces, SOURCES  # noqa: E402

    traces = read_traces()

    # Suizenji sits in the side table (columns D-F) of the Coordinates sheet
    import openpyxl
    from make_basecase import WORKBOOK
    ws = openpyxl.load_workbook(WORKBOOK, data_only=True)["Coordinates"]
    suizenji, active = [], False
    for row in ws.iter_rows(min_row=2, values_only=True):
        name, lat, lon = (row[4], row[5], row[6]) if len(row) > 6 else (None,) * 3
        if isinstance(name, str) and name.strip():
            active = name.strip().upper() == "SUIZENJI"
        if active and isinstance(lat, (int, float)) and isinstance(lon, (int, float)):
            suizenji.append((float(lon), float(lat)))
        elif active and suizenji and lat is None:
            break
    assert len(suizenji) == 3, suizenji

    # sens1/sens3 source: Suizenji (5.4 km, T 5, dip 60 SW, M5.8, 23.3e-5)
    write_model(d / "source_model_sens1_suizenji.xml",
                "IAEA Kumamoto sensitivity 1/3: Suizenji fault",
                [source_block("SUIZENJI", "Suizenji Fault", suizenji,
                              dip=60, lsd=5.0, rake=0.0, min_mag=5.8,
                              rates=[23.3e-5])])
    write_smlt(d / "smlt_sens1.xml", "lt_sm_kumamoto_sens1",
               [("source_model_sens1_suizenji.xml", 1.0, "Suizenji fault")])

    # sens4: the four rupture sources at their middle (w 0.6) branches
    middle = {"UTO": (6.5, 18.9e-5), "FUT_UTO": (6.9, 1.28e-5),
              "UTO_UHN": (7.0, 3.53e-5), "FUT_UTO_UHN": (7.2, 1.28e-5)}
    blocks = []
    for key, sid, dip, _, _ in SOURCES:
        mag, rate = middle[sid]
        blocks.append(source_block(f"{sid}_S4", key.title(),
                                   traces[key.upper()], dip=dip, lsd=14.0,
                                   rake=0.0, min_mag=mag, rates=[rate]))
    write_model(d / "source_model_sens4.xml",
                "IAEA Kumamoto sensitivity 4: four rupture sources"
                " (middle branches)", blocks)
    write_smlt(d / "smlt_sens4.xml", "lt_sm_kumamoto_sens4",
               [("source_model_sens4.xml", 1.0, "4 sources, middle branches")])


# ---------------------------------------------------------------- Le Teil

def leteil_traces():
    import shapefile

    sf = shapefile.Reader(str(
        INPUT_ROOT / "Le Teil Case - Reverse" / "Shapefile"
        / "Cevenne_Fault_Segments"))
    segs = {sr.record["Name"]: [tuple(p) for p in sr.shape.points]
            for sr in sf.shapeRecords()}

    def sw_ne(pts):  # order along strike SW -> NE (by latitude)
        return pts if pts[0][1] <= pts[-1][1] else pts[::-1]

    def joined(*names):
        """Concatenate segments SW->NE, merging near-coincident junctions."""
        pts = []
        for name in names:
            for p in sw_ne(segs[name]):
                if pts and (abs(p[0] - pts[-1][0]) + abs(p[1] - pts[-1][1])
                            < 0.002):  # ~150 m junction: midpoint-merge
                    pts[-1] = ((p[0] + pts[-1][0]) / 2,
                               (p[1] + pts[-1][1]) / 2)
                else:
                    pts.append(p)
        return pts

    return {
        "LRF2_5.8": segs["LRF2_Extended5.8km"],
        "LRF2_4.6": segs["LRF2"],
        "LRF2_LRF1": joined("LRF2", "LRF1"),
        "LRF3_LRF2": joined("LRF3", "LRF2"),
        "LRF3_LRF2_LRF1": joined("LRF3", "LRF2", "LRF1"),
        "MRF3": sw_ne(segs["MRF3"]),
        "PCF2": sw_ne(segs["PCF2Reduced3.2km"]),
    }


def leteil():
    d = HERE / "le_teil"
    tr = leteil_traces()

    # Sensitivity Case 1: four fault-source combinations (w 0.25; LRF_2
    # itself has rupture-length branches 5.8 km w 0.4 / 4.6 km w 0.6)
    sens1 = [
        # (tag, trace key, thickness, mag, rate, weight, comment)
        ("L5p8", "LRF2_5.8", 4.0, 5.52, 4.6e-5, 0.25 * 0.4,
         "LRF_2, rupture length 5.8 km"),
        ("L4p6", "LRF2_4.6", 2.0, 4.9, 1.55e-4, 0.25 * 0.6,
         "LRF_2, rupture length 4.6 km"),
        ("LRF21", "LRF2_LRF1", 6.85, 6.1, 2.38e-5, 0.25,
         "LRF_2 + LRF_1 (13 km)"),
        ("LRF23", "LRF3_LRF2", 4.0, 5.66, 3.92e-5, 0.25,
         "LRF_2 + LRF_3 (8 km)"),
        ("LRF213", "LRF3_LRF2_LRF1", 10.0, 6.35, 1.81e-5, 0.25,
         "LRF_2 + LRF_1 + LRF_3 (16 km)"),
    ]
    branches = []
    for tag, key, lsd, mag, rate, w, comment in sens1:
        fname = f"source_model_sens1_{tag}.xml"
        write_model(d / fname, f"IAEA Le Teil sensitivity 1: {comment}",
                    [source_block(f"LRF_{tag}", comment, tr[key], dip=45,
                                  lsd=lsd, rake=90.0, min_mag=mag,
                                  rates=[rate])])
        branches.append((fname, w, comment))
    assert abs(sum(w for _, w, _ in branches) - 1.0) < 1e-9
    write_smlt(d / "smlt_sens1.xml", "lt_sm_le_teil_sens1", branches)

    # Sensitive Case 4: three separate faults summed
    blocks = [
        source_block("LRF2_S4", "La Rouviere LRF_2", tr["LRF2_5.8"],
                     dip=45, lsd=4.0, rake=90.0, min_mag=5.52,
                     rates=[4.6e-5]),
        source_block("MRF3_S4", "MRF_3", tr["MRF3"],
                     dip=80, lsd=6.0, rake=90.0, min_mag=5.6,
                     rates=[4.22e-5]),
        source_block("PCF2_S4", "PCF_2", tr["PCF2"],
                     dip=70, lsd=2.0, rake=90.0, min_mag=4.9,
                     rates=[8.64e-5]),
    ]
    write_model(d / "source_model_sens4.xml",
                "IAEA Le Teil sensitivity 4: LRF_2 + MRF_3 + PCF_2", blocks)
    write_smlt(d / "smlt_sens4.xml", "lt_sm_le_teil_sens4",
               [("source_model_sens4.xml", 1.0, "3 faults summed")])


# ----------------------------------------------------------------- Norcia

# NFS trace and GR-like MFD from the workbook (Coordinates /
# "Sensitivity Case 3" sheets)
NFS_TRACE = [(13.0302, 42.8982), (13.0626, 42.8701), (13.0943, 42.8566),
             (13.1204, 42.8253), (13.1206, 42.8003), (13.1733, 42.7064),
             (13.1602, 42.6739)]
NFS_RATES = [1.1582e-3, 9.1998e-4, 7.3077e-4, 5.8047e-4, 4.6108e-4,
             3.6625e-4, 2.9092e-4, 2.3109e-4, 1.8356e-4, 1.4581e-4,
             1.1582e-4, 9.1998e-5, 7.3077e-5, 5.8047e-5]  # Mw 5.5..6.8
MVFS_RATES = [4.2864e-5, 5.6589e-5, 6.6852e-5, 7.0671e-5, 6.6852e-5,
              5.6589e-5, 4.2864e-5]  # Mw 6.4..7.0


def norcia():
    d = HERE / "norcia"
    # sens3: MVFS + NFS in one model (MVFS block reused from the base model)
    mvfs_xml = (d / "source_model_mvfs.xml").read_text()
    start = mvfs_xml.index("      <characteristicFaultSource")
    end = mvfs_xml.index("</characteristicFaultSource>") + \
        len("</characteristicFaultSource>")
    mvfs_block = mvfs_xml[start:end] + "\n"
    nfs_block = source_block("NFS", "Norcia Fault System", NFS_TRACE,
                             dip=50, lsd=12.0, rake=-90.0, min_mag=5.5,
                             rates=NFS_RATES)
    write_model(d / "source_model_sens3.xml",
                "IAEA Norcia sensitivity 3: MVFS + NFS",
                [mvfs_block, nfs_block])
    write_smlt(d / "smlt_sens3.xml", "lt_sm_norcia_sens3",
               [("source_model_sens3.xml", 1.0, "MVFS + NFS")])


if __name__ == "__main__":
    kumamoto()
    leteil()
    norcia()
