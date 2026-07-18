"""Logic-tree map-mode validation script.

One-fault validation mirroring ``../verify.py`` (curve mode). Runs the
three sibling map-mode configs::

    single_bilinear_map       - weight-1.0 tree on Petersen2011PrimaryFD (version=bilinear)
    single_elliptical_map     - weight-1.0 tree on Petersen2011PrimaryFD (version=elliptical)
    blend_50_50_map           - 50/50 LT over both models

through ``FdhaLogicTree`` (if outputs are missing) and checks:

1. Grid consistency: the three runs share the same N_sites, site lons/lats,
   D0 axis and return period.

2. Per-site arithmetic in RATE space (test 7, the map-mode analogue of
   the curve-mode ``verify.py`` post-processing):

       lambda_blend(i, j) == 0.5 * lambda_A(i, j) + 0.5 * lambda_B(i, j)

   on every site i and every D0 level j, within 1e-12. Here lambda_A and
   lambda_B come from the two weight-1.0 single-model runs and
   lambda_blend from the LT framework output - this is "compute each
   model individually, post-process the 50/50 average, compare to LT".

3. Inversion consistency in DISPLACEMENT space (test 8): manually invert
   rates_mean.h5 at the configured return period using
   ``openquake.fdha.calc.utils.interpolation.get_map_from_curves`` and
   compare to displacement_map_mean.csv within 1e-10 m.

4. End-to-end post-processing check: manually invert the analytical blend
       rates_analytical = 0.5 * rates_A + 0.5 * rates_B
   at the same return period, then compare to the LT framework's
   ``blend_50_50_map/out/aggregate/displacement_map_mean.csv`` within
   1e-10 m. This is the DISPLACEMENT-space version of step (2) - it
   proves the LT framework yields the same hazard map as the standalone
   per-model runs + external weighted averaging.

Exit code 0 on success, 1 on any failure.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
RUNS = {
    "single_bilinear_map": HERE / "single_bilinear_map",
    "single_elliptical_map": HERE / "single_elliptical_map",
    "blend_50_50_map": HERE / "blend_50_50_map",
}
OUT_SUBDIR = "out"


class CheckFailed(SystemExit):
    def __init__(self, msg: str) -> None:
        super().__init__(f"[FAIL] {msg}")


def _banner(title: str) -> None:
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def _ensure_runs() -> None:
    from openquake.fdha.logic_tree.driver import FdhaLogicTree

    for name, run_dir in RUNS.items():
        out = run_dir / OUT_SUBDIR / "aggregate" / "rates_mean.h5"
        if out.is_file():
            print(f"  [skip run] {name}: found existing {out.relative_to(HERE)}")
            continue
        ini = run_dir / "job.ini"
        print(f"  [run     ] {name}: {ini.relative_to(HERE)}")
        lt = FdhaLogicTree.from_ini(ini)
        lt.run(outdir=run_dir / OUT_SUBDIR)


def _read_rates_mean_h5(path: Path):
    import h5py

    with h5py.File(path, "r") as f:
        rates = np.asarray(f["rates_mean"][()], dtype=float)
        d0 = np.asarray(f["d0"][()], dtype=float)
        lons = np.asarray(f["site_lons"][()], dtype=float)
        lats = np.asarray(f["site_lats"][()], dtype=float)
    return rates, d0, lons, lats


def _read_displacement_map_csv(path: Path):
    rp = None
    rows: list[tuple[int, float, float, int, float]] = []
    with path.open() as f:
        for line in f:
            if line.startswith("#"):
                if "return_period" in line:
                    rp = float(line.split("=")[1].strip())
                continue
            break
        f.seek(0)
        reader = csv.reader(
            (r for r in f if not r.startswith("#"))
        )
        header = next(reader)
        # columns: site_id, lon, lat, is_trace, displ_*
        for row in reader:
            rows.append(
                (int(row[0]), float(row[1]), float(row[2]),
                 int(row[3]), float(row[4]))
            )
    return rp, header, rows


def _check_grid_consistency(data):
    ref_name, (ref_rates, ref_d0, ref_lons, ref_lats) = next(iter(data.items()))
    n_sites, n_d0 = ref_rates.shape
    print(f"  reference: {ref_name}  shape={ref_rates.shape}  "
          f"n_sites={n_sites} n_d0={n_d0}")
    for name, (r, d0, lons, lats) in data.items():
        if r.shape != ref_rates.shape:
            raise CheckFailed(
                f"shape mismatch: {name} has {r.shape}, expected {ref_rates.shape}"
            )
        if not np.array_equal(d0, ref_d0):
            raise CheckFailed(f"D0 mismatch between {ref_name} and {name}")
        if not np.allclose(lons, ref_lons) or not np.allclose(lats, ref_lats):
            raise CheckFailed(f"Site grid mismatch between {ref_name} and {name}")
    print("  OK: all three runs share grid, D0, site ordering.")


def _check_per_site_arithmetic(data, atol=1e-12):
    rA = data["single_bilinear_map"][0]
    rB = data["single_elliptical_map"][0]
    rC = data["blend_50_50_map"][0]
    expected = 0.5 * rA + 0.5 * rB
    diff = np.abs(rC - expected)
    max_abs = float(diff.max())
    print(f"  max|blend - 0.5*(A+B)| per-site = {max_abs:.3e}  (tol {atol})")
    if max_abs > atol:
        i_site, i_d0 = np.unravel_index(np.argmax(diff), diff.shape)
        raise CheckFailed(
            f"Per-site arithmetic FAILED at site={i_site}, D0_idx={i_d0}: "
            f"|blend-0.5(A+B)| = {max_abs:.3e} > {atol}."
        )
    print("  OK: per-site arithmetic holds to 1e-12.")


def _check_end_to_end_post_processing(data, atol=1e-10):
    """Step 4: compare LT blend displacement map against analytical inversion.

    Works the same way as the curve-mode ``verify.py`` except the
    post-processing is applied in rate space and then inverted, because
    averaging hazard MAPS directly would be physically wrong. Rate-space
    aggregation followed by inversion is the ONLY mathematically valid
    per-site path from (lambda_A, lambda_B) to the 50/50 blend hazard map.
    """
    from openquake.fdha.calc.utils.interpolation import get_map_from_curves

    rA, d0A, lonsA, latsA = data["single_bilinear_map"]
    rB, _, _, _ = data["single_elliptical_map"]

    blend_dir = RUNS["blend_50_50_map"] / OUT_SUBDIR / "aggregate"
    rp, header, rows = _read_displacement_map_csv(
        blend_dir / "displacement_map_mean.csv"
    )
    if rp is None:
        raise CheckFailed("blend_50_50_map displacement_map_mean.csv lacks return_period header")
    lt_lons = np.array([r[1] for r in rows])
    lt_lats = np.array([r[2] for r in rows])
    lt_displ = np.array([r[4] for r in rows])
    if not (np.allclose(lt_lons, lonsA) and np.allclose(lt_lats, latsA)):
        raise CheckFailed("LT blend displacement map site ordering differs from component runs")

    # Post-processing step: analytical blend rates -> invert at RP.
    rates_analytical = 0.5 * rA + 0.5 * rB
    pex = 1.0 / rp
    displ_analytical = get_map_from_curves(d0A, rates_analytical, pex)

    diff = np.abs(lt_displ - displ_analytical)
    max_abs = float(diff.max())
    print(f"  RP = {rp} yr  (pex = {pex:.3e})")
    print(
        f"  max|displ_LT_blend - displ_analytical_blend| = {max_abs:.3e} m  "
        f"(tol {atol})"
    )
    if max_abs > atol:
        i = int(np.argmax(diff))
        raise CheckFailed(
            f"End-to-end post-processing FAILED at site {i}: "
            f"LT={lt_displ[i]:.6f} m, analytical={displ_analytical[i]:.6f} m, "
            f"|diff|={max_abs:.3e} > {atol}."
        )
    nonzero = int((lt_displ > 0).sum())
    assert nonzero > 0, "All sites zero - test is vacuous."
    print(
        f"  OK: LT blend hazard map matches per-model runs + analytical "
        f"weighted average (non-zero sites: {nonzero}/{len(lt_displ)})."
    )


def _check_inversion_consistency(run_dir: Path, atol=1e-10):
    from openquake.fdha.calc.utils.interpolation import get_map_from_curves

    rates_path = run_dir / OUT_SUBDIR / "aggregate" / "rates_mean.h5"
    map_path = run_dir / OUT_SUBDIR / "aggregate" / "displacement_map_mean.csv"
    rates, d0, lons, lats = _read_rates_mean_h5(rates_path)
    rp, header, rows = _read_displacement_map_csv(map_path)
    if rp is None:
        raise CheckFailed(f"No '# return_period' header in {map_path}")

    # Compare site ordering.
    csv_lons = np.array([r[1] for r in rows])
    csv_lats = np.array([r[2] for r in rows])
    csv_displ = np.array([r[4] for r in rows])
    if not (np.allclose(csv_lons, lons) and np.allclose(csv_lats, lats)):
        raise CheckFailed("rates_mean.h5 and displacement_map_mean.csv site order mismatch")

    pex = 1.0 / rp
    manual = get_map_from_curves(d0, rates, pex)
    diff = np.abs(manual - csv_displ)
    max_abs = float(diff.max())
    print(f"  return_period  = {rp} yr  (pex={pex:.3e})")
    print(f"  max|manual - csv| = {max_abs:.3e}  (tol {atol})")
    if max_abs > atol:
        i = int(np.argmax(diff))
        raise CheckFailed(
            f"Inversion consistency FAILED at site {i}: "
            f"manual={manual[i]:.6f}, csv={csv_displ[i]:.6f}, "
            f"|diff|={max_abs:.3e} > {atol}."
        )
    print("  OK: displacement_map_mean.csv matches manual log-log inversion.")


def main() -> int:
    _banner("STEP 1  Ensure all three map-mode runs are present")
    _ensure_runs()

    _banner("STEP 2  Grid / D0 / site-ordering consistency")
    data = {
        name: _read_rates_mean_h5(run_dir / OUT_SUBDIR / "aggregate" / "rates_mean.h5")
        for name, run_dir in RUNS.items()
    }
    _check_grid_consistency(data)

    _banner("STEP 3  Per-site arithmetic (test 7): blend = 0.5 A + 0.5 B on every site")
    _check_per_site_arithmetic(data)

    _banner("STEP 4  Inversion consistency (test 8): rates_mean.h5 -> displacement_map_mean.csv")
    # Check all three scenarios' inversion; the log-log interpolation must
    # agree to 1e-10 regardless of which scenario we pick.
    for name, run_dir in RUNS.items():
        print(f"  scenario: {name}")
        _check_inversion_consistency(run_dir)

    _banner(
        "STEP 5  End-to-end post-processing check: run A alone, run B alone, "
        "blend them externally, compare to LT hazard map"
    )
    _check_end_to_end_post_processing(data)

    _banner("RESULT")
    print("  ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except CheckFailed as e:
        print(str(e), file=sys.stderr)
        sys.exit(1)
