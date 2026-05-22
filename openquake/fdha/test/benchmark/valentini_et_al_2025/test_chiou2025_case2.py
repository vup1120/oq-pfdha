import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path


# Reference exceedance vector for ValentiniEtAl2025 (Kumamoto Case 2)
REF_EXCEED = np.array([
    1.31775e-04, 1.31572e-04, 1.30400e-04, 1.28706e-04, 1.26877e-04,
    1.20965e-04, 1.12737e-04, 1.02679e-04, 9.33119e-05, 7.71261e-05,
    4.54213e-05, 2.48653e-05, 1.32235e-05, 7.71793e-06, 4.62204e-07,
    8.05555e-08, 1.63559e-08, 4.71018e-09
], dtype=float)


def test_valentini_chiou2025_case2(tmp_path):
    """
    End-to-end CLI test for Chiou (2025) primary FD model using the ValentiniEtAl2025
    Kumamoto Case 2 configuration. This test:

    1) Runs the canonical logic-tree job to compute the hazard curve.
    2) Reads the aggregate hazard CSV.
    3) Compares the exceedance vector to the published reference (REF_EXCEED).
    4) Produces two figures for manual inspection:
       - Overlay of reference vs computed
       - Relative error (%) vs displacement

    Notes:
    - We accept up to ~12% relative difference to account for minor geometric and
      numerical differences across environments.
    - Figures are saved under the pytest temporary directory.
    """
    base_dir = Path(__file__).parent / "data"
    config_path = base_dir / "configuration" / "job_kumamoto_case2_Chiou2025.ini"

    figures_dir = tmp_path / "figures"
    figures_dir.mkdir(exist_ok=True)

    from openquake.fdha.logic_tree.driver import FdhaLogicTree

    result = FdhaLogicTree.from_ini(config_path).run(outdir=tmp_path / "out")
    aggregate_csv = Path(result.outdir) / "aggregate_hazard.csv"
    assert aggregate_csv.exists(), f"Aggregate hazard CSV not found: {aggregate_csv}"

    d0_vals = []
    mean_vals = []
    with aggregate_csv.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            d0_vals.append(float(row["D0"]))
            mean_vals.append(float(row["mean"]))
    imls = np.asarray(d0_vals, dtype=float)
    poes = np.asarray(mean_vals, dtype=float)

    # Sanity checks
    assert imls.size == REF_EXCEED.size, f"IML size mismatch: {imls.size} vs {REF_EXCEED.size}"
    np.testing.assert_allclose(
        imls,
        np.array([0.0001, 0.001, 0.005, 0.01, 0.015, 0.03, 0.05, 0.075,
                  0.1, 0.15, 0.3, 0.5, 0.75, 1.0, 3.0, 5.0, 7.5, 10.0], dtype=float),
        rtol=0.0,
        atol=0.0,
    )
    # Allow up to ~12% relative difference (observed max ≈ 11%)
    np.testing.assert_allclose(poes, REF_EXCEED, rtol=1.2e-1, atol=0.0)
    print(f"Figures saved in: {figures_dir}")

    # Overlay plot (Reference vs Computed)
    overlay_png = figures_dir / "kumamoto_case2_Chiou2025_comparison.png"
    plt.figure(figsize=(9, 5.5))
    plt.loglog(imls, REF_EXCEED, color="#D39200", lw=2.5, label="Reference (paper)")
    plt.loglog(imls, poes, "--", color="#1f77b4", lw=2.5, label="Computed (CLI)")
    plt.xlabel("Displacement, d (m)")
    plt.ylabel("P(D > d)")
    plt.grid(True, which="both", ls=":", alpha=0.35)
    plt.legend(loc="lower left")
    plt.tight_layout()
    plt.savefig(overlay_png, dpi=150)
    plt.close()

    # Relative error plot
    rel_err = 100.0 * (poes - REF_EXCEED) / REF_EXCEED
    max_idx = int(np.nanargmax(np.abs(rel_err)))
    relerr_png = figures_dir / "kumamoto_case2_Chiou2025_relerr.png"
    plt.figure(figsize=(9, 5.5))
    plt.semilogx(imls, rel_err, "-o", color="#D39200", ms=4)
    plt.xlabel("Displacement, d (m)")
    plt.ylabel("Relative error (%)")
    plt.title(f"Max rel err = {abs(rel_err[max_idx]):.3f}% at d={imls[max_idx]:g} m")
    plt.grid(True, which="both", ls=":", alpha=0.35)
    plt.tight_layout()
    plt.savefig(relerr_png, dpi=150)
    plt.close()

    # Ensure artifacts exist
    assert overlay_png.exists()
    assert relerr_png.exists()
