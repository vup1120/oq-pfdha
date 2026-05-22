from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np


def write_branch_rates_csv(
    path: Path,
    d0: list[float],
    rates: Any,
    site_lons: list[float],
    site_lats: list[float],
) -> None:
    """
    rates: (n_sites, n_d0)
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    n_sites = len(site_lons)
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        if n_sites == 1:
            w.writerow(["D0", "annual_rate"])
            for j, dd in enumerate(d0):
                w.writerow([dd, float(rates[0][j])])
        else:
            w.writerow(["site_id", "lon", "lat", "D0", "annual_rate"])
            for i in range(n_sites):
                for j, dd in enumerate(d0):
                    w.writerow([i, float(site_lons[i]), float(site_lats[i]), dd, float(rates[i][j])])


def write_aggregate_csv(
    path: Path,
    d0: list[float],
    mean_rates: Any,
    fractiles: dict[float, Any],
    site_lons: list[float],
    site_lats: list[float],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    n_sites = len(site_lons)
    qs = [0.05, 0.16, 0.5, 0.84, 0.95]
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        if n_sites == 1:
            w.writerow(["D0", "mean", "p05", "p16", "p50", "p84", "p95"])
            for j, dd in enumerate(d0):
                row = [dd, float(mean_rates[0][j])]
                row.extend(float(fractiles[q][0][j]) for q in qs)
                w.writerow(row)
        else:
            w.writerow(["site_id", "lon", "lat", "D0", "mean", "p05", "p16", "p50", "p84", "p95"])
            for i in range(n_sites):
                for j, dd in enumerate(d0):
                    row = [i, float(site_lons[i]), float(site_lats[i]), dd, float(mean_rates[i][j])]
                    row.extend(float(fractiles[q][i][j]) for q in qs)
                    w.writerow(row)


def write_manifest(path: Path, manifest: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)


def write_validator_report(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def write_hazard_displacement_map_csv(
    path: Path,
    *,
    site_lons: Sequence[float],
    site_lats: Sequence[float],
    displ_mean: Sequence[float],
    displ_fractiles: dict[float, Sequence[float]] | None = None,
    target_return_period: float | None = None,
    site_is_trace: Sequence[bool] | None = None,
) -> None:
    """Write a per-site hazard-map CSV (one row per site).

    Columns: ``site_id, lon, lat, is_trace, displ_mean`` and optionally
    ``displ_p05, displ_p16, displ_p50, displ_p84, displ_p95`` when fractile
    maps are supplied. ``target_return_period`` is recorded as a
    ``# return_period = T`` header comment.
    """
    site_lons_a = [float(v) for v in site_lons]
    site_lats_a = [float(v) for v in site_lats]
    displ_mean_a = [float(v) for v in displ_mean]
    if len(site_lons_a) != len(site_lats_a) or len(site_lons_a) != len(displ_mean_a):
        raise ValueError("lons/lats/displ_mean length mismatch")
    n = len(site_lons_a)
    if site_is_trace is None:
        site_is_trace_a = [False] * n
    else:
        site_is_trace_a = [bool(v) for v in site_is_trace]
    qs = [0.05, 0.16, 0.5, 0.84, 0.95]
    frac_cols: list[tuple[str, list[float]]] = []
    if displ_fractiles:
        for q in qs:
            if q in displ_fractiles:
                frac_cols.append((f"displ_p{int(round(q*100)):02d}",
                                  [float(v) for v in displ_fractiles[q]]))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        if target_return_period is not None:
            f.write(f"# return_period = {float(target_return_period)}\n")
        w = csv.writer(f)
        header = ["site_id", "lon", "lat", "is_trace", "displ_mean"]
        header.extend(name for name, _ in frac_cols)
        w.writerow(header)
        for i in range(n):
            row = [i, site_lons_a[i], site_lats_a[i], int(site_is_trace_a[i]),
                   displ_mean_a[i]]
            row.extend(col[i] for _, col in frac_cols)
            w.writerow(row)


def write_branch_rate_cube_npz(
    path: Path,
    *,
    branch_fingerprints: Sequence[str],
    branch_weights: Sequence[float],
    d0: Sequence[float],
    site_lons: Sequence[float],
    site_lats: Sequence[float],
    rates: Any,
) -> None:
    """Store per-branch per-site annual rates as a compressed ``.npz`` file.

    Used by map-mode logic-tree execution. Single-site hazard-curve runs
    keep emitting per-branch CSVs via :func:`write_branch_rates_csv`.

    Shape contract:
        ``rates`` has shape ``(n_branches, n_sites, n_d0)``.
        ``branch_fingerprints`` / ``branch_weights`` length ``n_branches``.
        ``site_lons`` / ``site_lats`` length ``n_sites``.
        ``d0`` length ``n_d0``.

    The aggregate (mean + fractiles) and the derived displacement map are
    written as CSV elsewhere; only the branch rate cube uses this format.
    """
    rates_arr = np.asarray(rates, dtype=float)
    expected_shape = (len(branch_fingerprints), len(site_lons), len(d0))
    if rates_arr.shape != expected_shape:
        raise ValueError(
            f"rates shape {rates_arr.shape} does not match "
            f"(n_branches, n_sites, n_d0)={expected_shape}"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        branch_fingerprints=np.asarray(branch_fingerprints, dtype="U"),
        branch_weights=np.asarray(branch_weights, dtype=float),
        d0=np.asarray(d0, dtype=float),
        site_lons=np.asarray(site_lons, dtype=float),
        site_lats=np.asarray(site_lats, dtype=float),
        rates=rates_arr,
    )


# -----------------------------------------------------------------------
# v4 IO: one HDF5 per branch + aggregate HDF5 + per-fractile map CSVs
# -----------------------------------------------------------------------

# Fixed fractile order (contract with tests & downstream readers)
FRACTILE_QS: tuple[float, ...] = (0.05, 0.16, 0.5, 0.84, 0.95)
FRACTILE_LABELS: tuple[str, ...] = ("p05", "p16", "p50", "p84", "p95")


def write_branch_rates_h5(
    path: Path,
    *,
    rates: Any,
    d0: Sequence[float],
    site_lons: Sequence[float],
    site_lats: Sequence[float],
    weight: float,
    fingerprint: str,
    attrs: dict[str, Any] | None = None,
) -> None:
    """Write a single end-branch per-site rate grid to HDF5.

    Datasets:
        rates     shape (n_sites, n_d0), dtype float64
        d0        shape (n_d0,)
        site_lons shape (n_sites,)
        site_lats shape (n_sites,)

    Attributes (on the root group): ``weight``, ``fingerprint`` and any
    extras supplied via ``attrs`` (scalars or short strings only).
    """
    try:
        import h5py  # type: ignore
    except Exception as exc:  # pragma: no cover - optional dep
        raise NotImplementedError(
            "v4 map-mode IO requires h5py; install h5py to proceed."
        ) from exc

    rates_arr = np.asarray(rates, dtype=float)
    expected = (len(site_lons), len(d0))
    if rates_arr.shape != expected:
        raise ValueError(
            f"branch rates shape {rates_arr.shape} != (n_sites, n_d0)={expected}"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        f.create_dataset("rates", data=rates_arr, compression="gzip")
        f.create_dataset("d0", data=np.asarray(d0, dtype=float))
        f.create_dataset("site_lons", data=np.asarray(site_lons, dtype=float))
        f.create_dataset("site_lats", data=np.asarray(site_lats, dtype=float))
        f.attrs["weight"] = float(weight)
        f.attrs["fingerprint"] = str(fingerprint)
        if attrs:
            for k, v in attrs.items():
                try:
                    f.attrs[k] = v
                except TypeError:
                    f.attrs[k] = str(v)


def write_rates_mean_h5(
    path: Path,
    *,
    rates_mean: Any,
    d0: Sequence[float],
    site_lons: Sequence[float],
    site_lats: Sequence[float],
) -> None:
    """Write the LT-weighted per-site mean annual-rate grid as HDF5."""
    try:
        import h5py  # type: ignore
    except Exception as exc:  # pragma: no cover
        raise NotImplementedError(
            "v4 map-mode IO requires h5py; install h5py to proceed."
        ) from exc

    arr = np.asarray(rates_mean, dtype=float)
    expected = (len(site_lons), len(d0))
    if arr.shape != expected:
        raise ValueError(
            f"rates_mean shape {arr.shape} != (n_sites, n_d0)={expected}"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        f.create_dataset("rates_mean", data=arr, compression="gzip")
        f.create_dataset("d0", data=np.asarray(d0, dtype=float))
        f.create_dataset("site_lons", data=np.asarray(site_lons, dtype=float))
        f.create_dataset("site_lats", data=np.asarray(site_lats, dtype=float))


def write_rates_fractiles_h5(
    path: Path,
    *,
    rates_fractiles: Any,
    d0: Sequence[float],
    site_lons: Sequence[float],
    site_lats: Sequence[float],
) -> None:
    """Write per-site fractile annual-rate grids as HDF5.

    ``rates_fractiles`` shape ``(5, n_sites, n_d0)`` with quantile order
    fixed by :data:`FRACTILE_QS` = (0.05, 0.16, 0.50, 0.84, 0.95). The
    order is also recorded as an attribute ``quantiles`` and as a string
    dataset ``quantile_labels`` for self-description.
    """
    try:
        import h5py  # type: ignore
    except Exception as exc:  # pragma: no cover
        raise NotImplementedError(
            "v4 map-mode IO requires h5py; install h5py to proceed."
        ) from exc

    arr = np.asarray(rates_fractiles, dtype=float)
    expected = (len(FRACTILE_QS), len(site_lons), len(d0))
    if arr.shape != expected:
        raise ValueError(
            f"rates_fractiles shape {arr.shape} != "
            f"(n_quantiles, n_sites, n_d0)={expected}"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        f.create_dataset("rates_fractiles", data=arr, compression="gzip")
        f.create_dataset("quantiles", data=np.asarray(FRACTILE_QS, dtype=float))
        f.create_dataset(
            "quantile_labels",
            data=np.asarray(FRACTILE_LABELS, dtype="S"),
        )
        f.create_dataset("d0", data=np.asarray(d0, dtype=float))
        f.create_dataset("site_lons", data=np.asarray(site_lons, dtype=float))
        f.create_dataset("site_lats", data=np.asarray(site_lats, dtype=float))
        f.attrs["quantiles"] = np.asarray(FRACTILE_QS, dtype=float)


def write_displacement_map_csv(
    path: Path,
    *,
    site_lons: Sequence[float],
    site_lats: Sequence[float],
    displ: Sequence[float],
    target_return_period: float,
    label: str = "displ",
    site_is_trace: Sequence[bool] | None = None,
) -> None:
    """Write one displacement-map CSV for a single aggregation level.

    Columns: ``site_id, lon, lat, is_trace, <label>``. A header comment
    ``# return_period = T`` records the inversion target. Used by v4
    map-mode to emit one CSV per aggregation level
    (``displacement_map_mean.csv``, ``displacement_map_p05.csv``, ...).
    """
    lons = [float(v) for v in site_lons]
    lats = [float(v) for v in site_lats]
    vals = [float(v) for v in displ]
    if not (len(lons) == len(lats) == len(vals)):
        raise ValueError("lons/lats/displ length mismatch")
    n = len(lons)
    if site_is_trace is None:
        is_trace = [False] * n
    else:
        is_trace = [bool(v) for v in site_is_trace]
        if len(is_trace) != n:
            raise ValueError("site_is_trace length mismatch")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        f.write(f"# return_period = {float(target_return_period)}\n")
        w = csv.writer(f)
        w.writerow(["site_id", "lon", "lat", "is_trace", label])
        for i in range(n):
            w.writerow([i, lons[i], lats[i], int(is_trace[i]), vals[i]])


def write_branch_rate_cube_hdf5(
    path: Path,
    *,
    branch_fingerprints: Sequence[str],
    branch_weights: Sequence[float],
    d0: Sequence[float],
    site_lons: Sequence[float],
    site_lats: Sequence[float],
    rates: Any,
) -> None:
    """HDF5 variant of :func:`write_branch_rate_cube_npz`.

    Uses :mod:`h5py` when available; raises :class:`NotImplementedError` with a
    clear message otherwise. Shape contract matches the NPZ writer.
    """
    try:
        import h5py  # type: ignore
    except Exception as exc:  # pragma: no cover - optional dep
        raise NotImplementedError(
            "HDF5 branch rate cube requires h5py; install h5py or use "
            "write_branch_rate_cube_npz instead."
        ) from exc

    rates_arr = np.asarray(rates, dtype=float)
    expected_shape = (len(branch_fingerprints), len(site_lons), len(d0))
    if rates_arr.shape != expected_shape:
        raise ValueError(
            f"rates shape {rates_arr.shape} does not match "
            f"(n_branches, n_sites, n_d0)={expected_shape}"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        f.create_dataset("rates", data=rates_arr, compression="gzip")
        f.create_dataset(
            "branch_fingerprints",
            data=np.asarray(branch_fingerprints, dtype="S"),
        )
        f.create_dataset("branch_weights", data=np.asarray(branch_weights, dtype=float))
        f.create_dataset("d0", data=np.asarray(d0, dtype=float))
        f.create_dataset("site_lons", data=np.asarray(site_lons, dtype=float))
        f.create_dataset("site_lats", data=np.asarray(site_lats, dtype=float))

