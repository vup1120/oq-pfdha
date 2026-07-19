"""CSV/HDF5 writers and the run manifest for logic-tree results."""
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
    rates_principal: Any = None,
    rates_distributed: Any = None,
) -> None:
    """
    rates: (n_sites, n_d0)

    When both ``rates_principal`` and ``rates_distributed`` are given (same
    shape as ``rates``), two extra columns ``annual_rate_principal`` and
    ``annual_rate_distributed`` are appended after ``annual_rate`` so users
    can inspect the on-fault (principal) and off-fault (distributed)
    contributions separately; their sum equals ``annual_rate``.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    n_sites = len(site_lons)
    has_split = rates_principal is not None and rates_distributed is not None
    split_cols = ["annual_rate_principal", "annual_rate_distributed"] if has_split else []

    def _split_vals(i: int, j: int) -> list[float]:
        if not has_split:
            return []
        return [float(rates_principal[i][j]), float(rates_distributed[i][j])]

    with path.open("w", newline="") as f:
        w = csv.writer(f)
        if n_sites == 1:
            w.writerow(["D0", "annual_rate"] + split_cols)
            for j, dd in enumerate(d0):
                w.writerow([dd, float(rates[0][j])] + _split_vals(0, j))
        else:
            w.writerow(["site_id", "lon", "lat", "D0", "annual_rate"] + split_cols)
            for i in range(n_sites):
                for j, dd in enumerate(d0):
                    w.writerow(
                        [i, float(site_lons[i]), float(site_lats[i]), dd,
                         float(rates[i][j])] + _split_vals(i, j))


def write_aggregate_csv(
    path: Path,
    d0: list[float],
    mean_rates: Any,
    fractiles: dict[float, Any],
    site_lons: list[float],
    site_lats: list[float],
    qs: Sequence[float] | None = None,
    include_mean: bool = True,
    mean_principal: Any = None,
    mean_distributed: Any = None,
) -> None:
    """Write the aggregate hazard CSV (weighted mean and/or fractiles).

    ``qs`` selects the fractile columns (default :data:`FRACTILE_QS`); pass an
    empty sequence for no fractiles. ``include_mean`` toggles the ``mean``
    column. Fractile columns are labelled OpenQuake-style via
    :func:`quantile_label` (e.g. ``quantile-0.05``).

    When ``include_mean`` is true and both ``mean_principal`` and
    ``mean_distributed`` are given, ``mean_principal`` and
    ``mean_distributed`` columns follow ``mean``. The mean is linear in the
    branch rates, so the two component columns sum exactly to ``mean``.
    Fractiles are reported for the total only.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    n_sites = len(site_lons)
    qs = list(FRACTILE_QS) if qs is None else list(qs)
    frac_labels = [quantile_label(q) for q in qs]
    has_split = (
        include_mean and mean_principal is not None
        and mean_distributed is not None
    )
    mean_cols = ["mean"] if include_mean else []
    if has_split:
        mean_cols += ["mean_principal", "mean_distributed"]

    def _mean_vals(i: int, j: int) -> list[float]:
        vals: list[float] = []
        if include_mean:
            vals.append(float(mean_rates[i][j]))
        if has_split:
            vals.append(float(mean_principal[i][j]))
            vals.append(float(mean_distributed[i][j]))
        return vals

    with path.open("w", newline="") as f:
        w = csv.writer(f)
        if n_sites == 1:
            w.writerow(["D0"] + mean_cols + frac_labels)
            for j, dd in enumerate(d0):
                row: list[Any] = [dd]
                row.extend(_mean_vals(0, j))
                row.extend(float(fractiles[q][0][j]) for q in qs)
                w.writerow(row)
        else:
            w.writerow(["site_id", "lon", "lat", "D0"] + mean_cols + frac_labels)
            for i in range(n_sites):
                for j, dd in enumerate(d0):
                    row = [i, float(site_lons[i]), float(site_lats[i]), dd]
                    row.extend(_mean_vals(i, j))
                    row.extend(float(fractiles[q][i][j]) for q in qs)
                    w.writerow(row)


def write_manifest(path: Path, manifest: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)


def write_validator_report(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


# -----------------------------------------------------------------------
# v4 IO: one HDF5 per branch + aggregate HDF5 + per-fractile map CSVs
# -----------------------------------------------------------------------

def quantile_label(q: float) -> str:
    """Return an OpenQuake-style fractile label, e.g. ``0.05`` -> ``quantile-0.05``.

    ``%g`` formatting drops trailing zeros (``0.5`` -> ``quantile-0.5``) and
    is unambiguous for arbitrary quantiles (``0.025`` -> ``quantile-0.025``),
    matching the ``quantile-<q>`` naming used by the OpenQuake Engine.
    """
    return "quantile-%g" % float(q)


# Default fractile set when ``[output].quantiles`` is not specified in the job.
# (The effective set is configurable; see openquake.fdha.calc.config_loader.)
FRACTILE_QS: tuple[float, ...] = (0.05, 0.16, 0.5, 0.84, 0.95)
FRACTILE_LABELS: tuple[str, ...] = tuple(quantile_label(q) for q in FRACTILE_QS)


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
    rates_principal: Any = None,
    rates_distributed: Any = None,
) -> None:
    """Write a single end-branch per-site rate grid to HDF5.

    Datasets:
        rates     shape (n_sites, n_d0), dtype float64
        d0        shape (n_d0,)
        site_lons shape (n_sites,)
        site_lats shape (n_sites,)

    When both ``rates_principal`` and ``rates_distributed`` are supplied
    (same shape as ``rates``), they are stored as extra datasets
    ``rates_principal`` / ``rates_distributed``; their sum equals ``rates``.

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
    split_arrs: dict[str, np.ndarray] = {}
    if rates_principal is not None and rates_distributed is not None:
        for name, comp in (("rates_principal", rates_principal),
                           ("rates_distributed", rates_distributed)):
            comp_arr = np.asarray(comp, dtype=float)
            if comp_arr.shape != expected:
                raise ValueError(
                    f"{name} shape {comp_arr.shape} != (n_sites, n_d0)={expected}"
                )
            split_arrs[name] = comp_arr
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        f.create_dataset("rates", data=rates_arr, compression="gzip")
        for name, comp_arr in split_arrs.items():
            f.create_dataset(name, data=comp_arr, compression="gzip")
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
    rates_mean_principal: Any = None,
    rates_mean_distributed: Any = None,
) -> None:
    """Write the LT-weighted per-site mean annual-rate grid as HDF5.

    When both ``rates_mean_principal`` and ``rates_mean_distributed`` are
    supplied, they are stored as extra datasets of the same shape; the mean
    is linear so their sum equals ``rates_mean``.
    """
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
    split_arrs: dict[str, np.ndarray] = {}
    if rates_mean_principal is not None and rates_mean_distributed is not None:
        for name, comp in (("rates_mean_principal", rates_mean_principal),
                           ("rates_mean_distributed", rates_mean_distributed)):
            comp_arr = np.asarray(comp, dtype=float)
            if comp_arr.shape != expected:
                raise ValueError(
                    f"{name} shape {comp_arr.shape} != (n_sites, n_d0)={expected}"
                )
            split_arrs[name] = comp_arr
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        f.create_dataset("rates_mean", data=arr, compression="gzip")
        for name, comp_arr in split_arrs.items():
            f.create_dataset(name, data=comp_arr, compression="gzip")
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
    qs: Sequence[float] | None = None,
) -> None:
    """Write per-site fractile annual-rate grids as HDF5.

    ``rates_fractiles`` has shape ``(n_quantiles, n_sites, n_d0)`` with the
    quantile order given by ``qs`` (default :data:`FRACTILE_QS`). The order is
    also recorded as an attribute ``quantiles`` and as a string dataset
    ``quantile_labels`` (OpenQuake-style ``quantile-<q>``) for self-description.
    """
    try:
        import h5py  # type: ignore
    except Exception as exc:  # pragma: no cover
        raise NotImplementedError(
            "v4 map-mode IO requires h5py; install h5py to proceed."
        ) from exc

    qs = list(FRACTILE_QS) if qs is None else list(qs)
    labels = [quantile_label(q) for q in qs]
    arr = np.asarray(rates_fractiles, dtype=float)
    expected = (len(qs), len(site_lons), len(d0))
    if arr.shape != expected:
        raise ValueError(
            f"rates_fractiles shape {arr.shape} != "
            f"(n_quantiles, n_sites, n_d0)={expected}"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        f.create_dataset("rates_fractiles", data=arr, compression="gzip")
        f.create_dataset("quantiles", data=np.asarray(qs, dtype=float))
        f.create_dataset(
            "quantile_labels",
            data=np.asarray(labels, dtype="S"),
        )
        f.create_dataset("d0", data=np.asarray(d0, dtype=float))
        f.create_dataset("site_lons", data=np.asarray(site_lons, dtype=float))
        f.create_dataset("site_lats", data=np.asarray(site_lats, dtype=float))
        f.attrs["quantiles"] = np.asarray(qs, dtype=float)


def write_displacement_map_csv(
    path: Path,
    *,
    site_lons: Sequence[float],
    site_lats: Sequence[float],
    displ: Sequence[float],
    target_return_period: float,
    label: str = "displ",
    site_is_trace: Sequence[bool] | None = None,
    displ_principal: Sequence[float] | None = None,
    displ_distributed: Sequence[float] | None = None,
) -> None:
    """Write one displacement-map CSV for a single aggregation level.

    Columns: ``site_id, lon, lat, is_trace, <label>``. A header comment
    ``# return_period = T`` records the inversion target. Used by v4
    map-mode to emit one CSV per aggregation level
    (``displacement_map_mean.csv``, ``displacement_map_p05.csv``, ...).

    When both ``displ_principal`` and ``displ_distributed`` are given,
    ``<label>_principal`` and ``<label>_distributed`` columns follow
    ``<label>``. Each component is obtained by inverting that component's
    own rate curve at the target return period, so the components answer
    "what displacement has this return period considering only
    principal/distributed faulting" - the inversion is nonlinear, hence
    they do NOT sum to the total column.
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
    split_cols: list[tuple[str, list[float]]] = []
    if displ_principal is not None and displ_distributed is not None:
        for suffix, comp in (("principal", displ_principal),
                             ("distributed", displ_distributed)):
            comp_vals = [float(v) for v in comp]
            if len(comp_vals) != n:
                raise ValueError(f"displ_{suffix} length mismatch")
            split_cols.append((f"{label}_{suffix}", comp_vals))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        f.write(f"# return_period = {float(target_return_period)}\n")
        w = csv.writer(f)
        w.writerow(["site_id", "lon", "lat", "is_trace", label]
                   + [name for name, _ in split_cols])
        for i in range(n):
            w.writerow([i, lons[i], lats[i], int(is_trace[i]), vals[i]]
                       + [col[i] for _, col in split_cols])



