from __future__ import annotations

import glob
import hashlib
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import numpy as np

from openquake.fdha.calc.config_loader import load_config, validate_public_logic_tree_ini_file
from openquake.fdha.calc.utils.parsing import parse_source_model_faults
from openquake.fdha.logic_tree.aggregation import weighted_fractiles, weighted_mean
from openquake.fdha.logic_tree.config_builder import build_config, dump_config_to_ini
from openquake.fdha.logic_tree.enumerator import SourceInfo, enumerate_end_branches
from openquake.fdha.logic_tree.io import (
    FRACTILE_LABELS,
    FRACTILE_QS,
    write_aggregate_csv,
    write_branch_rate_cube_npz,
    write_branch_rates_csv,
    write_branch_rates_h5,
    write_displacement_map_csv,
    write_hazard_displacement_map_csv,
    write_manifest,
    write_rates_fractiles_h5,
    write_rates_mean_h5,
    write_validator_report,
)
from openquake.fdha.logic_tree.nrml_reader import parse as parse_nrml
from openquake.fdha.logic_tree.source_model_lt import (
    SourceModelBranch,
    apply_realization_to_sources,
    expand_branch_paths,
    load_source_model_branches,
)
from openquake.fdha.logic_tree.validators import validate_spec


@dataclass(frozen=True)
class LogicTreeResult:
    outdir: str
    d0: list[float]
    mean_rates: list[list[float]]
    fractiles: dict[float, list[list[float]]]
    mode: str = "hazard_curve"
    site_lons: Optional[list[float]] = None
    site_lats: Optional[list[float]] = None
    site_is_trace: Optional[list[bool]] = None
    displ_mean: Optional[list[float]] = None
    displ_fractiles: Optional[dict[float, list[float]]] = None
    target_return_period: Optional[float] = None


class FdhaLogicTree:
    def __init__(self, ini_path: str | Path):
        self.ini_path = str(ini_path)
        ini_path_res = Path(self.ini_path).resolve()
        self.base_config = load_config(self.ini_path)
        validate_public_logic_tree_ini_file(ini_path_res, self.base_config)
        self.config_dir = str(ini_path_res.parent)

        calc = self.base_config.get("calculation", {})
        if not isinstance(calc, dict):
            raise ValueError("[calculation] must be a mapping")
        singular = calc.get("fdha_logic_tree_file")
        if singular:
            if isinstance(singular, str):
                self.logic_tree_files = [singular.strip()] if singular.strip() else []
            elif isinstance(singular, (list, tuple)):
                self.logic_tree_files = [str(f).strip() for f in singular if str(f).strip()]
            else:
                self.logic_tree_files = [str(singular).strip()]
        else:
            raise ValueError(
                "Missing [calculation].fdha_logic_tree_file "
                "(canonical; plural fdha_logic_tree_files is no longer supported)."
            )
        if not self.logic_tree_files:
            raise ValueError(
                "Missing or empty [calculation].fdha_logic_tree_file (canonical)."
            )

        self._source_lt_file = (
            calc.get("source_model_logic_tree_file") if isinstance(calc, dict) else None
        )
        if not self._source_lt_file:
            raise ValueError(
                "Missing [calculation].source_model_logic_tree_file "
                "(canonical SMLT wrapper)."
            )

        self.source_model_branches: list[SourceModelBranch] = (
            load_source_model_branches(self.base_config, self.config_dir)
        )
        # Default ``source_model_paths`` to the first branch's files so any
        # legacy code path that reads ``self.source_model_paths`` still has
        # a valid value (the multi-source orchestration overrides this per
        # iteration).
        self.source_model_paths = expand_branch_paths(
            self.source_model_branches[0]
        )

    @classmethod
    def from_ini(cls, path: str | Path) -> "FdhaLogicTree":
        return cls(path)

    # ------------------------------------------------------------------ main
    def run(self, outdir: Optional[str | Path] = None) -> LogicTreeResult:
        outdir = str(outdir or Path(self.config_dir) / "out")
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        calc_type = _detect_calculation_type(self.base_config)

        # Multi-source orchestration: enumerate FDHA end-branches per source
        # model, scale weights by source-model weight, then aggregate.
        if calc_type == "hazard_map":
            return self._run_map_multi_source(outdir_path)
        return self._run_curve_or_map_multi_source(outdir_path, mode="hazard_curve")

    # --------------------------------------------------------- enumeration
    def _enumerate_and_validate(self, outdir_path: Path):
        converter_params = _converter_params_from_config(self.base_config)
        faults = parse_source_model_faults(
            self.source_model_paths, hdf5path="", **converter_params
        )
        sources: list[SourceInfo] = []
        for sid, src in faults.items():
            sources.append(SourceInfo(source_id=str(sid), rake=_infer_source_rake(src)))
        source_ids = {s.source_id for s in sources}

        specs = [parse_nrml(Path(self.config_dir) / f) for f in self.logic_tree_files]
        merged = _merge_specs(specs)
        merged_report = validate_spec(merged, source_ids=source_ids)
        all_reports = [(merged, merged_report)]
        _write_validator_files(outdir_path, all_reports)
        merged_report.raise_if_errors()

        end_branches = enumerate_end_branches(merged, sources)

        # Dedup by selection-only fingerprint (sum weights of duplicates).
        # When a source file contains multiple same-style sources, the
        # enumerator emits one end-branch per (source, selection) pair with
        # weight 1/N_selections each; the total across sources is N_sources.
        # Dedup collapses same-selection branches and keeps a consolidated
        # source_ids list so downstream rate aggregation stays consistent.
        end_branches = _dedup_end_branches(end_branches)
        return end_branches, all_reports, source_ids

    # -------------------------------------------------------- hazard-curve
    def _run_curve(self, outdir_path: Path, end_branches, all_reports) -> LogicTreeResult:
        branch_rates = []
        branch_weights = []
        d0 = None
        site_lons = None
        site_lats = None

        for idx, eb in enumerate(end_branches):
            tmp_ini = self._write_branch_ini(eb, outdir_path, idx)
            rates, d0_vals, lons, lats = _run_single(tmp_ini)
            if d0 is None:
                d0 = d0_vals
                site_lons = lons
                site_lats = lats
            branch_rates.append(np.asarray(rates, dtype=float))
            branch_weights.append(float(eb.weight))

            write_branch_rates_csv(
                outdir_path / "hazard_curves" / f"branch_{idx:04d}.csv",
                d0=d0_vals, rates=rates, site_lons=lons, site_lats=lats,
            )

        rates_arr = np.stack(branch_rates, axis=0)
        w = np.asarray(branch_weights, dtype=float)
        w = w / w.sum()
        mean = weighted_mean(rates_arr, w)
        fr = weighted_fractiles(rates_arr, w)

        write_aggregate_csv(
            outdir_path / "aggregate_hazard.csv",
            d0=d0, mean_rates=mean.tolist(),
            fractiles={q: fr[q].tolist() for q in fr},
            site_lons=site_lons, site_lats=site_lats,
        )

        manifest = _build_manifest(
            end_branches, branch_weights, outdir_path, all_reports, mode="hazard_curve",
        )
        write_manifest(outdir_path / "manifest.json", manifest)

        return LogicTreeResult(
            outdir=str(outdir_path),
            d0=list(d0),
            mean_rates=mean.tolist(),
            fractiles={q: fr[q].tolist() for q in fr},
            mode="hazard_curve",
            site_lons=list(map(float, site_lons)),
            site_lats=list(map(float, site_lats)),
        )

    # ------------------------------------------------------------ map mode
    def _run_map(
        self,
        outdir_path: Path,
        end_branches,
        all_reports,
        *,
        fault_sources: Optional[dict] = None,
    ) -> LogicTreeResult:
        """Run the hazard-map FDHA loop in ``outdir_path``.

        When ``fault_sources`` is provided (typically by the SMLT
        orchestrator), the calculator uses those pre-built sources directly
        and the site builder skips XML re-parsing.  This is the same
        ``fault_sources=`` plumbing used by hazard-curve mode and is what
        lets map mode honour ``source_model_logic_tree_file`` without
        re-implementing site/grid logic per realisation.
        """
        import logging as _lg
        log = _lg.getLogger(__name__)

        cfg = self.base_config
        geom = cfg.get("geometry", cfg.get("site_location", {}))
        region = geom.get("region")
        if not region:
            raise ValueError(
                "Map-mode logic tree requires [geometry].region (e.g. "
                "'lon1 lat1, lon2 lat2, ...')."
            )
        spacing = float(geom.get("region_grid_spacing", 0.1))
        max_dist = float(geom.get("max_distance_km", 10.0))
        vs30 = cfg.get("site_location", {}).get("vs30", None)

        para = cfg.get("parameters", {})
        calc_cfg = cfg.get("calculation", {})
        return_period = float(
            calc_cfg.get("return_period", para.get("return_period", 100000))
        )
        target_displ = np.asarray(para["target_displacement"], dtype=float)
        converter_params = _converter_params_from_config(cfg)

        from openquake.fdha.logic_tree.site_builder import build_hazard_map_sites
        sites = build_hazard_map_sites(
            region=region,
            spacing=spacing,
            max_distance_km=max_dist,
            vs30=vs30,
            source_model_paths=self.source_model_paths,
            hdf5path=None,
            fault_sources=fault_sources,
            **converter_params,
        )

        # Per-site lon/lat (active grid sites first, then trace sites).
        sl_lons = np.array([s.location.longitude for s in sites.combined_sitecol],
                           dtype=float)
        sl_lats = np.array([s.location.latitude for s in sites.combined_sitecol],
                           dtype=float)
        is_trace = np.zeros(len(sl_lons), dtype=bool)
        is_trace[sites.n_active_grid:] = True

        n_sites = len(sl_lons)
        n_d0 = len(target_displ)
        n_br = len(end_branches)
        log.info(
            "LT map mode: %d end-branches x %d sites x %d D0 levels "
            "(%d active-grid + %d trace)",
            n_br, n_sites, n_d0, sites.n_active_grid, sites.n_trace,
        )

        rates_cube = np.zeros((n_br, n_sites, n_d0), dtype=float)
        branch_weights: list[float] = []
        fingerprints: list[str] = []

        from openquake.fdha.calc.calculators import BaseFaultRuptureCalculator
        from openquake.fdha.calc.hazard import calculate_fdha_hazard

        # Restrict each end-branch to its own source(s) so the selected
        # models are only applied to the style they were chosen for.  When
        # ``fault_sources`` was passed in (SMLT path) we filter the dict;
        # otherwise we fall back to splitting the XML into per-source files
        # so the legacy single-SMLT path is byte-for-byte unchanged.
        per_source_xml: dict[str, str] = {}
        if fault_sources is None:
            per_source_xml = _build_per_source_subset_xmls(
                self.source_model_paths, outdir_path / "per_source_models"
            )

        branches_dir = outdir_path / "branches"
        branches_dir.mkdir(parents=True, exist_ok=True)
        expected_branch_shape = (n_sites, n_d0)

        for idx, eb in enumerate(end_branches):
            tmp_ini = self._write_branch_ini(eb, outdir_path, idx)
            # Each end-branch carries the first source_id (dedup may collapse
            # identical selections across same-style sources; we still emit
            # one physical run per original source_id listed).
            sid_list = [s for s in eb.source_id.split(",") if s]
            if fault_sources is not None:
                sub_sources = {
                    sid: fault_sources[sid] for sid in sid_list if sid in fault_sources
                }
                if not sub_sources:
                    raise ValueError(
                        f"End-branch {idx} references sources {sid_list!r} "
                        f"but none of them are present in the realisation's "
                        f"fault_sources (available: {sorted(fault_sources)})."
                    )
                branch_calc = BaseFaultRuptureCalculator(
                    tmp_ini, [], hdf5path="",
                    fault_sources=sub_sources, **converter_params,
                )
            else:
                source_xmls = [per_source_xml[sid] for sid in sid_list]
                branch_calc = BaseFaultRuptureCalculator(
                    tmp_ini, source_xmls, hdf5path="", **converter_params,
                )
            res = calculate_fdha_hazard(
                branch_calc, sites.combined_sitecol, show_progress=False,
            )
            branch_rates = np.asarray(res["annual_rate_total"], dtype=float)
            # Explicit pre-aggregation shape guard: bail out loudly when a
            # branch produces the wrong grid rather than silently broadcasting.
            if branch_rates.shape != expected_branch_shape:
                raise ValueError(
                    f"Branch {idx} (fingerprint "
                    f"{_fingerprint_end_branch(eb)}) produced rates of shape "
                    f"{branch_rates.shape}; expected {expected_branch_shape}. "
                    "Site grid or D0 axis is inconsistent across branches."
                )
            rates_cube[idx] = branch_rates
            branch_weights.append(float(eb.weight))
            fp = _fingerprint_end_branch(eb)
            fingerprints.append(fp)

            # v4: one HDF5 per end-branch under branches/branch_XXXX.h5
            write_branch_rates_h5(
                branches_dir / f"branch_{idx:04d}.h5",
                rates=branch_rates,
                d0=target_displ.tolist(),
                site_lons=sl_lons.tolist(),
                site_lats=sl_lats.tolist(),
                weight=float(eb.weight),
                fingerprint=fp,
                attrs={
                    "source_id": eb.source_id,
                    "style": str(eb.style or ""),
                    "index": int(idx),
                },
            )
            if (idx + 1) % max(1, n_br // 10) == 0 or idx == n_br - 1:
                log.info("  branch %d/%d done", idx + 1, n_br)

        # Physically correct aggregation for independent faults:
        #   total_mean(site, D0) = sum over sources of per-source LT-mean
        # where each per-source LT-mean uses weights normalised within its
        # own group (weights already sum to 1 per group by construction).
        source_ids_per_branch = [eb.source_id for eb in end_branches]
        mean_rates, per_source_means = _aggregate_multi_source_mean(
            rates_cube, branch_weights, source_ids_per_branch,
        )
        # Fractiles: computed per source and then SUM across sources (same
        # reasoning as the mean). For sources with a single end-branch the
        # fractile equals the mean; that case collapses to the deterministic
        # contribution as expected.
        frac_rates = _aggregate_multi_source_fractiles(
            rates_cube, branch_weights, source_ids_per_branch,
        )

        # Invert each site's mean rate curve at target_rate; likewise fractiles.
        from openquake.fdha.calc.utils.interpolation import get_map_from_curves
        target_rate = 1.0 / return_period
        displ_mean = get_map_from_curves(target_displ, mean_rates, target_rate)
        displ_frac: dict[float, np.ndarray] = {}
        for q, cube in frac_rates.items():
            displ_frac[q] = get_map_from_curves(target_displ, cube, target_rate)

        # v4 aggregate persistence layout:
        #   aggregate/rates_mean.h5
        #   aggregate/rates_fractiles.h5        (5, n_sites, n_d0)
        #   aggregate/displacement_map_mean.csv
        #   aggregate/displacement_map_p{05,16,50,84,95}.csv
        agg_dir = outdir_path / "aggregate"
        agg_dir.mkdir(parents=True, exist_ok=True)
        write_rates_mean_h5(
            agg_dir / "rates_mean.h5",
            rates_mean=mean_rates,
            d0=target_displ.tolist(),
            site_lons=sl_lons.tolist(),
            site_lats=sl_lats.tolist(),
        )
        fr_stack = np.stack([frac_rates[q] for q in FRACTILE_QS], axis=0)
        write_rates_fractiles_h5(
            agg_dir / "rates_fractiles.h5",
            rates_fractiles=fr_stack,
            d0=target_displ.tolist(),
            site_lons=sl_lons.tolist(),
            site_lats=sl_lats.tolist(),
        )
        write_displacement_map_csv(
            agg_dir / "displacement_map_mean.csv",
            site_lons=sl_lons.tolist(),
            site_lats=sl_lats.tolist(),
            displ=displ_mean.tolist(),
            target_return_period=return_period,
            label="displ_mean",
            site_is_trace=is_trace.tolist(),
        )
        for q, lab in zip(FRACTILE_QS, FRACTILE_LABELS):
            write_displacement_map_csv(
                agg_dir / f"displacement_map_{lab}.csv",
                site_lons=sl_lons.tolist(),
                site_lats=sl_lats.tolist(),
                displ=displ_frac[q].tolist(),
                target_return_period=return_period,
                label=f"displ_{lab}",
                site_is_trace=is_trace.tolist(),
            )

        manifest = _build_manifest(
            end_branches, branch_weights, outdir_path, all_reports,
            mode="hazard_map",
            extra={
                "grid": {
                    "region": region,
                    "spacing": spacing,
                    "max_distance_km": max_dist,
                    "n_sites_total": int(n_sites),
                    "n_active_grid": int(sites.n_active_grid),
                    "n_trace": int(sites.n_trace),
                },
                "return_period": return_period,
                "layout": {
                    "branches_dir": "branches/",
                    "aggregate_dir": "aggregate/",
                    "rates_mean_h5": "aggregate/rates_mean.h5",
                    "rates_fractiles_h5": "aggregate/rates_fractiles.h5",
                    "displacement_map_mean_csv": "aggregate/displacement_map_mean.csv",
                    "displacement_map_fractile_csvs": [
                        f"aggregate/displacement_map_{lab}.csv"
                        for lab in FRACTILE_LABELS
                    ],
                },
                "fractile_quantiles": list(FRACTILE_QS),
                "per_source_mean_sum_strategy": (
                    "total_mean = sum_over_sources(per_source_LT_weighted_mean); "
                    "each source's weights are normalised within its own group."
                ),
                "per_source_ids": sorted({eb.source_id for eb in end_branches}),
            },
        )
        write_manifest(outdir_path / "manifest.json", manifest)

        return LogicTreeResult(
            outdir=str(outdir_path),
            d0=target_displ.tolist(),
            mean_rates=mean_rates.tolist(),
            fractiles={q: frac_rates[q].tolist() for q in frac_rates},
            mode="hazard_map",
            site_lons=sl_lons.tolist(),
            site_lats=sl_lats.tolist(),
            site_is_trace=is_trace.tolist(),
            displ_mean=displ_mean.tolist(),
            displ_fractiles={q: displ_frac[q].tolist() for q in displ_frac},
            target_return_period=return_period,
        )

    # --------------------------------------------------- multi-source path
    def _run_curve_or_map_multi_source(
        self, outdir_path: Path, mode: str,
    ) -> LogicTreeResult:
        """Run the logic tree once per source-model branch and aggregate.

        For each ``SourceModelBranch`` we temporarily swap
        ``self.source_model_paths`` to that branch's XML file, re-enumerate the
        FDHA end-branches against it, scale every end-branch weight by the
        source-model weight, and accumulate into a single combined ensemble.
        Aggregation (mean, fractiles) is then computed over the combined
        realisations exactly as in the single-source case.
        """
        if mode not in ("hazard_curve", "hazard_map"):
            raise ValueError(f"Unsupported logic-tree mode: {mode}")

        # Persist the original list so we can restore it after the loop.  We
        # also write per-source subdirectories under outdir to keep manifests
        # readable when there are several source models.
        original_paths = list(self.source_model_paths)

        combined_rates: list[np.ndarray] = []
        combined_weights: list[float] = []
        combined_records: list[dict[str, Any]] = []
        d0_ref = None
        site_lons_ref = None
        site_lats_ref = None
        all_reports_combined: list = []

        per_source_outdirs: dict[str, Path] = {}

        # Cache: base parsed sources per (set of XML files) to avoid re-parsing
        # when several realisations share the same underlying ``sourceModel``
        # branch (e.g. ``sm_a`` × {bg0, bg1, mm0, mm1}).
        converter_params = _converter_params_from_config(self.base_config)
        base_sources_cache: dict[tuple[str, ...], dict] = {}

        def _base_sources(sm_files: list[str]) -> dict:
            key = tuple(sm_files)
            if key not in base_sources_cache:
                base_sources_cache[key] = parse_source_model_faults(
                    list(sm_files), hdf5path="", **converter_params,
                )
            return base_sources_cache[key]

        for sm_idx, sm_branch in enumerate(self.source_model_branches):
            sm_files = expand_branch_paths(sm_branch)
            self.source_model_paths = list(sm_files)

            sub_outdir = outdir_path / "source_model_branches" / (
                f"{sm_idx:02d}_{_sanitize_for_path(sm_branch.branch_id)}"
            )
            sub_outdir.mkdir(parents=True, exist_ok=True)
            per_source_outdirs[sm_branch.branch_id] = sub_outdir

            end_branches, all_reports, _ = self._enumerate_and_validate(sub_outdir)
            all_reports_combined.extend(all_reports)

            # Apply NRML uncertainties (bGRRelative, dip modifications, ...)
            # exactly once per realisation; pass the resulting source dict
            # straight to the FDHA calculator via ``fault_sources``.
            base = _base_sources(sm_files)
            modified_sources = apply_realization_to_sources(sm_branch, base)

            for fdha_idx, eb in enumerate(end_branches):
                rates, d0_vals, lons, lats = self._compute_branch_rates(
                    eb,
                    sub_outdir,
                    fdha_idx,
                    mode=mode,
                    fault_sources=modified_sources,
                )
                if d0_ref is None:
                    d0_ref = d0_vals
                    site_lons_ref = lons
                    site_lats_ref = lats

                combined_w = float(sm_branch.weight) * float(eb.weight)
                combined_rates.append(np.asarray(rates, dtype=float))
                combined_weights.append(combined_w)

                global_idx = len(combined_records)
                combined_records.append(
                    {
                        "global_index": global_idx,
                        "source_model_branch_id": sm_branch.branch_id,
                        "source_model_file": sm_branch.source_model_file,
                        "source_model_weight": float(sm_branch.weight),
                        "source_model_uncertainties": [
                            {
                                "uncertainty_type": u.uncertainty_type,
                                "value": _safe_jsonable(u.value),
                                "branch_set_id": u.branch_set_id,
                                "branch_id": u.branch_id,
                                "apply_to_sources": list(u.apply_to_sources)
                                if u.apply_to_sources else None,
                                "apply_to_branches": list(u.apply_to_branches)
                                if u.apply_to_branches else None,
                            }
                            for u in sm_branch.uncertainties
                        ],
                        "fdha_branch_id": _fdha_branch_path(eb),
                        "fdha_branch_weight": float(eb.weight),
                        "combined_branch_weight": combined_w,
                        "fdha_source_id": eb.source_id,
                        "fdha_style": eb.style,
                        "fdha_models": {
                            slot: eb.selections[slot].class_name
                            for slot in eb.selections
                        },
                    }
                )

                if mode == "hazard_curve":
                    write_branch_rates_csv(
                        outdir_path / "hazard_curves" / f"branch_{global_idx:04d}.csv",
                        d0=d0_vals, rates=rates, site_lons=lons, site_lats=lats,
                    )

        # Restore for any callers that read it later.
        self.source_model_paths = original_paths

        rates_arr = np.stack(combined_rates, axis=0)
        w = np.asarray(combined_weights, dtype=float)
        total_weight = float(w.sum())
        if total_weight <= 0:
            raise ValueError("Combined branch weights sum to zero")
        w_norm = w / total_weight

        mean = weighted_mean(rates_arr, w_norm)
        fr = weighted_fractiles(rates_arr, w_norm)

        if mode == "hazard_curve":
            write_aggregate_csv(
                outdir_path / "aggregate_hazard.csv",
                d0=d0_ref, mean_rates=mean.tolist(),
                fractiles={q: fr[q].tolist() for q in fr},
                site_lons=site_lons_ref, site_lats=site_lats_ref,
            )

        manifest = self._build_multi_source_manifest(
            combined_records=combined_records,
            outdir_path=outdir_path,
            all_reports=all_reports_combined,
            mode=mode,
            total_weight=total_weight,
        )
        write_manifest(outdir_path / "manifest.json", manifest)

        return LogicTreeResult(
            outdir=str(outdir_path),
            d0=list(d0_ref) if d0_ref is not None else [],
            mean_rates=mean.tolist(),
            fractiles={q: fr[q].tolist() for q in fr},
            mode=mode,
            site_lons=list(map(float, site_lons_ref)) if site_lons_ref is not None else [],
            site_lats=list(map(float, site_lats_ref)) if site_lats_ref is not None else [],
        )

    # ---------------------------------------------- multi-source map mode
    def _run_map_multi_source(self, outdir_path: Path) -> LogicTreeResult:
        """Run the hazard-map calculator once per source-model realisation
        and aggregate the per-realisation maps with the SMLT weights.

        For each ``SourceModelBranch`` we:

        1. parse the base sources (cached per file-set);
        2. apply the realisation's NRML uncertainties to produce a modified
           source dict;
        3. enumerate FDHA end-branches against that branch's sources and run
           the same map kernel as the single-SMLT path
           (:meth:`_run_map`), using ``fault_sources=`` so the calculator
           reuses the prepared sources without re-parsing XML.

        The final mean / fractile rate cubes are the SMLT-weighted
        arithmetic mean of the per-realisation cubes.  This keeps the
        OpenQuake hierarchy: source-model logic tree on the outside, FDHA
        logic tree on the inside, with combined weight ``sm_w * fdha_w``.

        The aggregation assumes every SMLT branch produces the same site
        grid (true whenever NRML uncertainties only change MFD/dip/etc.,
        which is the typical use-case).  A clear error is raised if the
        site count drifts between branches.
        """
        converter_params = _converter_params_from_config(self.base_config)
        base_cache: dict[tuple[str, ...], dict] = {}

        sm_runs: list[dict[str, Any]] = []
        combined_records: list[dict[str, Any]] = []
        all_reports_combined: list = []
        d0_ref: Optional[list[float]] = None
        sl_lons_ref: Optional[list[float]] = None
        sl_lats_ref: Optional[list[float]] = None
        is_trace_ref: Optional[list[bool]] = None
        return_period_ref: Optional[float] = None

        original_paths = list(self.source_model_paths)
        try:
            for sm_idx, sm_branch in enumerate(self.source_model_branches):
                sm_files = expand_branch_paths(sm_branch)
                sub_outdir = (
                    outdir_path / "source_model_branches"
                    / f"{sm_idx:02d}_{_sanitize_for_path(sm_branch.branch_id)}"
                )
                sub_outdir.mkdir(parents=True, exist_ok=True)

                key = tuple(sm_files)
                if key not in base_cache:
                    base_cache[key] = parse_source_model_faults(
                        list(sm_files), hdf5path="", **converter_params,
                    )
                modified = apply_realization_to_sources(sm_branch, base_cache[key])

                # Enumerate FDHA end-branches against this SMLT branch's
                # sources so each end-branch sees the right source IDs.
                self.source_model_paths = list(sm_files)
                end_branches, all_reports, _ = self._enumerate_and_validate(sub_outdir)
                all_reports_combined.extend(all_reports)

                sm_res = self._run_map(
                    sub_outdir, end_branches, all_reports,
                    fault_sources=modified,
                )

                sm_runs.append(
                    {
                        "branch": sm_branch,
                        "result": sm_res,
                        "sub_outdir": sub_outdir,
                        "n_fdha_branches": len(end_branches),
                    }
                )

                if d0_ref is None:
                    d0_ref = list(sm_res.d0)
                    sl_lons_ref = list(sm_res.site_lons or [])
                    sl_lats_ref = list(sm_res.site_lats or [])
                    is_trace_ref = list(sm_res.site_is_trace or [])
                    return_period_ref = sm_res.target_return_period

                # Validate consistent site grid across SMLT branches; if a
                # geometry-modifying uncertainty altered the grid, flag it
                # rather than producing a silently-wrong aggregate.
                if (
                    len(sm_res.site_lons or []) != len(sl_lons_ref or [])
                    or len(sm_res.d0) != len(d0_ref or [])
                ):
                    raise ValueError(
                        "Site grid or D0 axis differs across source-model "
                        "branches; aggregating SMLT-weighted maps requires a "
                        "consistent grid (typically MFD-only uncertainties)."
                    )

                for fdha_idx, eb in enumerate(end_branches):
                    combined_records.append(
                        {
                            "global_index": len(combined_records),
                            "source_model_branch_id": sm_branch.branch_id,
                            "source_model_file": sm_branch.source_model_file,
                            "source_model_weight": float(sm_branch.weight),
                            "source_model_uncertainties": [
                                {
                                    "uncertainty_type": u.uncertainty_type,
                                    "value": _safe_jsonable(u.value),
                                    "branch_set_id": u.branch_set_id,
                                    "branch_id": u.branch_id,
                                    "apply_to_sources": list(u.apply_to_sources)
                                    if u.apply_to_sources else None,
                                    "apply_to_branches": list(u.apply_to_branches)
                                    if u.apply_to_branches else None,
                                }
                                for u in sm_branch.uncertainties
                            ],
                            "fdha_branch_id": _fdha_branch_path(eb),
                            "fdha_branch_weight": float(eb.weight),
                            "combined_branch_weight": (
                                float(sm_branch.weight) * float(eb.weight)
                            ),
                            "fdha_source_id": eb.source_id,
                            "fdha_style": eb.style,
                            "fdha_models": {
                                slot: eb.selections[slot].class_name
                                for slot in eb.selections
                            },
                            "branch_subdir": str(
                                sub_outdir.relative_to(outdir_path)
                            ),
                            "branch_h5": (
                                f"{sub_outdir.relative_to(outdir_path)}/"
                                f"branches/branch_{fdha_idx:04d}.h5"
                            ),
                        }
                    )
        finally:
            self.source_model_paths = original_paths

        if not sm_runs:
            raise ValueError("No source-model realisations enumerated.")

        # SMLT-weighted aggregation of per-SMLT mean / fractile rate cubes.
        sm_weights = np.asarray([r["branch"].weight for r in sm_runs], dtype=float)
        if sm_weights.sum() <= 0:
            raise ValueError("Source-model branch weights sum to zero.")
        sm_weights_norm = sm_weights / sm_weights.sum()

        per_sm_mean = np.stack(
            [np.asarray(r["result"].mean_rates, dtype=float) for r in sm_runs], axis=0
        )
        mean_rates = np.tensordot(sm_weights_norm, per_sm_mean, axes=(0, 0))

        frac_rates: dict[float, np.ndarray] = {}
        for q in FRACTILE_QS:
            per_sm_q = np.stack(
                [
                    np.asarray(r["result"].fractiles[q], dtype=float)
                    for r in sm_runs
                ],
                axis=0,
            )
            frac_rates[q] = np.tensordot(sm_weights_norm, per_sm_q, axes=(0, 0))

        target_displ = np.asarray(d0_ref, dtype=float)
        sl_lons = np.asarray(sl_lons_ref, dtype=float)
        sl_lats = np.asarray(sl_lats_ref, dtype=float)
        is_trace = np.asarray(is_trace_ref, dtype=bool)
        return_period = float(return_period_ref or 0)

        from openquake.fdha.calc.utils.interpolation import get_map_from_curves
        target_rate = 1.0 / return_period if return_period > 0 else 1.0
        displ_mean = get_map_from_curves(target_displ, mean_rates, target_rate)
        displ_frac = {
            q: get_map_from_curves(target_displ, frac_rates[q], target_rate)
            for q in FRACTILE_QS
        }

        # Top-level aggregate outputs (mirror the single-SMLT layout).
        agg_dir = outdir_path / "aggregate"
        agg_dir.mkdir(parents=True, exist_ok=True)
        write_rates_mean_h5(
            agg_dir / "rates_mean.h5",
            rates_mean=mean_rates,
            d0=target_displ.tolist(),
            site_lons=sl_lons.tolist(),
            site_lats=sl_lats.tolist(),
        )
        fr_stack = np.stack([frac_rates[q] for q in FRACTILE_QS], axis=0)
        write_rates_fractiles_h5(
            agg_dir / "rates_fractiles.h5",
            rates_fractiles=fr_stack,
            d0=target_displ.tolist(),
            site_lons=sl_lons.tolist(),
            site_lats=sl_lats.tolist(),
        )
        write_displacement_map_csv(
            agg_dir / "displacement_map_mean.csv",
            site_lons=sl_lons.tolist(),
            site_lats=sl_lats.tolist(),
            displ=displ_mean.tolist(),
            target_return_period=return_period,
            label="displ_mean",
            site_is_trace=is_trace.tolist(),
        )
        for q, lab in zip(FRACTILE_QS, FRACTILE_LABELS):
            write_displacement_map_csv(
                agg_dir / f"displacement_map_{lab}.csv",
                site_lons=sl_lons.tolist(),
                site_lats=sl_lats.tolist(),
                displ=displ_frac[q].tolist(),
                target_return_period=return_period,
                label=f"displ_{lab}",
                site_is_trace=is_trace.tolist(),
            )

        # Manifest mirrors the curve multi-source layout: a flat list of
        # combined realisations + per-SMLT summary so consumers can find
        # both per-realisation HDF5s and the top-level aggregate maps.
        sm_summary = [
            {
                "branch_id": r["branch"].branch_id,
                "source_model_file": r["branch"].source_model_file,
                "weight": float(r["branch"].weight),
                "n_fdha_branches": r["n_fdha_branches"],
                "subdir": str(r["sub_outdir"].relative_to(outdir_path)),
            }
            for r in sm_runs
        ]
        total_weight = float(sum(r["combined_branch_weight"] for r in combined_records))
        manifest: dict[str, Any] = {
            "mode": "hazard_map",
            "source_model_logic_tree_file": str(self._source_lt_file),
            "source_model_branches": sm_summary,
            "branches": combined_records,
            "n_realizations": len(combined_records),
            "total_combined_weight": total_weight,
            "return_period": return_period,
            "layout": {
                "branches_per_smlt_dir": "source_model_branches/",
                "aggregate_dir": "aggregate/",
                "rates_mean_h5": "aggregate/rates_mean.h5",
                "rates_fractiles_h5": "aggregate/rates_fractiles.h5",
                "displacement_map_mean_csv": "aggregate/displacement_map_mean.csv",
                "displacement_map_fractile_csvs": [
                    f"aggregate/displacement_map_{lab}.csv"
                    for lab in FRACTILE_LABELS
                ],
            },
            "fractile_quantiles": list(FRACTILE_QS),
            "smlt_aggregation_strategy": (
                "total_mean = sum_over_smlt(sm_w * smlt_total_mean); "
                "fractiles weight per-SMLT fractiles by sm_w (approximation)."
            ),
            "advisory_warnings": [
                {"code": iss.code, "message": iss.message}
                for _, report in all_reports_combined
                for iss in report.warnings
            ],
        }
        write_manifest(outdir_path / "manifest.json", manifest)

        return LogicTreeResult(
            outdir=str(outdir_path),
            d0=target_displ.tolist(),
            mean_rates=mean_rates.tolist(),
            fractiles={q: frac_rates[q].tolist() for q in frac_rates},
            mode="hazard_map",
            site_lons=sl_lons.tolist(),
            site_lats=sl_lats.tolist(),
            site_is_trace=is_trace.tolist(),
            displ_mean=displ_mean.tolist(),
            displ_fractiles={q: displ_frac[q].tolist() for q in displ_frac},
            target_return_period=return_period,
        )

    def _compute_branch_rates(
        self,
        eb,
        sub_outdir: Path,
        fdha_idx: int,
        mode: str,
        fault_sources: Optional[dict] = None,
    ):
        """Run a single FDHA end-branch and return its rates / D0 / coords.

        When ``fault_sources`` is provided, those pre-built source objects
        (typically returned by :func:`apply_realization_to_sources`) are
        passed straight into the calculator so NRML uncertainties such as
        ``bGRRelative`` or ``simpleFaultDipRelative`` are honoured without
        having to round-trip through XML.
        """
        if mode != "hazard_curve":
            # Map mode is handled by ``_run_map_multi_source``; this helper
            # only ever runs single-site curve realisations.
            raise ValueError(
                f"_compute_branch_rates only handles hazard_curve mode; got {mode!r}"
            )
        tmp_ini = self._write_branch_ini(eb, sub_outdir, fdha_idx)
        rates, d0_vals, lons, lats = _run_single(
            tmp_ini, fault_sources=fault_sources,
        )
        return rates, d0_vals, lons, lats

    def _build_multi_source_manifest(
        self,
        combined_records: list[dict[str, Any]],
        outdir_path: Path,
        all_reports,
        mode: str,
        total_weight: float,
    ) -> dict[str, Any]:
        advisory: list[dict[str, str]] = []
        for _, report in all_reports:
            for iss in report.warnings:
                advisory.append({"code": iss.code, "message": iss.message})

        # Group records by source-model branch for human-readable output.
        sm_summary: dict[str, dict[str, Any]] = {}
        for rec in combined_records:
            sid = rec["source_model_branch_id"]
            entry = sm_summary.setdefault(
                sid,
                {
                    "branch_id": sid,
                    "source_model_file": rec["source_model_file"],
                    "weight": rec["source_model_weight"],
                    "n_fdha_branches": 0,
                },
            )
            entry["n_fdha_branches"] += 1

        manifest: dict[str, Any] = {
            "mode": mode,
            "source_model_logic_tree_file": str(self._source_lt_file),
            "source_model_branches": list(sm_summary.values()),
            "branches": [
                {
                    "global_index": rec["global_index"],
                    "source_model_branch_id": rec["source_model_branch_id"],
                    "source_model_file": rec["source_model_file"],
                    "source_model_weight": rec["source_model_weight"],
                    "source_model_uncertainties": rec.get(
                        "source_model_uncertainties", []
                    ),
                    "fdha_branch_id": rec["fdha_branch_id"],
                    "fdha_branch_weight": rec["fdha_branch_weight"],
                    "combined_branch_weight": rec["combined_branch_weight"],
                    "fdha_source_id": rec["fdha_source_id"],
                    "fdha_style": rec["fdha_style"],
                    "fdha_models": rec["fdha_models"],
                    "curve_file": (
                        f"hazard_curves/branch_{rec['global_index']:04d}.csv"
                        if mode == "hazard_curve"
                        else None
                    ),
                }
                for rec in combined_records
            ],
            "n_realizations": len(combined_records),
            "total_combined_weight": float(total_weight),
            "advisory_warnings": advisory,
        }
        return manifest

    # ------------------------------------------------------------- helpers
    def _write_branch_ini(self, eb, outdir_path: Path, idx: int) -> str:
        cfg = build_config(self.base_config, eb)
        if isinstance(cfg.get("calculation"), dict):
            cfg["calculation"]["source_model_file"] = (
                self.source_model_paths[0]
                if len(self.source_model_paths) == 1
                else list(self.source_model_paths)
            )
            traces_file = cfg["calculation"].get("rank1p5_traces_file")
            if traces_file:
                traces_path = Path(str(traces_file))
                if not traces_path.is_absolute():
                    cfg["calculation"]["rank1p5_traces_file"] = str(
                        (Path(self.config_dir) / traces_path).resolve()
                    )
        ini_text = dump_config_to_ini(cfg)
        (outdir_path / "branch_configs").mkdir(parents=True, exist_ok=True)
        branch_ini = outdir_path / "branch_configs" / f"branch_{idx:04d}.ini"
        branch_ini.write_text(ini_text)
        return str(branch_ini)


# ================================================================= helpers


def _merge_specs(specs):
    from openquake.fdha.logic_tree.types import LogicTreeSpec
    if not specs:
        raise ValueError("No logic tree specs")
    logic_tree_id = specs[0].logic_tree_id
    basepath = specs[0].basepath
    levels = []
    for sp in specs:
        levels.extend(list(sp.branching_levels))
    return LogicTreeSpec(
        logic_tree_id=logic_tree_id, branching_levels=tuple(levels), basepath=basepath
    )


def _resolve_source_models(config_path: str, config: dict[str, Any]) -> list[str]:
    return _find_source_model(config_path, config)


def _detect_calculation_type(config: dict[str, Any]) -> str:
    """Return ``hazard_map`` when the job geometry defines a region."""
    geometry = config.get("geometry", config.get("site_location", {}))
    if "region" in geometry:
        return "hazard_map"
    return "hazard_curve"


def _find_source_model(config_path: str, config: dict[str, Any]) -> list[str]:
    """Find source model file(s) from a materialised branch configuration."""
    config_dir = os.path.dirname(os.path.abspath(config_path))

    calc_cfg = config.get("calculation", config.get("parameters", {}))
    source_file = calc_cfg.get("source_model_file", calc_cfg.get("source_model", ""))

    if source_file:
        source_files = [source_file] if isinstance(source_file, str) else list(source_file)
        resolved = []
        for sf in source_files:
            if os.path.exists(sf):
                resolved.append(sf)
                continue
            rel_path = os.path.join(config_dir, sf)
            if os.path.exists(rel_path):
                resolved.append(rel_path)
                continue
            raise ValueError(f"Source model file not found: {sf}")
        return resolved

    for pattern in ["source_model*.xml", "*source*.xml"]:
        matches = glob.glob(os.path.join(config_dir, pattern))
        if matches:
            return matches

    raise ValueError(
        "No source model specified. Public jobs must use "
        "[calculation].source_model_logic_tree_file = ...; materialised branch "
        "configs are the only internal path that may provide source_model_file."
    )


def _infer_source_rake(src) -> float:
    try:
        rups = list(src.iter_ruptures())
        if rups:
            rup0 = rups[0]
            if hasattr(rup0, "rake"):
                return float(rup0.rake)
    except Exception:
        pass
    return 0.0


def _run_single(config_path: str, fault_sources: Optional[dict] = None):
    """Run one FDHA hazard-curve calculation and return its rates.

    ``fault_sources``, if given, replaces the XML-driven parsing inside the
    calculator so callers can inject NRML-uncertainty-modified sources
    produced by :func:`apply_realization_to_sources`.
    """
    from openquake.fdha.calc.calculators import FaultRuptureProbabilityCalculator

    cfg = load_config(config_path)
    converter_params = _converter_params_from_config(cfg)

    calculation_type = _detect_calculation_type(cfg)
    if calculation_type == "hazard_map":
        # Per-branch map execution is handled by FdhaLogicTree._run_map;
        # _run_single is the single-site curve path and should not be
        # reached in map mode.
        raise NotImplementedError(
            "Logic tree on hazard maps is not yet implemented; "
            "use hazard-curve mode or disable the logic tree."
        )
    if fault_sources is None:
        source_model_paths = _resolve_source_models(config_path, cfg)
        calc = FaultRuptureProbabilityCalculator(
            config_path, source_model_paths, **converter_params
        )
    else:
        # Source paths still need to be a list-like for legacy code paths;
        # the calculator does not actually re-parse them when fault_sources
        # is supplied, but ``[]`` would satisfy any iterator usage.
        calc = FaultRuptureProbabilityCalculator(
            config_path, [],
            fault_sources=fault_sources,
            **converter_params,
        )
    res = calc.run()
    return res["poes"], res["imls"], res["site_lons"], res["site_lats"]


def _converter_params_from_config(cfg: dict[str, Any]) -> dict[str, Any]:
    erf_cfg = cfg.get("erf", cfg.get("parameters", {}))
    rupture_mesh_spacing = float(erf_cfg.get("rupture_mesh_spacing", 1.0))
    width_of_mfd_bin = float(erf_cfg.get("width_of_mfd_bin", 0.1))
    complex_fault_mesh_spacing = erf_cfg.get("complex_fault_mesh_spacing")
    if complex_fault_mesh_spacing is not None:
        complex_fault_mesh_spacing = float(complex_fault_mesh_spacing)
    out = dict(
        rupture_mesh_spacing=rupture_mesh_spacing, width_of_mfd_bin=width_of_mfd_bin
    )
    if complex_fault_mesh_spacing is not None:
        out["complex_fault_mesh_spacing"] = complex_fault_mesh_spacing
    return out


def _safe_jsonable(value):
    """Coerce an uncertainty ``value`` to something json.dump can serialize."""
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, (list, tuple)):
        return [_safe_jsonable(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _safe_jsonable(v) for k, v in value.items()}
    return str(value)


def _sanitize_for_path(name: str) -> str:
    """Return a filesystem-safe slug of ``name``."""
    safe = []
    for ch in name:
        if ch.isalnum() or ch in ("-", "_", "."):
            safe.append(ch)
        else:
            safe.append("_")
    return "".join(safe) or "branch"


def _fdha_branch_path(eb) -> str:
    """Compose a stable path identifier for an FDHA end-branch.

    Joins the branch IDs of every model selection so the resulting string is
    informative across slots (e.g. ``"SR0|FD|SSR|SFD"``).
    """
    parts = [eb.selections[slot].branch_id for slot in sorted(eb.selections)]
    if not parts:
        return _fingerprint_end_branch(eb)
    return "|".join(parts)


def _fingerprint_end_branch(eb) -> str:
    payload = {
        "style": eb.style,
        "selections": {
            slot: {"class": ch.class_name, "params": ch.params, "branch_id": ch.branch_id}
            for slot, ch in eb.selections.items()
        },
    }
    s = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(s).hexdigest()[:16]


def _selection_fingerprint(eb) -> str:
    """Fingerprint ignoring source_id so cross-source duplicates collapse."""
    return _fingerprint_end_branch(eb)


def _dedup_end_branches(end_branches):
    """Collapse end-branches that share the same style + model selections.

    When multiple same-style sources live in one source file, the enumerator
    emits one end-branch per (source, selection) pair. Those branches run
    *identical* calculations (the selections are the same); deduplicating
    avoids doing the same work N_sources times and keeps weights sensible.
    The surviving branch carries a ``source_ids`` attribute listing every
    source the selection applies to.
    """
    if not end_branches:
        return end_branches
    from openquake.fdha.logic_tree.types import EndBranch

    buckets: dict[str, list] = {}
    order: list[str] = []
    for eb in end_branches:
        key = _selection_fingerprint(eb)
        if key not in buckets:
            buckets[key] = []
            order.append(key)
        buckets[key].append(eb)

    deduped: list = []
    for key in order:
        group = buckets[key]
        rep = group[0]
        source_ids = sorted({g.source_id for g in group})
        total_weight = sum(g.weight for g in group)
        collapsed = EndBranch(
            source_id=",".join(source_ids),
            style=rep.style,
            selections=rep.selections,
            weight=total_weight,
        )
        deduped.append(collapsed)
    return deduped


def _build_manifest(
    end_branches,
    branch_weights,
    outdir_path: Path,
    all_reports,
    *,
    mode: str = "hazard_curve",
    extra: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    branches = []
    for idx, eb in enumerate(end_branches):
        branches.append(
            {
                "index": idx,
                "fingerprint": _fingerprint_end_branch(eb),
                "weight": float(branch_weights[idx]),
                "source_id": eb.source_id,
                "style": eb.style,
                "models": {slot: eb.selections[slot].class_name for slot in eb.selections},
                "curve_file": (
                    f"hazard_curves/branch_{idx:04d}.csv"
                    if mode == "hazard_curve"
                    else f"branches/branch_{idx:04d}.h5"
                ),
            }
        )
    advisory = []
    for _, report in all_reports:
        for iss in report.warnings:
            advisory.append({"code": iss.code, "message": iss.message})
    manifest: dict[str, Any] = {
        "mode": mode,
        "branches": branches,
        "advisory_warnings": advisory,
    }
    if extra:
        manifest.update(extra)
    return manifest


def _build_per_source_subset_xmls(
    source_model_paths: list[str], outdir: Path
) -> dict[str, str]:
    """Split the source model(s) into one XML per ``source_id``.

    Returns a mapping ``{source_id: abs_path_to_single_source_xml}`` so the
    driver can run each end-branch against exactly the source it was
    enumerated for (i.e. apply style-specific models only to their style's
    faults).
    """
    import xml.etree.ElementTree as ET

    NS_NRML = "http://openquake.org/xmlns/nrml/0.4"
    NS_GML = "http://www.opengis.net/gml"
    ET.register_namespace("", NS_NRML)
    ET.register_namespace("gml", NS_GML)

    outdir.mkdir(parents=True, exist_ok=True)
    per_source: dict[str, str] = {}
    # Supported tags under sourceModel. We extract any element that exposes
    # an ``id`` attribute (characteristic / simple / complex / area / point /
    # nonParametric etc.). The XML structure is preserved verbatim for the
    # extracted element so the existing converter handles it unchanged.
    for path in source_model_paths:
        tree = ET.parse(path)
        root = tree.getroot()
        # Iterate every element that has an ``id`` attribute within a
        # sourceModel; this is the natural granularity the engine uses.
        for src in list(root.iter()):
            sid = src.get("id") if hasattr(src, "get") else None
            if not sid:
                continue
            tag = src.tag
            # Skip the nrml/sourceModel wrappers themselves.
            if tag.endswith("}nrml") or tag.endswith("}sourceModel"):
                continue
            # Build a minimal document wrapping only this source.
            new_root = ET.Element(f"{{{NS_NRML}}}nrml")
            sm = ET.SubElement(
                new_root, f"{{{NS_NRML}}}sourceModel",
                attrib={"name": f"subset_{sid}"},
            )
            sm.append(src)
            ET.indent(new_root, space="  ")
            out_path = outdir / f"source_{sid}.xml"
            ET.ElementTree(new_root).write(
                out_path, xml_declaration=True, encoding="utf-8"
            )
            # ElementTree may emit xmlns:gml twice on the root; sanitise.
            text = out_path.read_text()
            dup = ' xmlns:gml="http://www.opengis.net/gml" xmlns:gml="http://www.opengis.net/gml"'
            if dup in text:
                text = text.replace(
                    dup, ' xmlns:gml="http://www.opengis.net/gml"'
                )
                out_path.write_text(text)
            per_source[str(sid)] = str(out_path.resolve())
    return per_source


def _aggregate_multi_source_mean(
    rates_cube: np.ndarray,
    branch_weights: list[float],
    source_ids_per_branch: list[str],
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Group branches by source_id; per-group weighted mean; sum across groups.

    Returns the total mean ``(n_sites, n_d0)`` and a dict mapping source_id
    to that source's LT mean rate cube ``(n_sites, n_d0)``.
    """
    from collections import defaultdict

    groups: dict[str, list[int]] = defaultdict(list)
    for i, sid in enumerate(source_ids_per_branch):
        groups[sid].append(i)

    total = np.zeros(rates_cube.shape[1:], dtype=float)
    per_source_means: dict[str, np.ndarray] = {}
    for sid, indices in groups.items():
        w = np.asarray([branch_weights[i] for i in indices], dtype=float)
        wsum = w.sum()
        if wsum <= 0:
            continue
        w = w / wsum
        group_rates = rates_cube[indices]
        m = weighted_mean(group_rates, w)
        per_source_means[sid] = m
        total = total + m
    return total, per_source_means


def _aggregate_multi_source_fractiles(
    rates_cube: np.ndarray,
    branch_weights: list[float],
    source_ids_per_branch: list[str],
) -> dict[float, np.ndarray]:
    """Per-source fractiles then sum across sources.

    NOTE: This is an approximation of the true multi-source fractile
    distribution (a convolution across independent sources). For a single
    dominant source per site it is exact; for sites where multiple sources
    contribute it yields a conservative (wider) spread and is recorded in the
    manifest as such.
    """
    from collections import defaultdict

    groups: dict[str, list[int]] = defaultdict(list)
    for i, sid in enumerate(source_ids_per_branch):
        groups[sid].append(i)

    qs = [0.05, 0.16, 0.5, 0.84, 0.95]
    totals: dict[float, np.ndarray] = {q: np.zeros(rates_cube.shape[1:], dtype=float) for q in qs}
    for sid, indices in groups.items():
        w = np.asarray([branch_weights[i] for i in indices], dtype=float)
        wsum = w.sum()
        if wsum <= 0:
            continue
        w = w / wsum
        group_rates = rates_cube[indices]
        fr = weighted_fractiles(group_rates, w)
        for q in qs:
            totals[q] = totals[q] + np.asarray(fr[q])
    return totals


def _write_validator_files(outdir_path: Path, all_reports) -> None:
    lines: list[str] = []
    for spec, report in all_reports:
        lines.append(f"spec={Path(spec.basepath).name}")
        for i in report.issues:
            lines.append(f"{i.level.upper()} {i.code} {i.message}")
        lines.append("")
    write_validator_report(outdir_path / "validator_report.txt", lines)
