# -*- coding: utf-8 -*-
"""
Optimized Visini (2025) Secondary Rupture/Displacement Calculator.

Performance Optimizations:
1. Batch processing by (HW/FW, near/far) groups instead of per-site loops
2. Vectorized P_slice computation
3. Single MC call per group (leveraging cache in visini2025.py)
4. Vectorized FD model calls where possible
5. Precomputed near/far classification

Scientific Correctness:
- All mathematical operations are equivalent to the original
- Results are identical within floating-point precision
- Monte Carlo reproducibility maintained via fixed seed
"""

import numpy as np
from openquake.fdha.calc.utils import _to_sites_x_displ
from openquake.fdha.calc.decision_tree import choose_combinations, combine_probabilities


class VisiniSecondaryCalculator:
    """
    Optimized encapsulation of Visini (2025) secondary rupture/displacement computation.

    Key optimizations over original:
    - Sites processed in batches by (HW/FW, near/far) groups
    - P_slice computed vectorized over all sites
    - P_along cached per rupture (only 4 unique values per rupture)
    - FD model called with full site arrays
    """

    def __init__(
        self,
        base_sec_rup_params,
        base_sec_displ_params,
        case_label="case1",
        pixel_size=100,
        along_strike_width=None,
        near_far_threshold_km=0.2,
        rupture_traces=None,
        rank1p5_traces=None,
        combo_override=None,
        segment_sampling="truncated",
        distribution_type="uniform",
    ):
        self.base_sec_rup_params = base_sec_rup_params
        self.base_sec_displ_params = base_sec_displ_params
        self.pixel_size = int(pixel_size)
        # along_strike_width: if None, defaults to pixel_size (square site)
        self.along_strike_width = float(along_strike_width) if along_strike_width is not None else float(pixel_size)
        self.near_far_threshold_km = float(near_far_threshold_km)
        self.rupture_traces = rupture_traces or []
        self.rank1p5_traces = rank1p5_traces or []
        self.combos = combo_override or self._resolve_combos(case_label)
        self._a_diagnostics = None
        self._b_diagnostics = None

        # Monte Carlo parameters for segment sampling and distribution type
        self.segment_sampling = segment_sampling
        self.distribution_type = distribution_type

        # Precompute trace segment arrays for vectorized distance calculation
        self._trace_segments = self._precompute_trace_segments()

    def _precompute_trace_segments(self):
        """Per-trace (lon, lat) polylines for the vectorized distance
        calculation (one array per configured rank-1.5 trace)."""
        polylines = []
        for trace_name in self.rupture_traces:
            for trace_config in self.rank1p5_traces:
                if trace_config.get("name") != trace_name:
                    continue
                coords = trace_config.get("geometry", {}).get("coords", [])
                if coords and len(coords) >= 2:
                    polylines.append(np.asarray(coords, dtype=float))
        return polylines

    @staticmethod
    def _resolve_combos(case_label):
        # An unknown case label must fail the job: silently degrading to
        # combination A alone would under-count the distributed hazard.
        return choose_combinations(case_label)

    def _compute_secondary_distance_vectorized(self, r_km_arr, site_coords):
        """
        OPTIMIZED: Vectorized distance computation for all sites at once.

        Computes minimum distance (km) from each site to secondary rupture traces.

        Parameters:
        -----------
        r_km_arr : array
            Fallback distances (km) if no traces found
        site_coords : array of shape (n_sites, 2)
            Site coordinates as (lon, lat)

        Returns:
        --------
        array
            Minimum distance (km) for each site
        """
        r_km_arr = np.atleast_1d(r_km_arr)

        if not self._trace_segments or site_coords is None:
            return r_km_arr

        site_coords = np.asarray(site_coords, dtype=float)
        if site_coords.ndim == 1:
            site_coords = site_coords.reshape(1, 2)

        n_sites = site_coords.shape[0]
        min_dists = np.full(n_sites, np.inf)

        # Same polyline-distance machinery (and hazardlib
        # OrthographicProjection frame) as the 'segments' reference line -
        # one vectorized call per rank-1.5 trace, replacing a per-segment
        # Python loop in a hand-rolled 111.0 km/deg flat frame.
        from openquake.fdha.calc.utils.segments import _dist_to_polyline_km
        for poly in self._trace_segments:
            min_dists = np.minimum(min_dists, _dist_to_polyline_km(
                site_coords[:, 0], site_coords[:, 1], poly[:, 0], poly[:, 1]))

        # Fallback to main r_km where no trace found
        result = np.where(np.isfinite(min_dists), min_dists, r_km_arr[:n_sites])
        return result

    def compute(
        self,
        mag,
        r,
        rx,
        L,
        x_L,
        dip,
        target_displacements,
        sr_model,
        fd_model,
        s_sr_red_cfg,
        site_coords=None,
        style=None,
    ):
        """
        OPTIMIZED: Return secondary (distributed) contribution matrix of shape (n_sites, n_displ).

        ``style`` ('normal' | 'reverse') selects the Visini coefficient set
        for both the SR and FD models. The hazard pipeline always passes the
        resolved value (explicit model parameter, else derived from the
        rupture rake); when omitted (direct callers), the 'style' entry of
        ``base_sec_rup_params`` is used, else 'normal' for backward
        compatibility.

        Key optimizations:
        - Sites grouped by (HW/FW, near/far) for batched MC
        - P_slice vectorized over all sites
        - FD model called vectorized
        """
        # Convert inputs to arrays
        r_arr = np.atleast_1d(np.asarray(r, dtype=float))
        rx_arr = np.atleast_1d(np.asarray(rx, dtype=float))
        L_arr = np.atleast_1d(np.asarray(L, dtype=float))
        x_L_arr = np.atleast_1d(np.asarray(x_L, dtype=float))
        dip_arr = np.atleast_1d(np.asarray(dip, dtype=float))
        target_displacements = np.asarray(target_displacements, dtype=float)

        n_sites = len(r_arr)
        n_displ = len(target_displacements)

        # Distances in meters for model calls
        r_m = r_arr * 1000.0
        rx_m = rx_arr * 1000.0
        L_m = L_arr * 1000.0

        # For fault length, use the first value (constant per rupture)
        fault_length_m = float(L_m[0]) if L_m.size > 0 else 10000.0

        # Extract model kwargs
        rup_kwargs = {k: v for k, v in self.base_sec_rup_params.items() if k in {"style", "pixel_size"}}
        if style is None:
            style = rup_kwargs.get("style", "normal")
        pixel_size = rup_kwargs.get("pixel_size", self.pixel_size)

        # Site classifications are combination-invariant: compute once.
        # HW/FW from the rx sign, near/far from r against the threshold.
        hw_fw_arr = np.where(rx_arr < 0, 'FW', 'HW')
        near_far_arr = np.where(r_arr <= self.near_far_threshold_km, 'near', 'far')

        combo_probs = []

        for comb in self.combos:
            # Initialize r_sec_km to r_arr as fallback
            r_sec_km = r_arr.copy()

            # For combination B, compute secondary distances for all sites at once
            if comb == "B" and self._trace_segments:
                r_sec_km = self._compute_secondary_distance_vectorized(r_arr, site_coords)
                r_for_pslice_m = r_sec_km * 1000.0
            else:
                r_for_pslice_m = r_m

            # === OPTIMIZED: Vectorized P_slice for all sites ===
            if comb in {"A", "B"}:
                # Use the optimized vectorized method from the model
                if hasattr(sr_model, 'calculate_rank2_total_probability_vectorized'):
                    # Fully vectorized path
                    sr_result = sr_model.calculate_rank2_total_probability_vectorized(
                        mag=mag,
                        r_array=r_for_pslice_m,
                        rx_array=rx_m,
                        style=style,
                        across_strike_width=pixel_size,
                        combination=comb,
                        fault_length=fault_length_m,
                        along_strike_width=self.along_strike_width,
                        near_far_array=near_far_arr,
                        distribution_type=self.distribution_type,
                        segment_sampling=self.segment_sampling,
                    )
                    P_slice = np.atleast_1d(sr_result["P_slice"])
                    P_along = np.atleast_1d(sr_result["P_along_strike"])
                    P_total_sr = np.atleast_1d(sr_result["P_total"])

                    # Ensure arrays are the right size
                    if P_slice.size == 1:
                        P_slice = np.full(n_sites, float(P_slice[0]))
                    if P_along.size == 1:
                        P_along = np.full(n_sites, float(P_along[0]))
                    if P_total_sr.size == 1:
                        P_total_sr = np.full(n_sites, float(P_total_sr[0]))
                else:
                    # Fallback: Group by (hw_fw, near_far) - still much faster than per-site
                    P_slice = sr_model.get_prob_slice(mag, r_for_pslice_m, rx_m, style, pixel_size, comb)
                    P_slice = np.atleast_1d(P_slice)
                    if P_slice.size == 1:
                        P_slice = np.full(n_sites, float(P_slice[0]))

                    P_along = np.zeros(n_sites, dtype=float)
                    for hw_fw in ['HW', 'FW']:
                        for near_far in ['near', 'far']:
                            mask = (hw_fw_arr == hw_fw) & (near_far_arr == near_far)
                            if not np.any(mask):
                                continue

                            # Single MC call for this group (cached in model)
                            # Uses across_strike_width for F-ratio lookup,
                            # along_strike_width for Monte Carlo site window
                            p_along_group = sr_model.monte_carlo_rank2_occurrence(
                                fault_length_m,
                                self.pixel_size,  # across_strike_width
                                hw_fw,
                                style,
                                near_far,
                                distribution_type=self.distribution_type,
                                segment_sampling=self.segment_sampling,
                                along_strike_width=self.along_strike_width
                            )
                            P_along[mask] = p_along_group

                    P_total_sr = P_slice * P_along

                # Store diagnostics for first site
                if n_sites > 0:
                    if comb == "A":
                        self._a_diagnostics = {
                            "r_km": float(r_arr[0]),
                            "r_m": float(r_for_pslice_m[0]),
                            "rx_m": float(rx_m[0]),
                            "near_or_far": near_far_arr[0],
                            "P_slice": float(P_slice[0]) if hasattr(P_slice, '__len__') else float(P_slice),
                            "P_along": float(P_along[0]) if hasattr(P_along, '__len__') else float(P_along),
                            "P_total": float(P_total_sr[0]) if hasattr(P_total_sr, '__len__') else float(P_total_sr),
                        }
                    elif comb == "B":
                        self._b_diagnostics = {
                            "r_km": float(r_arr[0]),
                            "r_sec_km": float(r_sec_km[0]),
                            "r_m": float(r_for_pslice_m[0]),
                            "rx_m": float(rx_m[0]),
                            "near_or_far": near_far_arr[0],
                            "P_slice": float(P_slice[0]) if hasattr(P_slice, '__len__') else float(P_slice),
                            "P_along": float(P_along[0]) if hasattr(P_along, '__len__') else float(P_along),
                            "P_total": float(P_total_sr[0]) if hasattr(P_total_sr, '__len__') else float(P_total_sr),
                        }

                sr_prob_arr = np.atleast_1d(P_total_sr)
                if sr_prob_arr.size == 1:
                    sr_prob_arr = np.full(n_sites, float(sr_prob_arr[0]))
            else:
                # Combination C: direct call (already vectorized in model)
                c_kwargs = dict(rup_kwargs)
                c_kwargs["style"] = style
                c_kwargs.setdefault("pixel_size", pixel_size)
                sr_prob = sr_model.get_prob(
                    mag=mag, rx=rx_m, r=r_m, combination=comb, **c_kwargs
                )
                sr_prob_arr = np.atleast_1d(sr_prob)
                if sr_prob_arr.size == 1:
                    sr_prob_arr = np.full(n_sites, float(sr_prob_arr[0]))

            # === OPTIMIZED: Vectorized FD model call ===
            if comb == "B" and self._trace_segments:
                # r_sec_km was computed above for combination B
                s_for_fd = r_sec_km * 1000.0  # Use secondary distance
            else:
                s_for_fd = r_m

            # Call FD model vectorized over sites. The FD coefficients are
            # style-specific too: pass the resolved style unless the FD
            # parameters already pin one explicitly.
            fd_kwargs = dict(self.base_sec_displ_params)
            fd_kwargs.setdefault("style", style)
            fd_prob = fd_model.get_prob(
                d=target_displacements,
                mag=mag,
                s=s_for_fd,
                rx=rx_m,
                X_L_ratio=x_L_arr,
                dip=dip_arr,
                combination=comb,
                **fd_kwargs,
            )

            # Reshape SR probabilities: sr_prob_arr is already (n_sites,) per-site values
            # No MC reduction needed since we computed per-site values directly
            sr_mat = sr_prob_arr.reshape(n_sites, 1)

            fd_mat = _to_sites_x_displ(fd_prob, n_sites, n_displ, s_sr_red_cfg)

            combo_probs.append(sr_mat * fd_mat)

        return combine_probabilities(combo_probs) if combo_probs else np.zeros((n_sites, n_displ))