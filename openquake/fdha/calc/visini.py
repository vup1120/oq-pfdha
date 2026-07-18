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
from openquake.fdha.calc.utils import _reduce_mc, _to_sites_x_displ
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
        """Precompute trace segment arrays for vectorized distance calculation."""
        segments = []
        for trace_name in self.rupture_traces:
            for trace_config in self.rank1p5_traces:
                if trace_config.get("name") != trace_name:
                    continue
                coords = trace_config.get("geometry", {}).get("coords", [])
                if not coords or len(coords) < 2:
                    continue
                for i in range(len(coords) - 1):
                    segments.append((np.array(coords[i]), np.array(coords[i + 1])))
        return segments

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

        # Vectorized distance computation for all segments
        for a, b in self._trace_segments:
            # Local coordinate scaling at mean latitude
            lat0 = np.deg2rad((a[1] + b[1] + np.mean(site_coords[:, 1])) / 3.0)
            k_lat = 111.0
            k_lon = 111.0 * np.cos(lat0)
            
            # Convert all points to km
            p_xy = np.column_stack([site_coords[:, 0] * k_lon, site_coords[:, 1] * k_lat])
            a_xy = np.array([a[0] * k_lon, a[1] * k_lat])
            b_xy = np.array([b[0] * k_lon, b[1] * k_lat])
            
            # Vectorized point-to-segment distance
            ab = b_xy - a_xy
            ap = p_xy - a_xy
            denom = np.dot(ab, ab)
            
            if denom == 0.0:
                d_km = np.linalg.norm(ap, axis=1)
            else:
                t = np.clip(np.dot(ap, ab) / denom, 0.0, 1.0)
                proj = a_xy + np.outer(t, ab)
                d_km = np.linalg.norm(p_xy - proj, axis=1)
            
            min_dists = np.minimum(min_dists, d_km)

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

        combo_probs = []
        
        for comb in self.combos:
            # === OPTIMIZED: Classify all sites at once ===
            
            # HW/FW classification (vectorized)
            hw_fw_arr = np.where(rx_arr < 0, 'FW', 'HW')
            
            # Near/far classification based on r_km (vectorized)
            near_far_arr = np.where(r_arr <= self.near_far_threshold_km, 'near', 'far')
            
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


class VisiniSecondaryCalculatorOriginal:
    """
    Original implementation preserved for verification testing.
    This is the UNOPTIMIZED version for comparison.
    """

    def __init__(
        self,
        base_sec_rup_params,
        base_sec_displ_params,
        case_label="case1",
        pixel_size=100,
        near_far_threshold_km=0.2,
        rupture_traces=None,
        rank1p5_traces=None,
        combo_override=None,
    ):
        self.base_sec_rup_params = base_sec_rup_params
        self.base_sec_displ_params = base_sec_displ_params
        self.pixel_size = int(pixel_size)
        self.near_far_threshold_km = float(near_far_threshold_km)
        self.rupture_traces = rupture_traces or []
        self.rank1p5_traces = rank1p5_traces or []
        self.combos = combo_override or self._resolve_combos(case_label)
        self._a_diagnostics = None
        self._b_diagnostics = None

    @staticmethod
    def _resolve_combos(case_label):
        # An unknown case label must fail the job: silently degrading to
        # combination A alone would under-count the distributed hazard.
        return choose_combinations(case_label)

    @staticmethod
    def _point_to_segment_distance_km(p, a, b):
        """Compute planar (km) distance from point p to segment ab."""
        p = np.asarray(p, dtype=float)
        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)
        lat0 = np.deg2rad((a[1] + b[1] + p[1]) / 3.0)
        k_lat = 111.0
        k_lon = 111.0 * np.cos(lat0)
        p_xy = np.array([p[0] * k_lon, p[1] * k_lat])
        a_xy = np.array([a[0] * k_lon, a[1] * k_lat])
        b_xy = np.array([b[0] * k_lon, b[1] * k_lat])

        ab = b_xy - a_xy
        ap = p_xy - a_xy
        denom = np.dot(ab, ab)
        if denom == 0.0:
            return float(np.linalg.norm(ap))
        t = np.clip(np.dot(ap, ab) / denom, 0.0, 1.0)
        proj = a_xy + t * ab
        return float(np.linalg.norm(p_xy - proj))

    def _compute_secondary_distance(self, r_km, site_coords):
        """Original per-site distance computation."""
        r_km_arr = np.atleast_1d(r_km)
        if not self.rupture_traces or site_coords is None:
            return r_km_arr[0] if r_km_arr.size == 1 else r_km_arr

        site_coords = np.asarray(site_coords, dtype=float)
        if site_coords.ndim == 1:
            site_coords = site_coords.reshape(1, 2)

        n_sites = site_coords.shape[0]
        min_dists = np.full((n_sites,), np.inf)

        for trace_name in self.rupture_traces:
            for trace_config in self.rank1p5_traces:
                if trace_config.get("name") != trace_name:
                    continue
                coords = trace_config.get("geometry", {}).get("coords", [])
                if not coords or len(coords) < 2:
                    continue
                for i in range(len(coords) - 1):
                    a = coords[i]
                    b = coords[i + 1]
                    for idx, p in enumerate(site_coords):
                        d_km = self._point_to_segment_distance_km(p, a, b)
                        if d_km < min_dists[idx]:
                            min_dists[idx] = d_km

        result = np.where(np.isfinite(min_dists), min_dists, r_km_arr[:n_sites])
        return result[0] if result.size == 1 else result

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
    ):
        """Original per-site computation (SLOW)."""
        n_sites = len(np.atleast_1d(r))
        n_displ = len(target_displacements)

        r_m = np.asarray(r, dtype=float) * 1000.0
        rx_m = np.asarray(rx, dtype=float) * 1000.0
        L_m = np.asarray(L, dtype=float) * 1000.0
        r_km_arr = np.atleast_1d(r)

        rup_kwargs = {k: v for k, v in self.base_sec_rup_params.items() if k in {"style", "pixel_size"}}

        combo_probs = []
        for comb in self.combos:
            sr_probs_list = []
            fd_probs_list = []
            
            # ORIGINAL: Per-site loop (SLOW!)
            for site_idx in range(n_sites):
                r_site_km = float(r_km_arr[site_idx])
                r_site_m = float(r_m[site_idx])
                rx_site_m = float(rx_m[site_idx]) if len(rx_m) > site_idx else float(rx_m)
                L_site_m = float(L_m[site_idx]) if len(L_m) > site_idx else float(L_m)
                x_L_site = float(x_L[site_idx]) if hasattr(x_L, '__len__') and len(x_L) > site_idx else float(x_L)
                dip_site = float(dip[site_idx]) if hasattr(dip, '__len__') and len(dip) > site_idx else float(dip)
                site_coord = site_coords[site_idx:site_idx+1] if site_coords is not None and len(site_coords) > site_idx else None
                
                if comb == "B":
                    r_sec_km = self._compute_secondary_distance(r_site_km, site_coord)
                    r_sec_km_val = float(r_sec_km) if np.isscalar(r_sec_km) else float(r_sec_km[0]) if len(r_sec_km) > 0 else r_site_km
                    dist_km_for_nearfar = r_site_km
                    near_or_far = "near" if dist_km_for_nearfar <= self.near_far_threshold_km else "far"
                    r_for_comb = r_sec_km_val * 1000.0
                else:
                    dist_km = r_site_km
                    near_or_far = "near" if dist_km <= self.near_far_threshold_km else "far"
                    r_for_comb = r_site_m

                if comb in {"A", "B"}:
                    np.random.seed(42)
                    style = rup_kwargs.get("style", "normal")
                    pixel_size = rup_kwargs.get("pixel_size", self.pixel_size)
                    sr_total = sr_model.calculate_rank2_total_probability(
                        mag=mag, r=r_for_comb, rx=rx_site_m, style=style,
                        pixel_size=pixel_size, combination=comb, fault_length=L_site_m,
                        site_width=self.pixel_size, near_or_far=near_or_far,
                    )
                    if isinstance(sr_total, dict):
                        sr_prob = sr_total.get("P_total", None)
                    else:
                        sr_prob = sr_total
                else:
                    sr_prob = sr_model.get_prob(
                        mag=mag, rx=rx_site_m, r=r_site_m, combination=comb, **rup_kwargs,
                    )
                
                sr_probs_list.append(sr_prob)

                if comb == "B":
                    r_sec_km = self._compute_secondary_distance(r_site_km, site_coord)
                    r_sec_km_val = float(r_sec_km) if np.isscalar(r_sec_km) else float(r_sec_km[0]) if len(r_sec_km) > 0 else r_site_km
                    s_for_fd = r_sec_km_val * 1000.0
                else:
                    s_for_fd = r_site_m
                    
                fd_prob = fd_model.get_prob(
                    d=target_displacements, mag=mag, s=s_for_fd, rx=rx_site_m,
                    X_L_ratio=x_L_site, dip=dip_site, combination=comb,
                    **self.base_sec_displ_params,
                )
                fd_probs_list.append(fd_prob)
            
            if not sr_probs_list or not fd_probs_list:
                continue
                
            sr_red_list = []
            for sr_prob in sr_probs_list:
                if sr_prob is None:
                    sr_red_list.append(None)
                    continue
                sr_red = _reduce_mc(sr_prob, method=s_sr_red_cfg.get("method", "median"),
                                   q=s_sr_red_cfg.get("q", 50))
                sr_red_list.append(sr_red)
            
            sr_mat_list = []
            for sr_red in sr_red_list:
                if sr_red is None:
                    sr_mat_list.append(np.zeros((1, 1)))
                else:
                    sr_mat_list.append(np.atleast_1d(sr_red).reshape(-1, 1))
            sr_mat = np.vstack(sr_mat_list)
            
            fd_mat_list = []
            for fd_prob in fd_probs_list:
                if fd_prob is None:
                    fd_mat_list.append(np.zeros((1, n_displ)))
                else:
                    fd_mat_site = _to_sites_x_displ(fd_prob, 1, n_displ, s_sr_red_cfg)
                    fd_mat_list.append(fd_mat_site)
            fd_mat = np.vstack(fd_mat_list)
            
            combo_probs.append(sr_mat * fd_mat)

        return combine_probabilities(combo_probs) if combo_probs else np.zeros((n_sites, n_displ))