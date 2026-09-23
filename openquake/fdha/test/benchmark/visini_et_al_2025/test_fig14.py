"""
Benchmark test: Visini et al. (2025) Figure 14 reproduction

Reproduces Figure 14 from Visini et al. (2025) paper.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")  # Use non-interactive backend for headless environments
import matplotlib.pyplot as plt
import os
import pytest
from pathlib import Path
from scipy.stats import norm as _norm_ppf, truncnorm
from openquake.pfd.secondary_surf_rup.visini2025 import Visini2025SecondarySR

pytestmark = pytest.mark.benchmark



# Create Figure 14 reproduction with only Mw 6.5
def create_figure14_mw65():
    model = Visini2025SecondarySR()
    # Fix RNG seed for deterministic Monte Carlo components
    seed = int(os.environ.get("VISINI_SEED", "42"))
    np.random.seed(seed)
    
    # Parameters matching Figure 14
    magnitude = 6.5
    site_dims = [100, 500]  # meters
    num_sims = int(os.environ.get("VISINI_NUM_SIMS", "2000"))  # reduce for speed by default
    
    # Create distance array
    distances = np.linspace(-3000, 3000, 61)  # -3km to +3km
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    
    # Normal faulting plots
    for site_dim in site_dims:
        probs = []
        
        for d in distances:
            rx = d  # rx = d for this setup
            r = abs(d)  # r is always positive
            # Determine near or far
            near_or_far = 'near' if r <= 200 else 'far'
            
            # Calculate probability - using Combination A
            result = model.calculate_rank2_total_probability(
                mag=magnitude,
                r=r,
                rx=rx,
                style='normal',
                pixel_size=site_dim,
                combination='A',
                fault_length=10000,  # 10 km
                site_width=site_dim,
                near_or_far=near_or_far,
                num_simulations=num_sims
            )
            
            probs.append(result["P_total"])
        
        # Plot
        if site_dim == 100:
            ax1.plot(distances/1000, probs, 'b-', linewidth=2.5, label=f'Site dim: {site_dim}m')
        else:
            ax2.plot(distances/1000, probs, 'b-', linewidth=2.5, label=f'Site dim: {site_dim}m')

    # Diagnostics: decompose P_total into P_slice and P_along for normal, site_dim=500
    try:
        site_dim = 500
        sample_ds = np.array([-3000, -2000, -1000, 0, 1000, 2000, 3000], dtype=float)
        P_slice_list, P_along_list, P_total_list = [], [], []
        for d in sample_ds:
            rx = d
            r = abs(d)
            near_or_far = 'near' if r <= 200 else 'far'
            res = model.calculate_rank2_total_probability(
                mag=magnitude,
                r=r,
                rx=rx,
                style='normal',
                pixel_size=site_dim,
                combination='A',
                fault_length=10000,
                site_width=site_dim,
                near_or_far=near_or_far,
                num_simulations=num_sims
            )
            P_slice_list.append(res["P_slice"])   # depends on r, rx
            P_along_list.append(res["P_along_strike"])  # depends on near/far + HW/FW
            P_total_list.append(res["P_total"])   # product
        print("\n[Diagnostics] Normal, site=500m:")
        print("  distances (m):", sample_ds.tolist())
        print("  P_slice:", [round(x, 4) for x in P_slice_list])
        print("  P_along:", [round(x, 4) for x in P_along_list])
        print("  P_total:", [round(x, 4) for x in P_total_list])
    except Exception as e:
        print(f"Diagnostics failed: {e}")
    
    # Reverse faulting plots
    for site_dim in site_dims:
        probs = []
        
        for d in distances:
            rx = d
            r = abs(d)
            near_or_far = 'near' if r <= 200 else 'far'
            
            result = model.calculate_rank2_total_probability(
                mag=magnitude,
                r=r,
                rx=rx,
                style='reverse',
                pixel_size=site_dim,
                combination='A',
                fault_length=10000,
                site_width=site_dim,
                near_or_far=near_or_far,
                num_simulations=num_sims
            )
            
            probs.append(result["P_total"])
        
        if site_dim == 100:
            ax3.plot(distances/1000, probs, 'r-', linewidth=2.5, label=f'Site dim: {site_dim}m')
        else:
            ax4.plot(distances/1000, probs, 'r-', linewidth=2.5, label=f'Site dim: {site_dim}m')
    
    # Overlay reference curve for normal faulting, 500 m site (from CSV)
    # Reference data path (place in reference_data subdirectory or update path)
    try:
        ref_path = str(Path(__file__).parent / 'reference_data' / 'visini2025_fig14_normal_500m.csv')
        ref_arr = np.loadtxt(ref_path, delimiter=',')
        ref_x_km = ref_arr[:, 0]
        ref_y = ref_arr[:, 1]
        ax2.plot(ref_x_km, ref_y, 'k--', linewidth=2.0, label='Reference (paper) 500m')
        ax2.legend()
        print(f"Overlayed reference from: {ref_path}")
    except Exception as e:
        print(f"Failed to overlay reference CSV: {e}")
    
    # Format plots
    for ax, title in zip([ax1, ax2, ax3, ax4], 
                         ['Normal faulting - Site dim 100m', 'Normal faulting - Site dim 500m', 
                          'Reverse faulting - Site dim 100m', 'Reverse faulting - Site dim 500m']):
        ax.set_xlabel('Distance (km)')
        ax.set_ylabel('Probability of Distributed Surface Rupture')
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.legend()
        ax.set_xlim(-3, 3)
        ax.set_ylim(0, 0.60)
        
        # Add vertical line at x=0
        ax.axvline(x=0, color='black', linestyle='--', alpha=0.5, linewidth=1)
        
        # Shade near-fault zone
        ax.axvspan(-0.2, 0.2, alpha=0.2, color='gray', label='Near fault (±200m)')
        
        # Add text for magnitude
        ax.text(0.02, 0.98, f'Mw = {magnitude}', transform=ax.transAxes, 
                verticalalignment='top', fontsize=12, 
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.suptitle('Figure 14 Reproduction - Visini et al. (2025)\nMw 6.5 - Combination A', fontsize=14)
    plt.tight_layout()
    # Save plot under test/figures directory for consistency
    fig_dir = Path(__file__).resolve().parents[1] / 'figures'
    fig_dir.mkdir(parents=True, exist_ok=True)
    out_png = fig_dir / 'figure_14_visini2025_style.png'
    plt.savefig(out_png, bbox_inches='tight', dpi=300)
    print(f"Saved figure to: {out_png}")
    
    # Print some key values for verification
    print("\nKey probability values at selected distances (Mw 6.5, Combination A):")
    print("\nNormal faulting, 100m site:")
    for d in [0, 200, 500, 1000, 2000]:
        rx = d
        r = abs(d)
        near_or_far = 'near' if r <= 200 else 'far'
        result = model.calculate_rank2_total_probability(
            mag=magnitude, r=r, rx=rx, style='normal', pixel_size=100,
            combination='A', fault_length=10000, site_width=100,
            near_or_far=near_or_far, num_simulations=num_sims
        )
        print(f"  Distance {d}m: P_total = {result['P_total']:.4f}")

def print_diagnostics_normal_500m():
    model = Visini2025SecondarySR()
    np.random.seed(42)
    mag = 6.5
    site_dim = 500
    fault_length = 10000
    sample_ds = np.array([0, 500, 1000, 1500, 2000, 2500, 3000], dtype=float)
    print("\nd(m) | P_slice | P_along | P_total")
    for d in sample_ds:
        rx = d
        r  = abs(d)
        near_or_far = 'near' if r <= 200 else 'far'
        res = model.calculate_rank2_total_probability(
            mag=mag, r=r, rx=rx, style='normal', pixel_size=site_dim,
            combination='A', fault_length=fault_length, site_width=site_dim,
            near_or_far=near_or_far, num_simulations=2000
        )
        print(f"{int(d):>4} | {res['P_slice']:.4f} | {res['P_along_strike']:.4f} | {res['P_total']:.4f}")

def _compute_pal_along_normal(fault_length, site_width, hw_or_fw, mechanism, near_or_far, num_simulations=2000, seed=42):
    import numpy as np
    rng = np.random.default_rng(seed)

    # F-ratio (97.5th pctl) same as model
    _f_ratio = {
        "normal": {
            "HW": {"near": {10:0.02052,20:0.03171,50:0.05298,100:0.07413,200:0.09771,500:0.13565},
                    "far":  {10:0.00201,20:0.00398,50:0.01035,100:0.02005,200:0.03689,500:0.04013}},
            "FW": {"near": {10:0.00800,20:0.01197,50:0.02148,100:0.02941,200:0.03619,500:0.05050},
                    "far":  {10:0.00082,20:0.00153,50:0.00388,100:0.00797,200:0.01539,500:0.01712}},
        },
        "reverse": {}
    }

    # Lognormal params (to derive t1/t2 as 16th/84th, then build Normal(muN, sigmaN))
    logn_params = {
        "Normal": {"HW": (3.546, 1.358), "FW": (3.390, 1.472)},
        "Reverse": {"HW": (3.622, 1.589), "FW": (3.801, 1.542)}
    }

    L = int(fault_length)
    closest = 500  # this function used only for site_width=500 in diagnostics
    F = _f_ratio[mechanism.lower()][hw_or_fw][near_or_far][closest]
    total_DR_length = fault_length * F

    # Derive [t1,t2] from log-normal params (16th,84th) then Normal mu/sigma per MATLAB SL=1
    mu_ln, sigma_ln = logn_params[mechanism][hw_or_fw]
    z16, z84 = -0.994457883, 0.994457883  # approx ppf(0.16), ppf(0.84)
    t1 = float(np.exp(mu_ln + sigma_ln * z16))
    t2 = float(np.exp(mu_ln + sigma_ln * z84))
    muN = 0.5 * (t1 + t2)
    sigmaN = 0.5 * (t2 - t1)
    a_std = (0 - muN) / sigmaN  # will be clipped by truncation at [t1,t2]
    # Use exact bounds for truncnorm in x-space via floor/ceil below

    site_lo = int(fault_length/2 - site_width/2)
    site_hi = int(fault_length/2 + site_width/2)

    hits_uniform = 0
    hits_exponential = 0

    for _ in range(num_simulations):
        segments = []
        total = 0.0
        # truncated Normal sampling on [t1,t2]
        a = (t1 - muN) / sigmaN
        b = (t2 - muN) / sigmaN
        while total < total_DR_length and len(segments) < 1000:
            seg_len = truncnorm.rvs(a, b, loc=muN, scale=sigmaN, random_state=rng)
            seg_len = float(np.clip(seg_len, t1, t2))
            segments.append(seg_len)
            total += seg_len
        if not segments:
            continue

        num_segments = len(segments)
        # UNIFORM placement with overlap fix (1..L)
        space_positions = rng.permutation(np.arange(1, L+1))[:num_segments]
        centro_unif = np.sort(space_positions)
        semi = np.asarray(segments) / 2.0
        cum_semi = np.cumsum(semi)
        for g in range(1, len(centro_unif)):
            thr = int(cum_semi[g-1])
            if (centro_unif[g] - centro_unif[g-1]) < thr:
                centro_unif[g] = centro_unif[g] + thr

        site_hit_unif = False
        for j, seg_len in enumerate(segments):
            center = centro_unif[j]
            semi_len = seg_len / 2.0
            start = max(1, int(np.floor(center - semi_len)))
            end = min(L, int(np.ceil(center + semi_len)))
            if end >= start and (end >= site_lo) and (start <= site_hi):
                site_hit_unif = True
                break
        if site_hit_unif:
            hits_uniform += 1

        # EXPONENTIAL clustered placement
        mean_distance = L / num_segments
        starting_pos = rng.integers(1, L+1)
        interdist = rng.exponential(mean_distance, num_segments)
        ini_s = starting_pos
        site_hit_exp = False
        for j in range(num_segments):
            end_s = ini_s + segments[j]
            if ini_s > L or end_s > L:
                ini_s = ini_s - L
                end_s = end_s - L
            start_i = int(np.floor(ini_s))
            end_i = int(np.ceil(end_s))
            if end_i >= start_i and (end_i >= site_lo) and (start_i <= site_hi):
                site_hit_exp = True
            ini_s = end_s + int(interdist[j])
        if site_hit_exp:
            hits_exponential += 1

    prob_uniform = hits_uniform / num_simulations if num_simulations > 0 else 0.0
    prob_exponential = hits_exponential / num_simulations if num_simulations > 0 else 0.0
    return 0.5 * (prob_uniform + prob_exponential)

def compare_along_distribution_normal_vs_current():
    print("\n[Compare P_along] normal(derived bounds) vs current(lognormal)")
    model = Visini2025SecondarySR()
    mag = 6.5
    site_dim = 500
    fault_length = 10000
    for d in [1000, 2000, 3000]:  # far-field HW
        rx = d; r = abs(d)
        pal_log = model.monte_carlo_rank2_occurrence(fault_length, r, site_dim, 'HW', 'Normal', 'far', 2000)
        pal_log_val = pal_log["P_uniform"] if isinstance(pal_log, dict) else pal_log
        pal_nrm = _compute_pal_along_normal(fault_length, site_dim, 'HW', 'Normal', 'far', 2000, seed=42)
        print(f"HW d={d}m -> lognormal: {pal_log_val:.4f}, normal: {pal_nrm:.4f}")
    for d in [-1000, -2000, -3000]:  # far-field FW
        rx = d; r = abs(d)
        pal_log = model.monte_carlo_rank2_occurrence(fault_length, r, site_dim, 'FW', 'Normal', 'far', 2000)
        pal_log_val = pal_log["P_uniform"] if isinstance(pal_log, dict) else pal_log
        pal_nrm = _compute_pal_along_normal(fault_length, site_dim, 'FW', 'Normal', 'far', 2000, seed=42)
        print(f"FW d={abs(d)}m -> lognormal: {pal_log_val:.4f}, normal: {pal_nrm:.4f}")

# Run the reproduction
if __name__ == '__main__':
    create_figure14_mw65()
    print_diagnostics_normal_500m()
    compare_along_distribution_normal_vs_current()