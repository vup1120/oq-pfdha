"""
Plot Moss and Ross (2011) Figure 6 comparisons against the current code.

Usage:
    python openquake/fdha/test/unit/plot_moss_ross2011_fig6_comparison.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from openquake.pfd.primary_surf_displ import MossRoss2011PrimaryFD


FIG6_DIGITIZED = {
    5.5: {
        "displacement_m": np.array([
            0.01, 0.02, 0.03, 0.05, 0.07, 0.10, 0.15, 0.20, 0.25,
            0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 1.00, 1.20, 1.50,
            1.70,
        ]),
        "prob_exceed": np.array([
            0.998, 0.980, 0.978, 0.974, 0.936, 0.881, 0.784, 0.685,
            0.586, 0.501, 0.367, 0.258, 0.183, 0.129, 0.090, 0.052,
            0.030, 0.012, 0.008,
        ]),
    },
    6.5: {
        "displacement_m": np.array([
            0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50, 0.70,
            1.00, 1.50, 2.00, 3.00, 5.00, 10.00,
        ]),
        "prob_exceed": np.array([
            1.00, 1.00, 0.99, 0.97, 0.88, 0.79, 0.60, 0.44,
            0.26, 0.11, 0.06, 0.02, 0.00, 0.00,
        ]),
    },
    7.5: {
        "displacement_m": np.array([
            0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50, 0.70,
            1.00, 1.50, 2.00, 3.00, 5.00, 10.00,
        ]),
        "prob_exceed": np.array([
            1.00, 1.00, 1.00, 1.00, 0.97, 0.94, 0.86, 0.78,
            0.65, 0.42, 0.28, 0.15, 0.04, 0.00,
        ]),
    },
}

STYLES = {
    5.5: {"color": "tab:green", "ls": "--", "label": "M=5.5"},
    6.5: {"color": "tab:orange", "ls": "-.", "label": "M=6.5"},
    7.5: {"color": "tab:blue", "ls": "-", "label": "M=7.5"},
}


def _configure_matplotlib():
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 12,
        "axes.linewidth": 1.2,
    })


def _compute_curves():
    mags = [5.5, 6.5, 7.5]
    model = MossRoss2011PrimaryFD()
    d_curve = np.logspace(-2, 1, 500)
    p_curve = {
        mag: model.get_prob(
            d=d_curve,
            X_L_ratio=[0.25],
            mag=mag,
            norm_disp_type="AD",
        )[:, 0]
        for mag in mags
    }
    residuals = {}
    for mag in mags:
        d = FIG6_DIGITIZED[mag]["displacement_m"]
        p = FIG6_DIGITIZED[mag]["prob_exceed"]
        code = model.get_prob(d=d, X_L_ratio=[0.25], mag=mag, norm_disp_type="AD")[:, 0]
        residuals[mag] = code - p
    return mags, d_curve, p_curve, residuals


def _save_paper_vs_code(outdir, mags, d_curve, p_curve):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True, constrained_layout=True)

    ax = axes[0]
    for mag in mags:
        digitized = FIG6_DIGITIZED[mag]
        style = STYLES[mag]
        ax.plot(
            digitized["displacement_m"],
            digitized["prob_exceed"],
            color="black",
            lw=2.0,
            ls=style["ls"],
            label=style["label"],
        )
        ax.scatter(digitized["displacement_m"], digitized["prob_exceed"], s=14, color="black", zorder=3)
    ax.set_title("Digitized paper Figure 6")
    ax.set_xlabel("Displacement (m)")
    ax.set_ylabel("Probability")
    ax.set_xscale("log")
    ax.set_xlim(0.01, 10)
    ax.set_ylim(0, 1)
    ax.set_xticks([0.01, 0.1, 1, 10])
    ax.set_xticklabels(["0.01", "0.1", "1", "10"])
    ax.set_yticks(np.linspace(0, 1, 11))
    ax.grid(True, which="major", alpha=0.18)
    ax.legend(frameon=False, loc="upper right")

    ax = axes[1]
    for mag in mags:
        digitized = FIG6_DIGITIZED[mag]
        style = STYLES[mag]
        ax.plot(d_curve, p_curve[mag], color=style["color"], lw=2.0, ls=style["ls"], label=style["label"])
        ax.scatter(
            digitized["displacement_m"],
            digitized["prob_exceed"],
            s=18,
            facecolor="none",
            edgecolor=style["color"],
            zorder=3,
        )
    ax.set_title("Code curves over digitized points")
    ax.set_xlabel("Displacement (m)")
    ax.set_xscale("log")
    ax.set_xlim(0.01, 10)
    ax.set_ylim(0, 1)
    ax.set_xticks([0.01, 0.1, 1, 10])
    ax.set_xticklabels(["0.01", "0.1", "1", "10"])
    ax.set_yticks(np.linspace(0, 1, 11))
    ax.grid(True, which="major", alpha=0.18)
    ax.legend(frameon=False, loc="upper right")

    fig.suptitle("Moss and Ross (2011) Figure 6: x/L=0.25, AD-normalized, gamma model", y=1.03)
    fig.savefig(outdir / "moss_ross_2011_fig6_paper_vs_code.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "moss_ross_2011_fig6_paper_vs_code.pdf", bbox_inches="tight")
    plt.close(fig)


def _save_overlay_residual(outdir, mags, d_curve, p_curve, residuals):
    fig, (ax, axr) = plt.subplots(
        2,
        1,
        figsize=(8.2, 7.2),
        sharex=True,
        gridspec_kw={"height_ratios": [3.0, 1.1]},
        constrained_layout=True,
    )

    print("Pointwise comparison at digitized displacement values:")
    for mag in mags:
        digitized = FIG6_DIGITIZED[mag]
        d = digitized["displacement_m"]
        p = digitized["prob_exceed"]
        diff = residuals[mag]
        style = STYLES[mag]

        ax.plot(d_curve, p_curve[mag], color=style["color"], ls=style["ls"], lw=2.4, label=f"code {style['label']}")
        ax.scatter(d, p, s=34, facecolor="white", edgecolor=style["color"], lw=1.4, label=f"digitized {style['label']}")
        axr.axhline(0.0, color="0.45", lw=0.8)
        axr.plot(d, diff, color=style["color"], ls=style["ls"], marker="o", ms=4, lw=1.6, label=style["label"])

        worst = int(np.argmax(np.abs(diff)))
        print(
            f"M={mag}: max_abs_diff={np.max(np.abs(diff)):.4f}, "
            f"mean_abs_diff={np.mean(np.abs(diff)):.4f}, rmse={np.sqrt(np.mean(diff**2)):.4f}"
        )
        print(
            f"  worst at d={d[worst]:.3g} m: paper={p[worst]:.4f}, "
            f"code={p[worst] + diff[worst]:.4f}, diff={diff[worst]:+.4f}"
        )

    ax.set_title("Moss and Ross (2011) Figure 6 comparison\nx/L=0.25, AD-normalized, gamma spatial model")
    ax.set_ylabel("Probability, P(D > d | m, slip)")
    ax.set_xscale("log")
    ax.set_xlim(0.01, 10)
    ax.set_ylim(0, 1)
    ax.set_yticks(np.linspace(0, 1, 11))
    ax.grid(True, which="major", alpha=0.22)
    ax.legend(frameon=False, ncols=2, loc="upper right", fontsize=10)

    axr.set_xlabel("Displacement (m)")
    axr.set_ylabel("code - digitized")
    axr.set_xscale("log")
    axr.set_xlim(0.01, 10)
    axr.set_ylim(-0.08, 0.08)
    axr.set_xticks([0.01, 0.1, 1, 10])
    axr.set_xticklabels(["0.01", "0.1", "1", "10"])
    axr.grid(True, which="major", alpha=0.22)
    axr.legend(frameon=False, ncols=3, loc="upper right", fontsize=10)

    fig.savefig(outdir / "moss_ross_2011_fig6_overlay_residual.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "moss_ross_2011_fig6_overlay_residual.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    _configure_matplotlib()
    outdir = Path("Figures")
    outdir.mkdir(exist_ok=True)
    mags, d_curve, p_curve, residuals = _compute_curves()
    _save_paper_vs_code(outdir, mags, d_curve, p_curve)
    _save_overlay_residual(outdir, mags, d_curve, p_curve, residuals)
    print(f"Wrote figures to {outdir.resolve()}")


if __name__ == "__main__":
    main()
