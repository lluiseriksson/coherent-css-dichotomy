"""Generate publication-quality figures for Paper 07 from exact JSON results."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RES_DIR = ROOT / "results"
FIG_DIR = ROOT / "paper" / "figures"


def style_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, linestyle="--", alpha=0.35, linewidth=0.7)


def make_fig1(exp2: dict) -> None:
    """Fig 1: Universal Even-Odd Parity Dichotomy (Odd-|L_Z| vs Even-|L_Z| curves)."""
    sweeps = exp2["sweeps"]
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.1), sharey=True)

    # Panel (a): Odd-|L_Z| codes (F_CPTP == 1 identically!)
    ax = axes[0]
    style_axes(ax)
    odd_specs = [
        ("shor_9_1_3", "Shor [[9,1,3]] (Pauli)", "#1f77b4", "-", "o"),
        ("steane_7_1_3", "Steane [[7,1,3]] (Pauli)", "#2ca02c", "-", "s"),
        ("qrm_15_1_3", "QRM [[15,1,3]] (Pauli)", "#d62728", "-", "^"),
        ("rotsurf_d3", "RotSurf d=3 [[9,1,3]] (Pauli)", "#9467bd", "-", "D"),
        ("swap_immune_8_1_2_odd", "Swap-Immune [[8,1,2]] (Pauli)", "#ff7f0e", "--", "v"),
    ]
    xs = [pt["theta_over_pi"] for pt in sweeps["steane_7_1_3"]["points"]]
    ax.plot(xs, [1.0] * len(xs), color="#111111", linewidth=2.6, linestyle="-",
            label=r"Optimal CPTP ($F_e^{\mathrm{CPTP}} \equiv 1$, all odd $|L_Z|$)", zorder=6)
    for key, lbl, col, ls, mk in odd_specs:
        pts = sweeps[key]["points"]
        ax.plot(
            [p["theta_over_pi"] for p in pts],
            [p["F_Pauli"] for p in pts],
            label=lbl,
            color=col,
            linestyle=ls,
            marker=mk,
            markersize=4.0,
            markevery=2,
            linewidth=1.8,
        )
    ax.axvline(0.25, color="#555555", linestyle=":", linewidth=1.1)
    ax.axhline(0.5, color="#888888", linestyle=":", linewidth=0.9)
    ax.annotate(
        r"Swap collapse $F_e^{\mathrm{Pauli}}(\pi/4)=\frac{1}{2}$",
        xy=(0.249, 0.502),
        xytext=(0.015, 0.528),
        fontsize=8.3,
        arrowprops=dict(arrowstyle="->", color="#444444", lw=0.95),
    )
    ax.set_title(r"(a) Odd $|L_Z|$: Unitary Recovery ($F_e^{\mathrm{CPTP}} \equiv 1$)", fontsize=10.5, fontweight="bold")
    ax.set_xlabel(r"Coherent rotation angle $\theta / \pi$", fontsize=10)
    ax.set_ylabel(r"Entanglement fidelity $F_e(\theta)$", fontsize=10)
    ax.set_xlim(-0.005, 0.262)
    ax.set_ylim(0.47, 1.035)
    # Place legend cleanly BELOW the axes so zero curves are obscured!
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=7.6, frameon=True, framealpha=0.95)

    # Panel (b): Even-|L_Z| codes (F_CPTP == F_Pauli identically!) — highlighting 90° orientation split!
    ax = axes[1]
    style_axes(ax)
    even_specs = [
        ("rotsurf_d2_Xmaj", r"RotSurf $d=2$ X-maj ($F_e^{\mathrm{CPTP}}=F_e^{\mathrm{Pauli}}=F_e^{\mathrm{twirl}}$)", "#1f77b4", "-", "o"),
        ("rotsurf_d2_Zmaj", r"RotSurf $d=2$ Z-maj (90$^\circ$ rot., swap-immune)", "#1f77b4", "--", "^"),
        ("rotsurf_d4_Xmaj", r"RotSurf $d=4$ X-maj ($F_e^{\mathrm{CPTP}}=F_e^{\mathrm{Pauli}}\to\frac{1}{2}$)", "#2ca02c", "-", "s"),
        ("rotsurf_d4_Zmaj", r"RotSurf $d=4$ Z-maj (90$^\circ$ rot., $F_e(\pi/4)=1$)", "#d62728", "--", "v"),
    ]
    for key, lbl, col, ls, mk in even_specs:
        pts = sweeps[key]["points"]
        ax.plot(
            [p["theta_over_pi"] for p in pts],
            [p["F_Pauli"] for p in pts],
            label=lbl,
            color=col,
            linestyle=ls,
            marker=mk,
            markersize=4.3,
            markevery=2,
            linewidth=1.9,
        )
    pts_d4 = sweeps["rotsurf_d4_Xmaj"]["points"]
    ax.plot(
        [p["theta_over_pi"] for p in pts_d4],
        [p["F_Twirl"] for p in pts_d4],
        label=r"RotSurf $d=4$ X-maj Pauli-twirled ($F_e^{\mathrm{twirl}}$)",
        color="#2ca02c",
        linestyle=":",
        linewidth=1.5,
        alpha=0.85,
    )
    ax.axvline(0.25, color="#555555", linestyle=":", linewidth=1.1)
    ax.axhline(0.5, color="#888888", linestyle=":", linewidth=0.9)
    ax.annotate(
        r"90$^\circ$ Z-maj witness $g^*$: $F_e(\pi/4)=1$",
        xy=(0.249, 0.998),
        xytext=(0.075, 0.945),
        fontsize=8.3,
        arrowprops=dict(arrowstyle="->", color="#b22222", lw=1.0),
    )
    ax.set_title(r"(b) Even $|L_Z|$: Pauli Optimality & 90$^\circ$ Orientation Split", fontsize=10.3, fontweight="bold")
    ax.set_xlabel(r"Coherent rotation angle $\theta / \pi$", fontsize=10)
    ax.set_xlim(-0.005, 0.262)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=1, fontsize=7.5, frameon=True, framealpha=0.95)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig1_parity_dichotomy_curves.png", dpi=260, bbox_inches="tight")
    plt.close(fig)


def make_fig2(exp2: dict) -> None:
    """Fig 2: Coherent vs Pauli-twirled gap (Shor gain & exact p=1/4 zero, d=2 cancellation, QRM penalty)."""
    sweeps = exp2["sweeps"]
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.0))

    specs = [
        ("shor_9_1_3", "Shor [[9,1,3]]", "#1f77b4", "-", "o"),
        ("steane_7_1_3", "Steane [[7,1,3]]", "#2ca02c", "-", "s"),
        ("rotsurf_d4_Xmaj", "RotSurf d=4 X-maj [[16,1,4]]", "#9467bd", "-", "D"),
        ("rotsurf_d3", "RotSurf d=3 [[9,1,3]]", "#8c564b", "--", "^"),
        ("rotsurf_d2_Xmaj", "RotSurf d=2 X-maj [[4,1,2]] (exact 0)", "#111111", ":", "x"),
        ("qrm_15_1_3", "QRM [[15,1,3]] (coherent penalty)", "#d62728", "-", "v"),
    ]

    # Panel (a): Delta = F_Pauli - F_Twirl (unclipped full range [-0.37, 0.49])
    ax = axes[0]
    style_axes(ax)
    ax.axhline(0.0, color="#333333", linewidth=1.0, linestyle="-", alpha=0.6)
    for key, lbl, col, ls, mk in specs:
        pts = sweeps[key]["points"]
        ax.plot(
            [p["p"] for p in pts],
            [p["delta_coh_minus_twirl"] for p in pts],
            label=lbl,
            color=col,
            linestyle=ls,
            marker=mk,
            markersize=4.0,
            markevery=2,
            linewidth=1.8,
        )
    ax.fill_between([0, 0.5], -0.38, 0.0, color="#d62728", alpha=0.05)
    ax.text(0.22, -0.25, "Coherence hurts Pauli decoding", fontsize=8.4, color="#a01818", style="italic")
    ax.set_title(r"(a) Exact Coherent-Minus-Twirled Gap $F_e^{\mathrm{Pauli}}(p) - F_e^{\mathrm{twirl}}(p)$", fontsize=10.1, fontweight="bold")
    ax.set_xlabel(r"Effective phase-flip parameter $p = \sin^2\theta$", fontsize=10)
    ax.set_ylabel(r"$\Delta_{\mathrm{coh}}(p) = F_e^{\mathrm{Pauli}}(p) - F_e^{\mathrm{twirl}}(p)$", fontsize=10)
    ax.set_xlim(-0.01, 0.51)
    ax.set_ylim(-0.37, 0.49)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=7.6, frameon=True, framealpha=0.95)

    # Panel (b): Infidelity ratio (1 - F_Pauli) / (1 - F_Twirl) on log scale for p in (0, 0.5)
    ax = axes[1]
    style_axes(ax)
    ax.axhline(1.0, color="#333333", linewidth=1.1, linestyle="--", alpha=0.7)
    for key, lbl, col, ls, mk in specs:
        pts = sweeps[key]["points"]
        xs, ys = [], []
        for p in pts:
            if 0.0 < p["p"] < 0.5 and (1.0 - p["F_Twirl"]) > 1e-12:
                ratio = (1.0 - p["F_Pauli"]) / (1.0 - p["F_Twirl"])
                if ratio > 0.0:
                    xs.append(p["p"])
                    ys.append(ratio)
                else:
                    # Exact zero infidelity (Shor at p=1/4): break line cleanly
                    xs.append(p["p"])
                    ys.append(float("nan"))
        ax.semilogy(
            xs,
            ys,
            label=lbl,
            color=col,
            linestyle=ls,
            marker=mk,
            markersize=4.0,
            markevery=2,
            linewidth=1.8,
        )
    ax.axvline(0.25, color="#1f77b4", linestyle=":", linewidth=1.0, alpha=0.7)
    ax.annotate(
        r"Shor exact zero infidelity" + "\n" + r"$F_e^{\mathrm{Pauli}}(1/4) = 1$ ($\theta=\pi/6$)",
        xy=(0.25, 0.00018),
        xytext=(0.015, 0.0006),
        fontsize=8.1,
        color="#1f4e78",
        arrowprops=dict(arrowstyle="->", color="#1f4e78", lw=0.95),
    )
    ax.set_title(r"(b) Coherent-to-Twirled Infidelity Ratio $(1-F_e^{\mathrm{Pauli}})/(1-F_e^{\mathrm{twirl}})$", fontsize=10.1, fontweight="bold")
    ax.set_xlabel(r"Effective phase-flip parameter $p = \sin^2\theta$", fontsize=10)
    ax.set_ylabel(r"Infidelity ratio (lower = stronger coherent gain)", fontsize=10)
    ax.set_xlim(0.0, 0.50)
    ax.set_ylim(1.2e-4, 12.0)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=7.6, frameon=True, framealpha=0.95)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig2_coherent_vs_twirled_gap.png", dpi=260, bbox_inches="tight")
    plt.close(fig)


def make_fig3(exp3: dict, exp4: dict) -> None:
    """Fig 3: Theorem 2.2 0-1 Swap-Point Census + Surface/Repetition Parity Dichotomy."""
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.4))

    # Panel (a): Exact Theorem 2.2 0-1 Binary Law Census (129 codes, 4 quadrants)
    ax = axes[0]
    style_axes(ax)
    qc = exp3["quadrant_counts"]
    labels = [
        "Odd $|L_Z|$\n$L_Z \\perp \\mathcal{K}_Z$\n($F_P=\\frac{1}{2}, F_C=1$)",
        "Even $|L_Z|$\n$L_Z \\perp \\mathcal{K}_Z$\n($F_P=F_C=\\frac{1}{2}$)",
        "Odd $|L_Z|$\n$L_Z \\cdot \\mathcal{K}_Z \\neq 0$\n($F_P=F_C=1$)",
        "Even $|L_Z|$\n$L_Z \\cdot \\mathcal{K}_Z \\neq 0$\n($F_P=F_C=1$)",
    ]
    counts = [
        qc["odd_lz_collapse_half"],
        qc["even_lz_collapse_half"],
        qc["odd_lz_immune_unit"],
        qc["even_lz_immune_unit"],
    ]
    colors = ["#1f77b4", "#9467bd", "#ff7f0e", "#2ca02c"]
    bars = ax.bar(range(4), counts, color=colors, width=0.58, edgecolor="#222222", linewidth=0.9)
    for bar, cnt in zip(bars, counts):
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            cnt + 1.2,
            f"{cnt}",
            ha="center",
            va="bottom",
            fontsize=9.5,
            fontweight="bold",
        )
    ax.set_xticks(range(4))
    ax.set_xticklabels(labels, fontsize=8.2)
    ax.set_ylabel("Number of certified CSS codes", fontsize=10)
    ax.set_ylim(0, max(counts) + 12)
    ax.set_title(
        f"(a) Theorem 2.2 Binary Census at $\\theta=\\pi/4$ ({exp3['total_codes_tested']} codes, 0 intermediate)",
        fontsize=10.0,
        fontweight="bold",
    )

    # Panel (b): Even-Odd Distance Parity Dichotomy across RotSurf (d=2..4) and Rep (d=2..9) at p=1/5
    ax = axes[1]
    style_axes(ax)
    rep = exp4["phase_repetition"]
    surf = [r for r in exp4["rotated_surface"] if r["orientation"] == "X_majority"]
    surf_z = [r for r in exp4["rotated_surface"] if r["orientation"] == "Z_majority"]
    rep_d = [r["d"] for r in rep]
    rep_pauli = [r["evaluations"]["1/5"]["F_Pauli"] for r in rep]
    rep_cptp = [r["evaluations"]["1/5"]["F_CPTP"] for r in rep]
    surf_d = [r["d"] for r in surf]
    surf_pauli = [r["evaluations"]["1/5"]["F_Pauli"] for r in surf]
    surf_cptp = [r["evaluations"]["1/5"]["F_CPTP"] for r in surf]
    surfz_d = [r["d"] for r in surf_z]
    surfz_pauli = [r["evaluations"]["1/5"]["F_Pauli"] for r in surf_z]

    ax.plot(rep_d, rep_cptp, color="#1f77b4", marker="o", linewidth=1.9, markersize=6,
            label=r"Phase-Repetition $F_e^{\mathrm{CPTP}}(p=1/5)$ (sawtooth)")
    ax.plot(rep_d, rep_pauli, color="#1f77b4", marker="o", linestyle="--", linewidth=1.4, markersize=4.5, alpha=0.75,
            label=r"Phase-Repetition $F_e^{\mathrm{Pauli}}(p=1/5)$")
    ax.plot(surf_d, surf_cptp, color="#d62728", marker="s", linewidth=2.2, markersize=6.8,
            label=r"RotSurf X-maj $F_e^{\mathrm{CPTP}}(p=1/5)$ ($d=2,3,4$)")
    ax.plot(surf_d, surf_pauli, color="#d62728", marker="s", linestyle=":", linewidth=1.6, markersize=5.2, alpha=0.85,
            label=r"RotSurf X-maj $F_e^{\mathrm{Pauli}}(p=1/5)$ ($d=2,3,4$)")
    ax.scatter(surfz_d, surfz_pauli, color="#ff7f0e", marker="^", s=58, zorder=6,
               label=r"RotSurf Z-maj $F_e^{\mathrm{CPTP}}=F_e^{\mathrm{Pauli}}$ ($d=2,4$)")

    ax.set_xticks(range(2, 10))
    ax.set_xlabel(r"Code distance $d$ (odd $d \Rightarrow F_e^{\mathrm{CPTP}}=1$; even $d \Rightarrow F_e^{\mathrm{CPTP}}=F_e^{\mathrm{Pauli}}$)", fontsize=9.3)
    ax.set_ylabel(r"Exact entanglement fidelity at $p = 1/5$", fontsize=10)
    ax.set_ylim(0.615, 1.035)
    ax.set_title(r"(b) Even-Odd Distance Sawtooth Dichotomy ($p = 1/5$)", fontsize=10.0, fontweight="bold")
    ax.legend(loc="lower right", fontsize=7.6, frameon=True, framealpha=0.94)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig3_swap_point_binary_law.png", dpi=260, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    exp2 = json.loads((RES_DIR / "exp2_angle_sweeps.json").read_text(encoding="utf-8"))
    exp3 = json.loads((RES_DIR / "exp3_swap_point_binary_law.json").read_text(encoding="utf-8"))
    exp4 = json.loads((RES_DIR / "exp4_surface_and_rep_parity.json").read_text(encoding="utf-8"))

    make_fig1(exp2)
    make_fig2(exp2)
    make_fig3(exp3, exp4)
    print(f"Generated 3 figures in {FIG_DIR}")


if __name__ == "__main__":
    main()
