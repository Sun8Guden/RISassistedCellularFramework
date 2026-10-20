"""Render the README figures from the CSV files written by the validation and benchmark programs.

Usage: python3 scripts/make_figures.py [data_dir] [figure_dir]
Defaults: docs/data and docs/figures. Each figure is written in a light and a dark variant
(<name>.png and <name>-dark.png) for GitHub's prefers-color-scheme <picture> switch.
Requires numpy and matplotlib. Run `make figures` to regenerate the data first.
"""
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

THEMES = {
    "light": {
        "surface": "#fcfcfb", "ink": "#0b0b0b", "secondary": "#52514e", "muted": "#898781",
        "grid": "#e1e0d9", "axis": "#c3c2b7", "series": ["#2a78d6", "#eb6834"],
    },
    "dark": {
        "surface": "#1a1a19", "ink": "#f0efec", "secondary": "#c3c2b7", "muted": "#898781",
        "grid": "#2c2c2a", "axis": "#383835", "series": ["#3987e5", "#d95926"],
    },
}
LINE_WIDTH = 1.6    # ~2 px at the export resolution
MARKER_SIZE = 6.5   # >= 8 px diameter
RING_WIDTH = 1.4    # ~2 px ring in the surface colour


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def styled_axes(theme, title, subtitle):
    t = THEMES[theme]
    fig, ax = plt.subplots(figsize=(7.2, 4.2), dpi=160)
    fig.patch.set_facecolor(t["surface"])
    ax.set_facecolor(t["surface"])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(t["axis"])
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=t["muted"], labelcolor=t["secondary"], labelsize=9, width=0.8)
    ax.grid(True, which="major", color=t["grid"], linewidth=0.8, linestyle="-")
    ax.set_axisbelow(True)
    fig.text(0.015, 0.965, title, color=t["ink"], fontsize=12.5, fontweight="semibold", ha="left", va="top")
    fig.text(0.015, 0.905, subtitle, color=t["secondary"], fontsize=9.5, ha="left", va="top")
    fig.subplots_adjust(left=0.1, right=0.97, bottom=0.13, top=0.8)
    return fig, ax, t


def label_axes(ax, t, xlabel, ylabel):
    ax.set_xlabel(xlabel, color=t["secondary"], fontsize=9.5)
    ax.set_ylabel(ylabel, color=t["secondary"], fontsize=9.5)


def legend(ax, t, loc):
    leg = ax.legend(loc=loc, frameon=False, fontsize=9, labelcolor=t["ink"], handlelength=1.8)
    return leg


def save(fig, figure_dir, name, theme):
    suffix = "" if theme == "light" else "-dark"
    path = os.path.join(figure_dir, f"{name}{suffix}.png")
    fig.savefig(path, facecolor=fig.get_facecolor())
    plt.close(fig)
    return path


def marker_style(t, colour):
    return dict(marker="o", markersize=MARKER_SIZE, markerfacecolor=colour,
                markeredgecolor=t["surface"], markeredgewidth=RING_WIDTH, linestyle="none")


def coverage_figure(data_dir, figure_dir, theme):
    rows = read_csv(os.path.join(data_dir, "coverage_validation.csv"))
    x = [float(r["threshold_db"]) for r in rows]
    fig, ax, t = styled_axes(theme, "Analytical model vs. 4 million Monte Carlo samples",
                             "Coverage probability P(SINR > threshold) of a user 100 m from its base station")
    blue, orange = t["series"]
    ax.plot(x, [float(r["analytic_with_ris"]) for r in rows], color=orange, linewidth=LINE_WIDTH,
            solid_capstyle="round", label="With RISs: analytical")
    ax.plot(x, [float(r["mc_with_ris"]) for r in rows], label="With RISs: Monte Carlo", **marker_style(t, orange))
    ax.plot(x, [float(r["analytic_without_ris"]) for r in rows], color=blue, linewidth=LINE_WIDTH,
            solid_capstyle="round", label="Without RISs: analytical")
    ax.plot(x, [float(r["mc_without_ris"]) for r in rows], label="Without RISs: Monte Carlo", **marker_style(t, blue))
    ax.set_ylim(0, 1.02)
    ax.set_xlim(min(x) - 0.5, max(x) + 0.5)
    label_axes(ax, t, "SINR threshold [dB]", "Coverage probability")
    legend(ax, t, "upper right")
    return save(fig, figure_dir, "coverage_validation", theme)


def accuracy_figure(data_dir, figure_dir, theme):
    rows = read_csv(os.path.join(data_dir, "accuracy_vs_time.csv"))
    mc = [r for r in rows if r["method"] == "monte_carlo"]
    analytic = [r for r in rows if r["method"] == "analytic"][0]
    p = float(analytic["probability"])
    seconds_per_sample = float(mc[-1]["seconds"]) / float(mc[-1]["samples"])
    a_time, a_err = float(analytic["seconds"]), float(analytic["absolute_error"])

    fig, ax, t = styled_axes(theme, "Monte Carlo would need days for the accuracy the engine gets in 0.3 s",
                             "Coverage with RISs at 0 dB: error vs. single-thread compute time")
    blue, orange = t["series"]
    # Expected Monte Carlo error (one standard error) from 100 samples down to the analytical accuracy.
    target_samples = p * (1 - p) / a_err ** 2
    n = [100 * (target_samples / 100) ** (i / 200) for i in range(201)]
    ax.plot([k * seconds_per_sample for k in n], [(p * (1 - p) / k) ** 0.5 for k in n], color=blue,
            linewidth=LINE_WIDTH, solid_capstyle="round", label="Monte Carlo: expected error (1 standard error)")
    observed = [r for r in mc if float(r["absolute_error"]) > 0]
    ax.plot([float(r["seconds"]) for r in observed], [float(r["absolute_error"]) for r in observed],
            label="Monte Carlo: observed error", **dict(marker_style(t, blue), markersize=4.5))
    ax.plot([a_time], [a_err], label="Analytical engine (integration tolerance)", **marker_style(t, orange))
    mc_seconds = target_samples * seconds_per_sample
    ax.annotate(f"{a_time:.2f} s", (a_time, a_err), xytext=(0, -16), textcoords="offset points", ha="center",
                color=t["ink"], fontsize=9)
    ax.annotate(f"~{mc_seconds / 86400:.0f} days", (mc_seconds, a_err), xytext=(0, -16), textcoords="offset points",
                ha="center", color=t["ink"], fontsize=9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    lowest = min([a_err] + [float(r["absolute_error"]) for r in observed])
    ax.set_ylim(lowest / 3, 1)
    label_axes(ax, t, "Compute time [s]", "Absolute error in coverage probability")
    legend(ax, t, "upper right")
    return save(fig, figure_dir, "accuracy_vs_time", theme)


def cubature_figure(data_dir, figure_dir, theme):
    rows = read_csv(os.path.join(data_dir, "cubature_convergence.csv"))
    fig, ax, t = styled_axes(theme, "Adaptive cubature reaches 1e-13 where Monte Carlo is still at 1e-2",
                             "Relative error on a 2D Genz Gaussian test integral vs. number of integrand evaluations")
    blue, orange = t["series"]
    for method, colour, name in (("cubature", blue, "Adaptive cubature (degree-7 Genz-Malik)"),
                                 ("monte_carlo", orange, "Plain Monte Carlo (mean of 20 runs)")):
        sel = [r for r in rows if r["method"] == method and float(r["relative_error"]) > 0]
        ax.plot([float(r["evaluations"]) for r in sel], [float(r["relative_error"]) for r in sel], color=colour,
                linewidth=LINE_WIDTH, solid_joinstyle="round", label=name)
    ax.set_xscale("log")
    ax.set_yscale("log")
    label_axes(ax, t, "Integrand evaluations", "Relative error")
    legend(ax, t, "lower left")
    return save(fig, figure_dir, "cubature_convergence", theme)


def smile_figure(data_dir, figure_dir, theme):
    rows = read_csv(os.path.join(data_dir, "heston_smile.csv"))
    k = [float(r["strike"]) for r in rows]
    fig, ax, t = styled_axes(theme, "The same inversion engine prices options",
                             "Heston implied-volatility skew: Gil-Pelaez inversion vs. Monte Carlo (95% intervals smaller than the markers)")
    blue, orange = t["series"]
    ax.plot(k, [100 * float(r["inversion_implied_vol"]) for r in rows], color=blue, linewidth=LINE_WIDTH,
            solid_capstyle="round", label="Gil-Pelaez inversion")
    low = [100 * float(r["mc_implied_vol_low"]) for r in rows]
    high = [100 * float(r["mc_implied_vol_high"]) for r in rows]
    mid = [(a + b) / 2 for a, b in zip(low, high)]
    ax.errorbar(k, mid, yerr=[[m - a for m, a in zip(mid, low)], [b - m for m, b in zip(mid, high)]],
                fmt="none", ecolor=orange, elinewidth=LINE_WIDTH, capsize=0)
    ax.plot(k, mid, label="Monte Carlo, 800k paths", **marker_style(t, orange))
    label_axes(ax, t, "Strike (spot = 100)", "Implied volatility [%]")
    legend(ax, t, "upper right")
    return save(fig, figure_dir, "heston_smile", theme)


def main():
    data_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join("docs", "data")
    figure_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join("docs", "figures")
    os.makedirs(figure_dir, exist_ok=True)
    for theme in THEMES:
        for figure in (coverage_figure, accuracy_figure, cubature_figure, smile_figure):
            print(figure(data_dir, figure_dir, theme))


if __name__ == "__main__":
    main()
