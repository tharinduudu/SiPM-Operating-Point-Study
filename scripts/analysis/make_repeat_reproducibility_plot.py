#!/usr/bin/env python3
"""Compare two identical SiPM breakdown-voltage scans."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


SIGNALS = [
    ("Triangle_PCB_CH3_ScopeA", "Triangle SiPM", "#1f77b4"),
    ("Star_PCB_CH2_ScopeB", "Star SiPM", "#ff7f0e"),
]


def load_fit(root: Path, signal: str) -> dict:
    path = root / "results" / signal / "area" / "vbr_fit.json"
    with path.open() as handle:
        return json.load(handle)


def fit_line(fit: dict, x: np.ndarray) -> np.ndarray:
    return fit["slope_per_V"] * x + fit["intercept"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("first_run", type=Path)
    parser.add_argument("repeat_run", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    fits = {
        signal: (load_fit(args.first_run, signal), load_fit(args.repeat_run, signal))
        for signal, _, _ in SIGNALS
    }

    plt.style.use("default")
    fig = plt.figure(figsize=(15, 9), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[2.2, 1.25])

    summary = {}
    for column, (signal, label, color) in enumerate(SIGNALS):
        ax = fig.add_subplot(grid[0, column])
        first, repeat = fits[signal]
        for run_label, fit, marker, filled in [
            ("First run", first, "o", True),
            ("Sanity-check repeat", repeat, "s", False),
        ]:
            x = np.asarray(fit["bias_points_V"], dtype=float)
            y = np.asarray(fit["pe_spacings"], dtype=float)
            face = color if filled else "white"
            ax.scatter(
                x,
                y,
                s=58,
                marker=marker,
                facecolors=face,
                edgecolors=color,
                linewidths=1.8,
                zorder=3,
                label=(
                    f"{run_label}: Vbr = {fit['breakdown_voltage_V']:.3f} "
                    f"+/- {fit['breakdown_voltage_unc_V']:.3f} V"
                ),
            )
            line_x = np.linspace(min(x) - 0.05, max(x) + 0.05, 200)
            ax.plot(
                line_x,
                fit_line(fit, line_x),
                color=color,
                linewidth=2,
                linestyle="-" if filled else "--",
                alpha=0.9,
            )

        delta = repeat["breakdown_voltage_V"] - first["breakdown_voltage_V"]
        delta_unc = math.hypot(
            first["breakdown_voltage_unc_V"], repeat["breakdown_voltage_unc_V"]
        )
        z_score = delta / delta_unc
        summary[signal] = {
            "label": label,
            "first_vbr_V": first["breakdown_voltage_V"],
            "first_vbr_unc_V": first["breakdown_voltage_unc_V"],
            "repeat_vbr_V": repeat["breakdown_voltage_V"],
            "repeat_vbr_unc_V": repeat["breakdown_voltage_unc_V"],
            "repeat_minus_first_V": delta,
            "difference_unc_V": delta_unc,
            "difference_sigma": z_score,
            "first_r_squared": first["r_squared"],
            "repeat_r_squared": repeat["r_squared"],
        }
        ax.set_title(label, fontsize=15, weight="bold")
        ax.set_xlabel("Effective bias voltage (V)", fontsize=12)
        ax.set_ylabel("1 p.e. pulse-area spacing (mV ns)", fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=10, loc="upper left")
        ax.text(
            0.98,
            0.05,
            f"Repeat - first = {delta:+.3f} +/- {delta_unc:.3f} V\n"
            f"Difference = {abs(z_score):.2f} sigma",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=11,
            bbox={"facecolor": "white", "edgecolor": "0.7", "alpha": 0.9},
        )

    ax = fig.add_subplot(grid[1, :])
    positions = np.arange(len(SIGNALS), dtype=float)
    offsets = [-0.11, 0.11]
    for run_index, (run_label, marker, filled) in enumerate(
        [("First run", "o", True), ("Sanity-check repeat", "s", False)]
    ):
        for index, (signal, _, color) in enumerate(SIGNALS):
            fit = fits[signal][run_index]
            ax.errorbar(
                positions[index] + offsets[run_index],
                fit["breakdown_voltage_V"],
                yerr=fit["breakdown_voltage_unc_V"],
                fmt=marker,
                markersize=9,
                markerfacecolor=color if filled else "white",
                markeredgecolor=color,
                markeredgewidth=1.8,
                ecolor=color,
                elinewidth=1.8,
                capsize=5,
                label=run_label if index == 0 else None,
            )
    ax.set_xticks(positions, [label for _, label, _ in SIGNALS])
    ax.set_ylabel("Fitted breakdown voltage (V)", fontsize=12)
    ax.set_title("Breakdown-voltage reproducibility", fontsize=14, weight="bold")
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(fontsize=11)

    fig.suptitle(
        "Repeated SiPM breakdown-voltage measurement",
        fontsize=19,
        weight="bold",
    )
    png = args.output / "vbr_repeat_reproducibility.png"
    fig.savefig(png, dpi=180)
    plt.close(fig)

    with (args.output / "vbr_repeat_reproducibility.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
    print(png)


if __name__ == "__main__":
    main()
