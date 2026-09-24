#!/usr/bin/env python3
"""Compare two SiPM breakdown-voltage fits from one matched clean scan."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DEVICES = [
    ("PCB_CH3_ScopeA", "SiPM 1 (△)", "#1f77b4", "o"),
    ("PCB_CH2_ScopeB", "SiPM 2 (★)", "#ff7f0e", "s"),
]


def resolve_signal(run_root: Path, suffix: str) -> str:
    matches = list((run_root / "results").glob(f"*_{suffix}/area/vbr_fit.json"))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one signal ending in {suffix}, found {len(matches)}")
    return matches[0].parents[1].name


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    fits = {}
    signals = {}
    manifest = json.loads((args.run_root / "scan_manifest.json").read_text())
    for suffix, label, color, marker in DEVICES:
        signal = resolve_signal(args.run_root, suffix)
        signals[suffix] = signal
        path = args.run_root / "results" / signal / "area" / "vbr_fit.json"
        fits[suffix] = json.loads(path.read_text())

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.8), gridspec_kw={"width_ratios": [1.75, 1]})
    x_line = np.linspace(50.2, 57.2, 400)
    for suffix, label, color, marker in DEVICES:
        signal = signals[suffix]
        fit = fits[suffix]
        x = np.asarray(fit["bias_points_V"], dtype=float)
        y = np.asarray(fit["pe_spacings"], dtype=float)
        yerr = np.asarray(
            [
                json.loads(
                    (
                        args.run_root
                        / "analysis"
                        / point["tag"]
                        / signal
                        / "point_summary.json"
                    ).read_text()
                )["one_pe_area_spacing_unc_mV_ns"]
                for point in manifest["points"]
            ],
            dtype=float,
        )
        axes[0].errorbar(x, y, yerr=yerr, fmt=marker, capsize=3, color=color, label=f"{label} data")
        axes[0].plot(
            x_line,
            fit["slope_per_V"] * (x_line - fit["breakdown_voltage_V"]),
            color=color,
            linewidth=1.8,
            label=(
                f"{label} fit: Vbr={fit['breakdown_voltage_V']:.3f} +/- "
                f"{fit['breakdown_voltage_unc_V']:.3f} V, R2={fit['r_squared']:.4f}"
            ),
        )
        axes[0].axvline(fit["breakdown_voltage_V"], color=color, linestyle="--", linewidth=1.2)

    axes[0].axhline(0, color="0.35", linewidth=1)
    axes[0].set(
        xlabel="Effective SiPM bias voltage (V)",
        ylabel="One-photoelectron area spacing (mV ns)",
        title="Matched pulse-area fits",
        xlim=(50.2, 57.2),
        ylim=(-150, 3600),
    )
    axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3)

    labels, values, errors, colors = [], [], [], []
    for suffix, label, color, _ in DEVICES:
        fit = fits[suffix]
        labels.append(label)
        values.append(fit["breakdown_voltage_V"])
        errors.append(fit["breakdown_voltage_unc_V"])
        colors.append(color)
    positions = np.arange(len(labels))
    axes[1].errorbar(positions, values, yerr=errors, fmt="none", ecolor=colors[0], capsize=7, linewidth=2)
    for position, value, error, color, marker in zip(positions, values, errors, colors, ["o", "s"]):
        axes[1].scatter(position, value, s=90, marker=marker, color=color, zorder=3)
        axes[1].text(position, value + error + 0.025, f"{value:.3f} +/- {error:.3f} V", ha="center", fontsize=10)
    difference = values[1] - values[0]
    difference_unc = float(np.hypot(*errors))
    significance = difference / difference_unc
    axes[1].text(
        0.5,
        min(values) - 0.22,
        f"SiPM 2 (★) - SiPM 1 (△) = {difference:.3f} +/- {difference_unc:.3f} V\n"
        f"Difference = {significance:.2f} sigma",
        ha="center",
        va="center",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "0.7"},
    )
    axes[1].set_xticks(positions, labels)
    axes[1].set(ylabel="Breakdown voltage (V)", title="Breakdown-voltage comparison", xlim=(-0.5, 1.5))
    axes[1].set_ylim(min(values) - 0.3, max(values) + 0.3)
    axes[1].grid(True, axis="y", alpha=0.3)

    fig.suptitle("Two-SiPM comparison: matched 25 mV acquisition", fontsize=17)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(args.output / "clean_two_sipm_vbr_comparison.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    summary = {
        "sipm1_vbr_V": values[0],
        "sipm1_vbr_unc_V": errors[0],
        "sipm2_vbr_V": values[1],
        "sipm2_vbr_unc_V": errors[1],
        "sipm2_minus_sipm1_V": difference,
        "difference_unc_V": difference_unc,
        "difference_significance_sigma": significance,
    }
    (args.output / "clean_two_sipm_vbr_comparison.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
