#!/usr/bin/env python3
"""Plot pulse-area spectra for SiPM 1 and SiPM 2 from a matched scan."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEVICES = [
    ("PCB_CH3_ScopeA", "SiPM 1 (△)", "#1f77b4", "sipm1_pulse_area_spectra.png"),
    ("PCB_CH2_ScopeB", "SiPM 2 (★)", "#ff7f0e", "sipm2_pulse_area_spectra.png"),
]
LIMITS = (0.0, 12000.0)


def fd_edges(values: np.ndarray) -> np.ndarray:
    values = values[np.isfinite(values)]
    values = values[(values >= LIMITS[0]) & (values <= LIMITS[1])]
    if values.size < 2:
        return np.linspace(*LIMITS, 81)
    q25, q75 = np.percentile(values, [25, 75])
    width = 2.0 * (q75 - q25) / np.cbrt(values.size)
    if not np.isfinite(width) or width <= 0:
        width = (LIMITS[1] - LIMITS[0]) / 120.0
    bins = int(np.clip(np.ceil((LIMITS[1] - LIMITS[0]) / width), 60, 350))
    return np.linspace(LIMITS[0], LIMITS[1], bins + 1)


def as_bool(series: pd.Series) -> np.ndarray:
    if series.dtype == bool:
        return series.to_numpy(bool)
    return series.astype(str).str.lower().isin({"true", "1", "yes"}).to_numpy(bool)


def signal_directory(point_root: Path, suffix: str) -> Path:
    matches = list(point_root.glob(f"*_{suffix}"))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one signal ending in {suffix}, found {len(matches)}")
    return matches[0]


def plot_device(scan_root: Path, output_root: Path, suffix: str, label: str,
                color: str, filename: str) -> None:
    manifest = json.loads((scan_root / "scan_manifest.json").read_text())
    figure, axes = plt.subplots(2, 4, figsize=(18, 9), sharex=True)

    for axis, point in zip(axes.flat, manifest["points"]):
        directory = signal_directory(scan_root / "analysis" / point["tag"], suffix)
        summary = json.loads((directory / "point_summary.json").read_text())
        table = pd.read_csv(directory / "pulse_measurements.csv")
        selected = as_bool(table["accepted"])
        if "overflow" in table:
            selected &= ~as_bool(table["overflow"])
        values = table.loc[selected, "pulse_area_mV_ns"].to_numpy(float)

        axis.hist(values, bins=fd_edges(values), histtype="step", color=color, linewidth=1.25)
        for center in summary.get("area_selected_peak_centers_mV_ns", []):
            axis.axvline(center, color="#d62728", linewidth=1.15)
        accepted = int(summary["accepted_events"])
        total = int(summary["event_files"])
        axis.set_title(
            f"{point['effective_bias_v']:.3f} V\n"
            f"accepted {accepted:,}/{total:,} ({100.0 * accepted / total:.2f}%)"
        )
        axis.set_yscale("log")
        axis.set_xlim(*LIMITS)
        axis.grid(True, alpha=0.25)

    for axis in axes.flat[len(manifest["points"]):]:
        axis.axis("off")
    for axis in axes[:, 0]:
        axis.set_ylabel("Events per bin")
    for axis in axes[-1, :]:
        if axis.axison:
            axis.set_xlabel("Pulse area (mV ns)")

    figure.suptitle(f"{label}: pulse-area spectra", fontsize=17)
    figure.tight_layout(rect=(0, 0, 1, 0.95))
    figure.savefig(output_root / filename, dpi=180, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("scan_root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    for suffix, label, color, filename in DEVICES:
        plot_device(args.scan_root, args.output, suffix, label, color, filename)


if __name__ == "__main__":
    main()
