#!/usr/bin/env python3
"""Summarize representative raw waveforms without copying the full raw dataset."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyze_point import read_waveform, robust_baseline


COLORS = {"SIPM1": "#1f77b4", "SIPM2": "#ff7f0e"}
SIGNALS = {"SIPM1": "B", "SIPM2": "A"}
TAGS = ("bias_53p400", "bias_55p500", "bias_56p900")


def load_ensemble(root: Path, tag: str, sipm: str, maximum_events: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    analysis = root / "analysis_aligned" / tag / sipm
    summary = json.loads((analysis / "point_summary.json").read_text(encoding="utf-8"))
    table = pd.read_csv(analysis / "pulse_measurements.csv")
    table = table[table["accepted"].fillna(False)].copy()
    selected_centers = summary.get("area_selected_peak_centers_mV_ns", [])
    gap_value = summary.get("one_pe_area_spacing_mV_ns")
    if selected_centers and isinstance(gap_value, (int, float)) and np.isfinite(gap_value):
        first_center = float(selected_centers[0])
        gap = float(gap_value)
    else:
        area_values = table["pulse_area_mV_ns"].to_numpy(float)
        low, high = np.percentile(area_values, [0.5, 99.5])
        counts, edges = np.histogram(area_values, bins="fd", range=(low, high))
        first_center = float(0.5 * (edges[np.argmax(counts)] + edges[np.argmax(counts) + 1]))
        gap = max(0.5 * abs(first_center), 500.0)
    selected = table[np.abs(table["pulse_area_mV_ns"] - first_center) <= 0.24 * gap]
    if len(selected) < 40:
        selected = table.iloc[np.argsort(np.abs(table["pulse_area_mV_ns"] - first_center))[:maximum_events]]
    if len(selected) > maximum_events:
        selected = selected.iloc[np.linspace(0, len(selected) - 1, maximum_events, dtype=int)]

    grid = np.arange(-60.0, 301.0, 1.0)
    aligned = []
    signal_column = 1 if SIGNALS[sipm] == "A" else 2
    for row in selected.itertuples(index=False):
        waveform = read_waveform(root / "raw" / tag / sipm / row.file)
        time_ns = waveform[:, 0]
        signal_mv = waveform[:, signal_column]
        baseline, _ = robust_baseline(signal_mv[time_ns < -30.0])
        corrected = signal_mv - baseline
        aligned.append(np.interp(grid, time_ns - float(row.peak_time_ns), corrected))
    traces = np.asarray(aligned, float)
    return grid, traces, selected["pulse_area_mV_ns"].to_numpy(float), float(summary["bias_v"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--events", type=int, default=240)
    args = parser.parse_args()
    root = args.run_root.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({"font.size": 10.5, "axes.grid": True, "grid.alpha": 0.22, "figure.dpi": 170, "savefig.dpi": 220, "savefig.bbox": "tight"})
    ensembles: dict[tuple[str, str], tuple[np.ndarray, np.ndarray, np.ndarray, float]] = {}
    summary_rows = []
    for sipm in ("SIPM1", "SIPM2"):
        for tag in TAGS:
            grid, traces, areas, bias = load_ensemble(root, tag, sipm, args.events)
            ensembles[(sipm, tag)] = (grid, traces, areas, bias)
            median = np.median(traces, axis=0)
            full_region = (grid >= -20.0) & (grid <= 300.0)
            gate_region = (grid >= -20.0) & (grid <= 200.0)
            full_charge = float(np.trapezoid(median[full_region], grid[full_region]))
            gate_charge = float(np.trapezoid(median[gate_region], grid[gate_region]))
            summary_rows.append({"sipm": sipm, "tag": tag, "bias_v": bias, "waveforms": len(traces), "median_full_charge_mV_ns": full_charge, "median_gate_charge_mV_ns": gate_charge, "gate_fraction": gate_charge / full_charge if full_charge else np.nan})

    figure, axes = plt.subplots(2, 3, figsize=(15.0, 8.2), sharex=True)
    for row, sipm in enumerate(("SIPM1", "SIPM2")):
        for column, tag in enumerate(TAGS):
            axis = axes[row, column]
            grid, traces, _, bias = ensembles[(sipm, tag)]
            median = np.median(traces, axis=0)
            lower, upper = np.percentile(traces, [16, 84], axis=0)
            axis.fill_between(grid, lower, upper, color=COLORS[sipm], alpha=0.22, label="16-84% envelope")
            axis.plot(grid, median, color=COLORS[sipm], linewidth=1.8, label="median selected population")
            axis.axvspan(-20.0, 200.0, color="0.5", alpha=0.08, label="charge gate" if row == 0 and column == 0 else None)
            axis.axvline(0, color="0.35", linestyle=":", linewidth=1.0)
            axis.set_title(f"{sipm}, {bias:.3f} V")
            axis.set_xlabel("Time relative to pulse maximum (ns)")
            axis.set_ylabel("Amplifier output (mV)")
            if row == 0 and column == 0:
                axis.legend()
    figure.suptitle("First resolved photoelectron-population waveform shape", fontsize=15)
    figure.tight_layout()
    figure.savefig(output / "single_pe_waveform_shapes.png")
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(12.5, 5.0), sharey=True)
    for axis, sipm in zip(axes, ("SIPM1", "SIPM2")):
        for tag in TAGS:
            grid, traces, _, bias = ensembles[(sipm, tag)]
            median = np.median(traces, axis=0)
            region = (grid >= -20.0) & (grid <= 300.0)
            time = grid[region]
            values = median[region]
            cumulative = np.concatenate([[0.0], np.cumsum(0.5 * (values[:-1] + values[1:]) * np.diff(time))])
            final = cumulative[-1]
            axis.plot(time, 100.0 * cumulative / final, linewidth=1.7, label=f"{bias:.3f} V")
        axis.axvline(200.0, color="black", linestyle="--", linewidth=1.1, label="gate end")
        axis.axhline(100.0, color="0.45", linewidth=0.8)
        axis.set_title(sipm)
        axis.set_xlabel("Time after pulse maximum (ns)")
        axis.set_ylabel("Cumulative charge (% of -20 to +300 ns)")
        axis.legend()
    figure.suptitle("Charge-integration convergence", fontsize=15)
    figure.tight_layout()
    figure.savefig(output / "charge_gate_convergence.png")
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(12.5, 5.0), sharex=True)
    for axis, sipm in zip(axes, ("SIPM1", "SIPM2")):
        grid, traces, _, bias = ensembles[(sipm, "bias_55p500")]
        for trace in traces[:120]:
            axis.plot(grid, trace, color=COLORS[sipm], alpha=0.035, linewidth=0.7)
        axis.plot(grid, np.median(traces, axis=0), color="black", linewidth=1.8, label="median")
        axis.axvspan(-20.0, 200.0, color="0.5", alpha=0.08)
        axis.set_title(f"{sipm}, {bias:.3f} V")
        axis.set_xlabel("Time relative to pulse maximum (ns)")
        axis.set_ylabel("Amplifier output (mV)")
        axis.legend()
    figure.suptitle("Individual waveforms in the first resolved p.e. population", fontsize=15)
    figure.tight_layout()
    figure.savefig(output / "waveform_overlay_55p5V.png")
    plt.close(figure)

    pd.DataFrame(summary_rows).to_csv(output / "waveform_gate_summary.csv", index=False)
    print(pd.DataFrame(summary_rows).to_string(index=False))


if __name__ == "__main__":
    main()
