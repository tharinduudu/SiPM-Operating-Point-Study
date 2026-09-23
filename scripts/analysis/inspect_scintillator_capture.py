#!/usr/bin/env python3
"""Summarize four-channel scintillator waveforms before physics selection."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


CHANNELS = {"A": 1, "B": 2, "C": 3, "D": 4}
COLORS = {"A": "#1f77b4", "B": "#ff7f0e", "C": "#2ca02c", "D": "#d62728"}


def read_waveform(path: Path) -> np.ndarray:
    waveform = np.loadtxt(path, delimiter=",", skiprows=3, dtype=float)
    if waveform.ndim != 2 or waveform.shape[1] < 5:
        raise ValueError(f"{path} is not a four-channel waveform")
    return waveform


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--baseline-stop-ns", type=float, default=-600.0)
    parser.add_argument("--search-start-ns", type=float, default=-100.0)
    parser.add_argument("--search-stop-ns", type=float, default=200.0)
    args = parser.parse_args()

    files = sorted(args.run_dir.glob("event_[0-9]*.csv"))
    if not files:
        raise FileNotFoundError(f"No waveform CSV files in {args.run_dir}")
    args.output.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, float | int | str]] = []
    traces: dict[str, list[np.ndarray]] = {channel: [] for channel in CHANNELS}
    time_ns: np.ndarray | None = None
    for event, path in enumerate(files):
        waveform = read_waveform(path)
        time_ns = waveform[:, 0]
        baseline_mask = time_ns < args.baseline_stop_ns
        search_mask = (time_ns >= args.search_start_ns) & (time_ns <= args.search_stop_ns)
        if np.sum(baseline_mask) < 10 or not np.any(search_mask):
            raise ValueError("waveform does not cover the requested baseline and search windows")
        row: dict[str, float | int | str] = {"event": event, "file": path.name}
        for channel, column in CHANNELS.items():
            baseline = float(np.median(waveform[baseline_mask, column]))
            corrected = waveform[:, column] - baseline
            search = corrected[search_mask]
            search_time = time_ns[search_mask]
            maximum_index = int(np.argmax(search))
            minimum_index = int(np.argmin(search))
            row[f"{channel}_baseline_mV"] = baseline
            row[f"{channel}_positive_excursion_mV"] = float(search[maximum_index])
            row[f"{channel}_positive_time_ns"] = float(search_time[maximum_index])
            row[f"{channel}_negative_excursion_mV"] = float(search[minimum_index])
            row[f"{channel}_negative_time_ns"] = float(search_time[minimum_index])
            traces[channel].append(corrected)
        rows.append(row)

    table = pd.DataFrame(rows)
    table.to_csv(args.output / "event_excursions.csv", index=False)

    summaries: list[dict[str, float | int | str]] = []
    quantiles = [0.1, 0.5, 0.9, 0.99]
    for channel in CHANNELS:
        positive = table[f"{channel}_positive_excursion_mV"].to_numpy(float)
        negative = table[f"{channel}_negative_excursion_mV"].to_numpy(float)
        pos_q = np.quantile(positive, quantiles)
        neg_q = np.quantile(negative, quantiles)
        summaries.append(
            {
                "channel": channel,
                "events": len(table),
                "baseline_median_mV": float(np.median(table[f"{channel}_baseline_mV"])),
                **{f"positive_q{int(q * 100):02d}_mV": float(v) for q, v in zip(quantiles, pos_q)},
                **{f"negative_q{int(q * 100):02d}_mV": float(v) for q, v in zip(quantiles, neg_q)},
                "positive_time_median_ns": float(np.median(table[f"{channel}_positive_time_ns"])),
                "negative_time_median_ns": float(np.median(table[f"{channel}_negative_time_ns"])),
            }
        )
    summary = pd.DataFrame(summaries)
    summary.to_csv(args.output / "channel_summary.csv", index=False)

    assert time_ns is not None
    figure, axes = plt.subplots(2, 2, figsize=(13.5, 8.5), sharex=True)
    for axis, channel in zip(axes.flat, CHANNELS):
        values = np.asarray(traces[channel], dtype=float)
        low, median, high = np.percentile(values, [10, 50, 90], axis=0)
        for trace in values[:: max(1, len(values) // 20)]:
            axis.plot(time_ns, trace, color=COLORS[channel], alpha=0.08, linewidth=0.7)
        axis.fill_between(time_ns, low, high, color=COLORS[channel], alpha=0.18, label="10-90% band")
        axis.plot(time_ns, median, color=COLORS[channel], linewidth=1.5, label="median")
        axis.axvspan(args.search_start_ns, args.search_stop_ns, color="0.4", alpha=0.08)
        axis.axhline(0.0, color="0.25", linewidth=0.8)
        axis.set_title(f"Scope {channel}")
        axis.set_ylabel("Baseline-subtracted voltage (mV)")
        axis.grid(alpha=0.22)
        axis.legend()
    axes[1, 0].set_xlabel("Time relative to trigger (ns)")
    axes[1, 1].set_xlabel("Time relative to trigger (ns)")
    figure.suptitle("Four-channel trigger pilot: waveform timing and polarity")
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    figure.savefig(args.output / "waveform_timing_polarity.png", dpi=220)
    plt.close(figure)

    figure, axes = plt.subplots(2, 2, figsize=(13.5, 8.5))
    for axis, channel in zip(axes.flat, CHANNELS):
        positive = table[f"{channel}_positive_excursion_mV"].to_numpy(float)
        negative_magnitude = -table[f"{channel}_negative_excursion_mV"].to_numpy(float)
        axis.hist(positive, bins=30, histtype="step", linewidth=1.5, color=COLORS[channel], label="positive")
        axis.hist(negative_magnitude, bins=30, histtype="step", linewidth=1.5, color="0.35", label="negative magnitude")
        axis.set_title(f"Scope {channel}")
        axis.set_xlabel("Largest excursion in search window (mV)")
        axis.set_ylabel("Events per bin")
        axis.grid(alpha=0.22)
        axis.legend()
    figure.suptitle("Four-channel trigger pilot: pulse-excursion distributions")
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    figure.savefig(args.output / "pulse_excursion_histograms.png", dpi=220)
    plt.close(figure)

    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
