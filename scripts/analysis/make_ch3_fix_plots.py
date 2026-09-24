#!/usr/bin/env python3
"""Plot the CH3 trigger-threshold diagnosis and corrected bias scan."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COLOR = "#1f77b4"
RED = "#d62728"


def fd_edges(values: np.ndarray, limits: tuple[float, float]) -> np.ndarray:
    values = values[np.isfinite(values)]
    values = values[(values >= limits[0]) & (values <= limits[1])]
    if values.size < 2:
        return np.linspace(*limits, 81)
    q25, q75 = np.percentile(values, [25, 75])
    width = 2.0 * (q75 - q25) / np.cbrt(values.size)
    if not np.isfinite(width) or width <= 0:
        width = (limits[1] - limits[0]) / 120.0
    bins = int(np.clip(np.ceil((limits[1] - limits[0]) / width), 60, 350))
    return np.linspace(limits[0], limits[1], bins + 1)


def threshold_summary(threshold_root: Path, output_root: Path) -> None:
    rows = []
    for threshold in [10, 15, 20, 25, 30, 35, 40, 45]:
        directory = threshold_root / "analysis" / f"trigger_{threshold}mV"
        summary = json.loads((directory / "point_summary.json").read_text())
        data = pd.read_csv(directory / "pulse_measurements.csv")
        rows.append((threshold, summary, data))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2))
    thresholds = np.array([row[0] for row in rows])
    accepted = np.array([row[1]["accepted_events"] for row in rows])
    totals = np.array([row[1]["event_files"] for row in rows])
    fraction = accepted / totals
    error = 100.0 * np.sqrt(fraction * (1.0 - fraction) / totals)
    axes[0].errorbar(thresholds, 100.0 * fraction, yerr=error, marker="o", capsize=3, color=COLOR)
    axes[0].axvline(25, color=RED, linestyle="--", label="Selected: 25 mV")
    axes[0].set(xlabel="Trigger threshold (mV)", ylabel="Accepted events (%)", ylim=(-3, 104), title="CH3 pulse acceptance")
    axes[0].legend()

    area = np.array([row[1]["one_pe_area_spacing_mV_ns"] for row in rows], dtype=float)
    area_unc = np.array([row[1]["one_pe_area_spacing_unc_mV_ns"] for row in rows], dtype=float)
    valid = np.isfinite(area) & np.isfinite(area_unc)
    axes[1].errorbar(thresholds[valid], area[valid], yerr=area_unc[valid], marker="o", capsize=3, color=COLOR)
    axes[1].axvline(25, color=RED, linestyle="--", label="Selected: 25 mV")
    axes[1].set(xlabel="Trigger threshold (mV)", ylabel="1 p.e. area spacing (mV ns)", title="Recovered photoelectron spacing")
    axes[1].legend()
    for ax in axes:
        ax.grid(True, alpha=0.3)
    fig.suptitle("PCB CH3 trigger-threshold test at 56.401 V", fontsize=16)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(output_root / "ch3_threshold_selection.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(2, 4, figsize=(18, 9), sharex=True)
    for ax, (threshold, summary, data) in zip(axes.flat, rows):
        values = data["sampled_peak_height_mV"].to_numpy(float)
        edges = fd_edges(values, (0.0, 250.0))
        ax.hist(values, bins=edges, histtype="step", color=COLOR, linewidth=1.25)
        ax.axvline(threshold, color=RED, linewidth=1.2)
        count, total = summary["accepted_events"], summary["event_files"]
        ax.set_title(f"Trigger {threshold} mV\naccepted {count:,}/{total:,} ({100*count/total:.2f}%)")
        ax.set_yscale("log")
        ax.set_xlim(0, 250)
        ax.grid(True, alpha=0.25)
    for ax in axes[:, 0]:
        ax.set_ylabel("All triggered events per bin")
    for ax in axes[-1, :]:
        ax.set_xlabel("Sampled peak height (mV)")
    fig.suptitle("PCB CH3 triggered-event spectra", fontsize=17)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(output_root / "ch3_threshold_histogram_panels.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    pd.DataFrame(
        {
            "trigger_threshold_mV": thresholds,
            "accepted_events": accepted,
            "total_events": totals,
            "acceptance_percent": 100.0 * fraction,
            "area_spacing_mV_ns": area,
            "area_spacing_unc_mV_ns": area_unc,
        }
    ).to_csv(output_root / "ch3_threshold_summary.csv", index=False)


def corrected_scan_panel(scan_root: Path, output_root: Path, metric: str) -> None:
    settings = {
        "peak_height_mV": ("Corrected CH3 pulse-height spectra", "Pulse height (mV)", (0.0, 250.0), "selected_peak_centers_mV"),
        "pulse_area_mV_ns": ("Corrected CH3 pulse-area spectra", "Pulse area (mV ns)", (0.0, 12000.0), "area_selected_peak_centers_mV_ns"),
    }
    title, xlabel, limits, peak_key = settings[metric]
    manifest = json.loads((scan_root / "scan_manifest.json").read_text())
    fig, axes = plt.subplots(2, 4, figsize=(18, 9), sharex=True)
    rows = []
    for ax, point in zip(axes.flat, manifest["points"]):
        point_root = scan_root / "analysis" / point["tag"]
        matches = list(point_root.glob("*_PCB_CH3_ScopeA"))
        if len(matches) != 1:
            raise RuntimeError(
                f"Expected one PCB CH3/scope A signal in {point_root}, found {len(matches)}"
            )
        directory = matches[0]
        summary = json.loads((directory / "point_summary.json").read_text())
        data = pd.read_csv(directory / "pulse_measurements.csv")
        accepted = data.loc[data["accepted"].astype(bool) & ~data["overflow"].astype(bool), metric].to_numpy(float)
        ax.hist(accepted, bins=fd_edges(accepted, limits), histtype="step", color=COLOR, linewidth=1.25)
        for peak in summary.get(peak_key, []):
            ax.axvline(peak, color=RED, linewidth=1.15)
        count, total = summary["accepted_events"], summary["event_files"]
        ax.set_title(f"{point['effective_bias_v']:.3f} V\naccepted {count:,}/{total:,} ({100*count/total:.2f}%)")
        ax.set_yscale("log")
        ax.set_xlim(*limits)
        ax.grid(True, alpha=0.25)
        rows.append({"effective_bias_V": point["effective_bias_v"], "temperature_C": point["environment"]["temperature_C"], "accepted_events": count, "total_events": total})
    axes.flat[-1].axis("off")
    for ax in axes[:, 0]:
        ax.set_ylabel("Accepted events per bin")
    for ax in axes[-1, :]:
        if ax.axison:
            ax.set_xlabel(xlabel)
    fig.suptitle(f"SiPM 1 (△): PCB CH3, scope A, 25 mV trigger\n{title}", fontsize=17)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    output = "ch3_corrected_height_panels.png" if metric == "peak_height_mV" else "ch3_corrected_area_panels.png"
    fig.savefig(output_root / output, dpi=180, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(rows).to_csv(output_root / "ch3_corrected_scan_summary.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--threshold-root", type=Path, required=True)
    parser.add_argument("--scan-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    threshold_summary(args.threshold_root, args.output)
    corrected_scan_panel(args.scan_root, args.output, "peak_height_mV")
    corrected_scan_panel(args.scan_root, args.output, "pulse_area_mV_ns")


if __name__ == "__main__":
    main()
