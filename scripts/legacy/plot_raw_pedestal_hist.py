#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


RAIL_MV = 4000.0


def unit_scales_to_mV(unit_line: str) -> list[float]:
    parts = [part.strip().strip("()").lower() for part in unit_line.split(",")]
    scales = [1.0, 1.0, 1.0]
    for idx in [1, 2]:
        if idx < len(parts) and parts[idx] == "v":
            scales[idx] = 1000.0
        elif idx < len(parts) and parts[idx] == "uv":
            scales[idx] = 0.001
    return scales


def read_scope_csv(path: Path) -> np.ndarray | None:
    rows: list[list[float]] = []
    scales = [1.0, 1.0, 1.0]
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for idx, line in enumerate(fh):
            if idx < 3:
                if idx == 1:
                    scales = unit_scales_to_mV(line)
                continue
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            if len(parts) < 3:
                continue
            vals: list[float] = []
            ok = True
            for col, value in enumerate(parts[:3]):
                text = value.strip()
                if text in ("∞", "+∞", "inf", "+inf", "Inf", "+Inf", "Infinity", "+Infinity"):
                    vals.append(RAIL_MV)
                elif text in ("-∞", "-inf", "-Inf", "-Infinity"):
                    vals.append(-RAIL_MV)
                else:
                    try:
                        vals.append(float(text) * scales[col])
                    except ValueError:
                        ok = False
                        break
            if ok:
                rows.append(vals)
    if not rows:
        return None
    return np.asarray(rows, dtype=np.float64)


def fd_bins(values: np.ndarray, min_bins: int = 40, max_bins: int = 180) -> int:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 2:
        return 20
    q25, q75 = np.percentile(x, [25, 75])
    width = 2.0 * (q75 - q25) / np.cbrt(x.size)
    if not np.isfinite(width) or width <= 0:
        width = (np.nanmax(x) - np.nanmin(x)) / 60.0
    if not np.isfinite(width) or width <= 0:
        return 20
    return int(max(min_bins, min(max_bins, np.ceil((np.nanmax(x) - np.nanmin(x)) / width))))


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot raw voltage histograms without baseline subtraction.")
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--bias-v", required=True, type=float)
    parser.add_argument("--sipm", required=True)
    parser.add_argument("--detector-channel", required=True)
    parser.add_argument("--scope-channel", choices=["A", "B"], default="B")
    parser.add_argument("--max-files", type=int, default=20000)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    files = sorted(args.run_dir.glob("*.csv"))[: args.max_files]
    if not files:
        raise FileNotFoundError(f"No CSV files found in {args.run_dir}")

    col = 1 if args.scope_channel == "A" else 2
    pedestal_samples: list[np.ndarray] = []
    baseline_means: list[float] = []
    raw_peaks: list[float] = []
    raw_peak_times: list[float] = []
    raw_areas: list[float] = []

    for path in files:
        arr = read_scope_csv(path)
        if arr is None or arr.ndim != 2 or arr.shape[1] < 3:
            continue
        t = arr[:, 0]
        y = arr[:, col]
        finite = np.isfinite(y) & (np.abs(y) < 0.95 * RAIL_MV)
        pre = finite & (t < -30.0)
        if pre.sum() < 3:
            pre = finite & (t < -10.0)
        sig_win = finite & (t >= -10.0) & (t <= 260.0)
        area_win = finite & (t >= -5.0) & (t <= 260.0)
        if pre.any():
            pedestal_samples.append(y[pre])
            baseline_means.append(float(np.mean(y[pre])))
        if sig_win.any():
            yy = y[sig_win]
            tt = t[sig_win]
            i = int(np.argmax(yy))
            raw_peaks.append(float(yy[i]))
            raw_peak_times.append(float(tt[i]))
        if area_win.any():
            raw_areas.append(float(np.trapezoid(y[area_win], t[area_win])))

    pedestal = np.concatenate(pedestal_samples) if pedestal_samples else np.array([])
    baseline = np.asarray(baseline_means, dtype=float)
    peaks = np.asarray(raw_peaks, dtype=float)
    peak_times = np.asarray(raw_peak_times, dtype=float)
    areas = np.asarray(raw_areas, dtype=float)

    rows = {
        "run_name": args.run_dir.name,
        "bias_V": args.bias_v,
        "sipm": args.sipm,
        "detector_channel": args.detector_channel,
        "scope_channel": args.scope_channel,
        "files_used": len(files),
        "events_with_baseline": len(baseline),
        "events_with_peak": len(peaks),
        "pedestal_sample_count": len(pedestal),
        "pedestal_mean_mV": float(np.mean(pedestal)) if pedestal.size else np.nan,
        "pedestal_median_mV": float(np.median(pedestal)) if pedestal.size else np.nan,
        "pedestal_rms_mV": float(np.std(pedestal, ddof=1)) if pedestal.size > 1 else np.nan,
        "event_baseline_mean_mV": float(np.mean(baseline)) if baseline.size else np.nan,
        "event_baseline_std_mV": float(np.std(baseline, ddof=1)) if baseline.size > 1 else np.nan,
        "raw_peak_mean_mV": float(np.mean(peaks)) if peaks.size else np.nan,
        "raw_peak_median_mV": float(np.median(peaks)) if peaks.size else np.nan,
        "raw_peak_std_mV": float(np.std(peaks, ddof=1)) if peaks.size > 1 else np.nan,
        "raw_peak_time_median_ns": float(np.median(peak_times)) if peak_times.size else np.nan,
        "raw_area_mean_mVns": float(np.mean(areas)) if areas.size else np.nan,
        "raw_area_median_mVns": float(np.median(areas)) if areas.size else np.nan,
    }
    pd.DataFrame([rows]).to_csv(args.out_dir / "raw_voltage_pedestal_summary.csv", index=False)

    plt.style.use("default")
    fig, axes = plt.subplots(2, 1, figsize=(11, 8.2), dpi=160)
    fig.suptitle(f"{args.sipm.upper()} / {args.detector_channel.upper()} raw voltage at {args.bias_v:.2f} V")

    ax = axes[0]
    ax.hist(pedestal, bins=fd_bins(pedestal, 60, 220), color="C1", alpha=0.72, edgecolor="white", linewidth=0.4)
    ax.axvline(np.mean(pedestal), color="black", linestyle="--", linewidth=1.4, label=f"mean = {np.mean(pedestal):.2f} mV")
    ax.axvline(np.median(pedestal), color="0.35", linestyle=":", linewidth=1.4, label=f"median = {np.median(pedestal):.2f} mV")
    ax.set_title(f"Pedestal: raw pre-trigger samples on scope Channel {args.scope_channel}")
    ax.set_xlabel("Raw voltage (mV), no baseline subtraction")
    ax.set_ylabel("Samples per bin")
    ax.grid(True, alpha=0.28)
    ax.legend(fontsize=9)

    ax = axes[1]
    ax.hist(peaks, bins=fd_bins(peaks, 60, 220), color="C1", alpha=0.72, edgecolor="white", linewidth=0.4)
    ax.axvline(np.mean(peaks), color="black", linestyle="--", linewidth=1.4, label=f"mean = {np.mean(peaks):.2f} mV")
    ax.axvline(np.median(peaks), color="0.35", linestyle=":", linewidth=1.4, label=f"median = {np.median(peaks):.2f} mV")
    ax.set_title("Raw pulse-peak distribution")
    ax.set_xlabel("Raw pulse peak (mV), no baseline subtraction")
    ax.set_ylabel("Events per bin")
    ax.grid(True, alpha=0.28)
    ax.legend(fontsize=9)

    fig.tight_layout(rect=[0, 0, 1, 0.955])
    fig.savefig(args.out_dir / "raw_voltage_pedestal_hist.png")
    plt.close(fig)

    print(args.out_dir)
    print(pd.DataFrame([rows]).T.to_string())


if __name__ == "__main__":
    main()
