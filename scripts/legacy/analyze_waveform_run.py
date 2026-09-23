#!/usr/bin/env python3
from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.signal import find_peaks


RAIL_MV = 4000.0
RAIL_V = RAIL_MV / 1000.0


def unit_scales_to_mV(unit_line: str) -> list[float]:
    parts = [part.strip().strip("()").lower() for part in unit_line.split(",")]
    scales = [1.0, 1.0, 1.0]
    for idx in [1, 2]:
        if idx < len(parts) and parts[idx] == "v":
            scales[idx] = 1000.0
        elif idx < len(parts) and parts[idx] == "uv":
            scales[idx] = 0.001
    return scales


def read_scope_csv(path: Path) -> tuple[np.ndarray | None, list[int]]:
    rows: list[list[float]] = []
    overflow_cols = [0, 0, 0]
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
                    overflow_cols[col] += 1
                elif text in ("-∞", "-inf", "-Inf", "-Infinity"):
                    vals.append(-RAIL_MV)
                    overflow_cols[col] += 1
                else:
                    try:
                        vals.append(float(text) * scales[col])
                    except ValueError:
                        ok = False
                        break
            if ok:
                rows.append(vals)
    if not rows:
        return None, overflow_cols
    return np.asarray(rows, dtype=np.float64), overflow_cols


def fd_bins(values: np.ndarray, min_bins: int = 40, max_bins: int = 180) -> int:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 2:
        return 20
    q25, q75 = np.percentile(x, [25, 75])
    width = 2 * (q75 - q25) / np.cbrt(x.size)
    if not np.isfinite(width) or width <= 0:
        width = (np.nanmax(x) - np.nanmin(x)) / 60.0
    if not np.isfinite(width) or width <= 0:
        return 20
    return int(max(min_bins, min(max_bins, np.ceil((np.nanmax(x) - np.nanmin(x)) / width))))


def gaussian_bg(x: np.ndarray, amp: float, mu: float, sigma: float, bg: float) -> np.ndarray:
    return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2) + bg


def select_even_peak_sequence(fits: pd.DataFrame) -> tuple[list[int], float, float, float]:
    if len(fits) < 3:
        return [], np.nan, np.nan, np.nan
    mus = fits["mu_mV"].to_numpy(float)
    best: tuple[float, tuple[int, ...], float, float] | None = None
    for count in range(min(5, len(mus)), 2, -1):
        for combo in combinations(range(len(mus)), count):
            vals = mus[list(combo)]
            diffs = np.diff(vals)
            median_diff = np.median(diffs)
            if median_diff < 8 or median_diff > 120:
                continue
            residual = float(np.sqrt(np.mean((diffs - median_diff) ** 2)))
            score = residual / median_diff - 0.02 * count
            if best is None or score < best[0]:
                best = (float(score), combo, float(median_diff), residual)
        if best is not None:
            break
    if best is None:
        return [], np.nan, np.nan, np.nan
    _, combo, _, residual = best
    selected = list(combo)
    diffs = np.diff(mus[selected])
    spacing = float(np.mean(diffs))
    spacing_unc = float(np.std(diffs, ddof=1) / np.sqrt(len(diffs))) if len(diffs) > 1 else np.nan
    return selected, spacing, spacing_unc, residual


def fit_height_peaks(heights: np.ndarray, out_csv: Path) -> tuple[pd.DataFrame, float, float, float, str]:
    hist_bins = fd_bins(heights, 50, 220)
    counts, edges = np.histogram(heights, bins=hist_bins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    kernel = np.array([1, 2, 3, 2, 1], dtype=float)
    kernel /= kernel.sum()
    smooth = np.convolve(counts, kernel, mode="same")
    bin_width = float(np.mean(np.diff(centers)))
    prominence = max(6.0, 0.035 * np.nanmax(smooth))
    peak_indexes, _ = find_peaks(smooth, prominence=prominence, distance=max(2, int(8 / bin_width)))

    fit_rows: list[dict[str, float]] = []
    for idx in peak_indexes:
        mu0 = float(centers[idx])
        half_width = max(6.0, 4 * bin_width)
        mask = (centers >= mu0 - half_width) & (centers <= mu0 + half_width)
        if mask.sum() < 5:
            continue
        x = centers[mask]
        y = counts[mask]
        amp0 = max(1.0, float(counts[idx] - np.percentile(y, 10)))
        sigma0 = max(bin_width * 1.5, 1.5)
        bg0 = max(0.0, float(np.percentile(y, 10)))
        try:
            popt, pcov = curve_fit(
                gaussian_bg,
                x,
                y,
                p0=[amp0, mu0, sigma0, bg0],
                bounds=([0, mu0 - half_width, 0.3, 0], [np.inf, mu0 + half_width, 30.0, np.inf]),
                maxfev=10000,
            )
            perr = np.sqrt(np.diag(pcov))
            fit_rows.append({"mu_mV": popt[1], "mu_err_mV": perr[1], "sigma_mV": abs(popt[2]), "amp": popt[0], "bg": popt[3]})
        except Exception:
            fit_rows.append({"mu_mV": mu0, "mu_err_mV": np.nan, "sigma_mV": np.nan, "amp": counts[idx], "bg": np.nan})

    fits = pd.DataFrame(fit_rows).sort_values("mu_mV") if fit_rows else pd.DataFrame(columns=["mu_mV", "mu_err_mV", "sigma_mV", "amp", "bg"])
    fits["selected_for_1pe_spacing"] = False
    fits["selected_index"] = np.nan
    selected, spacing, spacing_unc, residual = select_even_peak_sequence(fits)
    for selected_index, fit_index in enumerate(selected, start=1):
        fits.loc[fits.index[fit_index], "selected_for_1pe_spacing"] = True
        fits.loc[fits.index[fit_index], "selected_index"] = selected_index
    selected_mus = ";".join(f"{fits.iloc[i]['mu_mV']:.3f}" for i in selected)
    fits.to_csv(out_csv, index=False)
    return fits, spacing, spacing_unc, residual, selected_mus


def add_standard_labels(ax, title: str, xlabel: str, ylabel: str) -> None:
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.28)


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze PicoScope waveform CSV folders.")
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--bias-v", required=True, type=float)
    parser.add_argument("--sipm", required=True)
    parser.add_argument("--detector-channel", required=True)
    parser.add_argument("--color", default="C0")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    files = sorted(args.run_dir.glob("*.csv"))
    if not files:
        raise FileNotFoundError(f"No CSV files found in {args.run_dir}")

    rows: list[dict[str, object]] = []
    wave_a: list[np.ndarray] = []
    wave_b: list[np.ndarray] = []
    time_ns: np.ndarray | None = None

    for event, path in enumerate(files, start=1):
        arr, overflow = read_scope_csv(path)
        if arr is None or arr.ndim != 2 or arr.shape[1] < 3:
            rows.append({"event": event, "file": path.name, "ok_shape": False, "accepted": False})
            continue

        t = arr[:, 0]
        if time_ns is None:
            time_ns = t
        pre = t < -30.0
        if not np.any(pre):
            pre = t < -10.0
        sig_win = (t >= -10.0) & (t <= 260.0)
        area_win = (t >= -5.0) & (t <= 260.0)
        if not np.any(sig_win):
            sig_win = np.ones_like(t, dtype=bool)
        if not np.any(area_win):
            area_win = sig_win

        rec: dict[str, object] = {
            "event": event,
            "file": path.name,
            "ok_shape": True,
            "n_samples": len(t),
            "t_min_ns": float(np.nanmin(t)),
            "t_max_ns": float(np.nanmax(t)),
            "A_overflow_samples": int(overflow[1]),
            "B_overflow_samples": int(overflow[2]),
            "any_overflow_samples": int(overflow[1] + overflow[2]),
        }
        for label, col in [("A", 1), ("B", 2)]:
            y = arr[:, col]
            baseline_values = y[pre]
            baseline_values = baseline_values[np.abs(baseline_values) < RAIL_MV * 0.95]
            if baseline_values.size:
                baseline = float(np.mean(baseline_values))
            else:
                finite = y[np.abs(y) < RAIL_MV * 0.95]
                baseline = float(np.median(finite)) if finite.size else 0.0
            yy_mV = y - baseline
            yy_v = yy_mV / 1000.0
            if label == "A":
                wave_a.append(yy_v.astype(np.float32))
            else:
                wave_b.append(yy_v.astype(np.float32))
            ysig = yy_mV[sig_win]
            tsig = t[sig_win]
            peak_idx = int(np.argmax(ysig))
            trough_idx = int(np.argmin(ysig))
            rms = float(np.std(yy_mV[pre], ddof=1)) if np.sum(pre) > 1 else np.nan
            peak = float(ysig[peak_idx])
            rec[f"{label}_baseline_mV"] = baseline
            rec[f"{label}_rms_mV"] = rms
            rec[f"{label}_peak_mV"] = peak
            rec[f"{label}_peak_time_ns"] = float(tsig[peak_idx])
            rec[f"{label}_trough_mV"] = float(ysig[trough_idx])
            rec[f"{label}_trough_time_ns"] = float(tsig[trough_idx])
            rec[f"{label}_area_mVns"] = float(np.trapezoid(yy_mV[area_win], t[area_win]))
            rec[f"{label}_snr"] = peak / rms if rms and np.isfinite(rms) else np.nan
        rows.append(rec)

    if time_ns is None:
        raise RuntimeError("No valid waveform files were parsed.")

    meas = pd.DataFrame(rows)
    valid = meas["ok_shape"].fillna(False).astype(bool)
    median_a = meas.loc[valid & meas["A_overflow_samples"].fillna(0).eq(0), "A_peak_mV"].median()
    median_b = meas.loc[valid & meas["B_overflow_samples"].fillna(0).eq(0), "B_peak_mV"].median()
    active = "A" if median_a >= median_b else "B"
    passive = "B" if active == "A" else "A"
    active_overflow = f"{active}_overflow_samples"
    active_peak = f"{active}_peak_mV"
    active_rms = f"{active}_rms_mV"
    active_time = f"{active}_peak_time_ns"

    meas["active_scope_channel"] = active
    meas["passive_scope_channel"] = passive
    meas["clipped"] = meas[active_overflow].fillna(0).gt(0) | meas["any_overflow_samples"].fillna(0).gt(0)
    meas["accepted"] = (
        meas["ok_shape"].fillna(False).astype(bool)
        & ~meas["clipped"]
        & meas[active_peak].gt(5 * meas[active_rms])
        & meas[active_time].between(-5, 170)
    )
    meas.to_csv(args.out_dir / "pulse_measurements.csv", index=False)

    accepted = meas["accepted"].fillna(False).astype(bool)
    heights = meas.loc[accepted, active_peak].to_numpy(float)
    areas = meas.loc[accepted, f"{active}_area_mVns"].to_numpy(float)
    rms_values = meas.loc[accepted, active_rms].to_numpy(float)
    peak_times = meas.loc[accepted, active_time].to_numpy(float)

    summary: dict[str, object] = {
        "run_name": args.run_dir.name,
        "bias_V": args.bias_v,
        "sipm": args.sipm,
        "detector_channel": args.detector_channel,
        "run_dir": str(args.run_dir),
        "analysis_dir": str(args.out_dir),
        "n_files": len(files),
        "n_events_measured": len(meas),
        "n_accepted": int(accepted.sum()),
        "active_scope_channel": active,
        "median_A_peak_mV": float(median_a),
        "median_B_peak_mV": float(median_b),
        "overflow_events_active": int(meas[active_overflow].fillna(0).gt(0).sum()),
        "overflow_samples_active": int(meas[active_overflow].fillna(0).sum()),
        "overflow_events_any": int(meas["any_overflow_samples"].fillna(0).gt(0).sum()),
        "overflow_samples_any": int(meas["any_overflow_samples"].fillna(0).sum()),
    }
    for channel in ["A", "B"]:
        for metric in ["peak_mV", "area_mVns", "rms_mV", "peak_time_ns", "snr"]:
            col = f"{channel}_{metric}"
            x = meas.loc[accepted, col].dropna().to_numpy(float) if col in meas else np.array([])
            if x.size:
                summary[f"{col}_mean"] = float(np.mean(x))
                summary[f"{col}_median"] = float(np.median(x))
                summary[f"{col}_std"] = float(np.std(x, ddof=1)) if x.size > 1 else np.nan

    fits, spacing, spacing_unc, residual, selected_mus = fit_height_peaks(heights, args.out_dir / "pulse_height_refined_peak_fits.csv")
    summary.update(
        {
            "fit_peak_count": int(len(fits)),
            "refined_selected_peak_count": int(fits["selected_for_1pe_spacing"].sum()) if len(fits) else 0,
            "refined_selected_peak_mus_mV": selected_mus,
            "refined_one_pe_height_mV": spacing,
            "refined_one_pe_height_unc_mV": spacing_unc,
            "refined_sequence_residual_rms_mV": residual,
        }
    )

    plt.style.use("default")
    wave_a_np = np.asarray(wave_a, dtype=np.float32)
    wave_b_np = np.asarray(wave_b, dtype=np.float32)
    active_waves = wave_a_np if active == "A" else wave_b_np
    active_color = args.color

    fig, ax = plt.subplots(figsize=(11, 6), dpi=160)
    for idx in np.linspace(0, len(meas) - 1, min(60, len(meas)), dtype=int):
        if idx < len(wave_a_np):
            ax.plot(time_ns, wave_a_np[idx] * 1000, color="C0", alpha=0.12, linewidth=0.55)
        if idx < len(wave_b_np):
            ax.plot(time_ns, wave_b_np[idx] * 1000, color="C1", alpha=0.08, linewidth=0.45)
    ax.plot([], [], color="C0", label="Channel A")
    ax.plot([], [], color="C1", label="Channel B")
    ax.axhline(0, color="0.35", linestyle="--", linewidth=0.9)
    add_standard_labels(ax, f"{args.run_dir.name}: example baseline-subtracted waveforms", "Time (ns)", "Voltage (mV)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.out_dir / "example_waveforms.png")
    plt.close(fig)

    accepted_mask = accepted.to_numpy()
    accepted_waves = active_waves[accepted_mask] if accepted_mask.any() else active_waves
    mean_wave = np.mean(accepted_waves, axis=0)
    p05 = np.percentile(accepted_waves, 5, axis=0)
    p95 = np.percentile(accepted_waves, 95, axis=0)
    fig, ax = plt.subplots(figsize=(11, 6), dpi=160)
    ax.fill_between(time_ns, p05 * 1000, p95 * 1000, color=active_color, alpha=0.18, label="5-95% band")
    ax.plot(time_ns, mean_wave * 1000, color="black", linewidth=1.5, label="mean waveform")
    ax.axhline(0, color="0.35", linestyle="--", linewidth=0.9)
    add_standard_labels(ax, f"{args.run_dir.name}: active {active} mean waveform", "Time (ns)", "Voltage (mV)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.out_dir / "mean_waveform.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 6.4), dpi=160)
    hist_bins = fd_bins(heights, 50, 220)
    ax.hist(heights, bins=hist_bins, color=active_color, alpha=0.72, edgecolor="white", linewidth=0.45, label=f"accepted = {len(heights)}")
    ax.axvline(np.mean(heights), color="black", linestyle="--", linewidth=1.5, label=f"mean = {np.mean(heights):.2f} mV")
    ax.axvline(np.median(heights), color="0.35", linestyle=":", linewidth=1.5, label=f"median = {np.median(heights):.2f} mV")
    peak_colors = ["C3", "C4", "C2", "C5", "C6"]
    for j, (_, row) in enumerate(fits[fits["selected_for_1pe_spacing"]].iterrows(), start=1):
        ax.axvline(row["mu_mV"], color=peak_colors[(j - 1) % len(peak_colors)], linewidth=1.2, label=f"peak {j}: {row['mu_mV']:.2f} mV")
    add_standard_labels(ax, f"{args.sipm.upper()} / {args.detector_channel.upper()} pulse-height histogram at {args.bias_v:.2f} V", f"Pulse height on scope Channel {active} (mV)", "Events per bin")
    if np.isfinite(spacing):
        title = f"1 p.e. spacing = {spacing:.2f} +/- {spacing_unc:.2f} mV" if np.isfinite(spacing_unc) else f"1 p.e. spacing = {spacing:.2f} mV"
        ax.legend(title=title, fontsize=9, title_fontsize=9)
    else:
        ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(args.out_dir / "pulse_height_hist.png")
    plt.close(fig)

    counts, edges = np.histogram(heights, bins=hist_bins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    kernel = np.array([1, 2, 3, 2, 1], dtype=float)
    kernel /= kernel.sum()
    smooth = np.convolve(counts, kernel, mode="same")
    fig, ax = plt.subplots(figsize=(11, 6), dpi=160)
    ax.step(centers, counts, where="mid", color=active_color, linewidth=1.1, alpha=0.85, label="histogram")
    ax.plot(centers, smooth, color="black", linewidth=1.2, label="smoothed")
    if len(fits):
        ax.scatter(fits["mu_mV"], np.interp(fits["mu_mV"], centers, smooth), marker="x", color="C3", s=45, label="fit peaks")
        for _, row in fits[fits["selected_for_1pe_spacing"]].iterrows():
            ax.axvline(row["mu_mV"], color="C2", linestyle="--", alpha=0.8)
    add_standard_labels(ax, f"{args.sipm.upper()} / {args.detector_channel.upper()} pulse-height peak search at {args.bias_v:.2f} V", "Pulse height (mV)", "Events per bin")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(args.out_dir / "pulse_height_peak_search.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 6), dpi=160)
    ax.hist(areas, bins=fd_bins(areas), color=active_color, alpha=0.72, edgecolor="white", linewidth=0.45)
    ax.axvline(np.mean(areas), color="black", linestyle="--", linewidth=1.3, label=f"mean = {np.mean(areas):.1f} mV ns")
    ax.axvline(np.median(areas), color="0.35", linestyle=":", linewidth=1.3, label=f"median = {np.median(areas):.1f} mV ns")
    add_standard_labels(ax, f"{args.sipm.upper()} / {args.detector_channel.upper()} pulse-area histogram at {args.bias_v:.2f} V", "Pulse area (mV ns)", "Events per bin")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(args.out_dir / "pulse_area_hist.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 6), dpi=160)
    ax.hist(rms_values, bins=fd_bins(rms_values), color=active_color, alpha=0.72, edgecolor="white", linewidth=0.45)
    ax.axvline(np.median(rms_values), color="black", linestyle="--", label=f"median = {np.median(rms_values):.2f} mV")
    add_standard_labels(ax, f"{args.sipm.upper()} / {args.detector_channel.upper()} baseline noise at {args.bias_v:.2f} V", "Baseline RMS (mV)", "Events per bin")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(args.out_dir / "baseline_noise_hist.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 6), dpi=160)
    ax.hist(peak_times, bins=fd_bins(peak_times), color=active_color, alpha=0.72, edgecolor="white", linewidth=0.45)
    ax.axvline(np.median(peak_times), color="black", linestyle="--", label=f"median = {np.median(peak_times):.2f} ns")
    add_standard_labels(ax, f"{args.sipm.upper()} / {args.detector_channel.upper()} peak-time histogram at {args.bias_v:.2f} V", "Peak time (ns)", "Events per bin")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(args.out_dir / "peak_time_hist.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 6.5), dpi=160)
    ax.scatter(areas, heights, s=6, alpha=0.22, color=active_color, edgecolors="none")
    add_standard_labels(ax, f"{args.sipm.upper()} / {args.detector_channel.upper()} pulse height vs area at {args.bias_v:.2f} V", "Pulse area (mV ns)", "Pulse height (mV)")
    fig.tight_layout()
    fig.savefig(args.out_dir / "area_vs_height.png")
    plt.close(fig)

    peaks_v_all = meas[f"{active}_peak_mV"].to_numpy(float) / 1000.0
    integrals_vns_all = meas[f"{active}_area_mVns"].to_numpy(float) / 1000.0
    troughs_v_all = meas[f"{active}_trough_mV"].to_numpy(float) / 1000.0
    mean_all = np.mean(active_waves, axis=0)
    p05_all = np.percentile(active_waves, 5, axis=0)
    p95_all = np.percentile(active_waves, 95, axis=0)
    overflow_events = int(meas[active_overflow].fillna(0).gt(0).sum())
    full_ymin = min(-0.20, float(np.percentile(troughs_v_all, 0.1) * 1.15))
    full_ymax = max(float(np.percentile(peaks_v_all, 99.9) * 1.10), float(np.max(peaks_v_all) * 1.03), 0.20)
    if overflow_events:
        full_ymax = max(full_ymax, 4.15)
    normal_peaks = peaks_v_all[peaks_v_all < 0.95 * RAIL_V]
    zoom_ymin = min(-0.05, float(np.percentile(troughs_v_all, 0.2) * 1.2))
    zoom_ymax = max(0.05, float(np.percentile(normal_peaks, 99.5) * 1.18) if normal_peaks.size else float(np.percentile(peaks_v_all, 99.5) * 1.18))

    fig, axes = plt.subplots(3, 1, figsize=(12.8, 10.0), dpi=160, gridspec_kw={"height_ratios": [1.0, 1.9, 1.25]})
    fig.suptitle(f"{args.sipm.upper()} / {args.detector_channel.upper()}: waveform saturation check at {args.bias_v:.2f} V", fontsize=15)
    ax = axes[0]
    ax.hist(integrals_vns_all, bins=fd_bins(integrals_vns_all), histtype="step", linewidth=1.4, color=active_color)
    ax.axvline(np.mean(integrals_vns_all), color="black", linestyle="--", linewidth=1.1, label=f"mean = {np.mean(integrals_vns_all):.2f} V ns")
    add_standard_labels(ax, f"Channel {active}: waveform integral", "Integrated area after baseline subtraction (V ns)", "Events")
    ax.legend(fontsize=8, loc="upper right")
    ax = axes[1]
    ax.plot(time_ns, active_waves.T, color=active_color, alpha=0.016, linewidth=0.25, rasterized=True)
    ax.fill_between(time_ns, p05_all, p95_all, color=active_color, alpha=0.16, linewidth=0, label="5-95% band")
    ax.plot(time_ns, mean_all, color="black", linewidth=1.3, label="mean waveform")
    ax.axhline(0, color="0.35", linestyle="--", linewidth=0.9)
    if overflow_events:
        ax.axhline(RAIL_V, color="crimson", linestyle=":", linewidth=1.1, label="scope over-range marker")
    add_standard_labels(ax, f"Channel {active}: baseline-subtracted waveforms ({len(active_waves)} plotted)", "Time (ns)", "Voltage (V)")
    ax.set_ylim(full_ymin, full_ymax)
    ax.legend(fontsize=8, loc="upper right")
    ax = axes[2]
    ax.scatter(integrals_vns_all, peaks_v_all, s=5, alpha=0.17, color=active_color, edgecolors="none")
    ax.axhline(float(np.max(peaks_v_all)), color="black", linestyle=":", linewidth=1.0, label=f"max peak = {np.max(peaks_v_all):.3f} V")
    if overflow_events:
        ax.axhline(RAIL_V, color="crimson", linestyle=":", linewidth=1.1, label="over-range level")
    add_standard_labels(ax, f"Channel {active}: peak height vs waveform integral", "Integrated area after baseline subtraction (V ns)", "Peak height (V)")
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    fig.savefig(args.out_dir / "saturation_check_panel.png")
    plt.close(fig)

    for suffix, ylim in [("full_range", (full_ymin, full_ymax)), ("zoom", (zoom_ymin, zoom_ymax))]:
        fig, ax = plt.subplots(figsize=(12.8, 5.8), dpi=160)
        ax.plot(time_ns, active_waves.T, color=active_color, alpha=0.016, linewidth=0.25, rasterized=True)
        ax.fill_between(time_ns, p05_all, p95_all, color=active_color, alpha=0.16, linewidth=0, label="5-95% band")
        ax.plot(time_ns, mean_all, color="black", linewidth=1.4, label="mean waveform")
        ax.axhline(0, color="0.35", linestyle="--", linewidth=0.9)
        if suffix == "full_range" and overflow_events:
            ax.axhline(RAIL_V, color="crimson", linestyle=":", linewidth=1.1, label="scope over-range marker")
        add_standard_labels(ax, f"{args.sipm.upper()} / {args.detector_channel.upper()}: all baseline-subtracted waveforms ({suffix.replace('_', ' ')})", "Time (ns)", "Voltage (V)")
        ax.set_ylim(*ylim)
        ax.legend(fontsize=9, loc="upper right")
        fig.tight_layout()
        fig.savefig(args.out_dir / f"all_waveforms_overlay_{suffix}.png")
        plt.close(fig)

    summary["peak_max_V_all"] = float(np.max(peaks_v_all))
    summary["peak_99p9_V_all"] = float(np.percentile(peaks_v_all, 99.9))
    summary["near_rail_saturation_candidates"] = overflow_events
    summary["flat_top_candidates"] = int((peaks_v_all >= 0.95 * RAIL_V).sum())
    pd.DataFrame([summary]).to_csv(args.out_dir / "summary.csv", index=False)
    print(args.out_dir)
    print(pd.DataFrame([summary]).T.to_string())


if __name__ == "__main__":
    main()
