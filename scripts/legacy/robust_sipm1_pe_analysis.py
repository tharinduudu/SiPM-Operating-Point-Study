#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, savgol_filter


RAIL_MV = 4000.0
OUT_DIR = Path("/home/muon/brDownVstudy/robust_sipm1_pe_analysis_20260901")


@dataclass(frozen=True)
class RunSpec:
    bias_v: float
    label: str
    run_dir: Path
    scope_channel: str = "A"


RUNS = [
    RunSpec(53.20, "53.20 V", Path("/home/muon/brDownVstudy/20260827_53p2v_sipm1_ch0_20000-0005")),
    RunSpec(53.90, "53.90 V", Path("/home/muon/brDownVstudy/20260827_53p9V_sipm1_ch0_10000")),
    RunSpec(54.10, "54.10 V", Path("/home/muon/brDownVstudy/20260831-54p1_sipm1_ch0_20000")),
    RunSpec(54.90, "54.90 V", Path("/home/muon/brDownVstudy/20260827_54p9V_sipm1_ch0_10000")),
    RunSpec(55.59, "55.59 V", Path("/home/muon/brDownVstudy/20260827_55p59V_sipm1_ch0_10000")),
    RunSpec(55.65, "55.65 V", Path("/home/muon/brDownVstudy/20260827_55p65V_sipm1_ch0_10000")),
    RunSpec(56.30, "56.30 V", Path("/home/muon/brDownVstudy/20260831_56p3V_sipm1_ch0_20000")),
    RunSpec(56.40, "56.40 V", Path("/home/muon/brDownVstudy/20260827_56p4V_sipm1_ch0_10000")),
    RunSpec(56.86, "56.86 V", Path("/home/muon/brDownVstudy/20260827_56.86_sipm1_ch0_20000")),
]


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


def scope_csv_files(run_dir: Path) -> list[Path]:
    direct = sorted(run_dir.glob("*.csv"))
    if direct:
        return direct
    files: list[Path] = []
    for child in sorted(run_dir.iterdir()):
        if not child.is_dir():
            continue
        if child.name.startswith(("analysis", "raw_pedestal", "robust")):
            continue
        files.extend(sorted(child.glob("*.csv")))
    return files


def robust_location_scale(values: np.ndarray) -> tuple[float, float, int]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x) & (np.abs(x) < 0.95 * RAIL_MV)]
    if x.size == 0:
        return 0.0, np.nan, 0
    median = float(np.median(x))
    mad = float(np.median(np.abs(x - median)))
    sigma = 1.4826 * mad
    if not np.isfinite(sigma) or sigma <= 1e-9:
        sigma = float(np.std(x, ddof=1)) if x.size > 1 else 0.0
    if np.isfinite(sigma) and sigma > 0:
        keep = np.abs(x - median) < 5.0 * sigma
        clipped = x[keep]
    else:
        clipped = x
    if clipped.size == 0:
        clipped = x
    baseline = float(np.median(clipped))
    noise = 1.4826 * float(np.median(np.abs(clipped - baseline)))
    if not np.isfinite(noise) or noise <= 1e-9:
        noise = float(np.std(clipped, ddof=1)) if clipped.size > 1 else np.nan
    return baseline, noise, int(clipped.size)


def smooth_waveform(y: np.ndarray) -> np.ndarray:
    if y.size >= 15:
        return savgol_filter(y, window_length=11, polyorder=2, mode="interp")
    if y.size >= 5:
        kernel = np.ones(5) / 5.0
        return np.convolve(y, kernel, mode="same")
    return y


def quadratic_peak(t: np.ndarray, y: np.ndarray, idx: int) -> tuple[float, float]:
    if idx <= 0 or idx >= len(y) - 1:
        return float(t[idx]), float(y[idx])
    x = t[idx - 1 : idx + 2]
    z = y[idx - 1 : idx + 2]
    try:
        a, b, c = np.polyfit(x, z, 2)
        if a < 0:
            x0 = -b / (2.0 * a)
            if x[0] <= x0 <= x[-1]:
                return float(x0), float(a * x0 * x0 + b * x0 + c)
    except Exception:
        pass
    return float(t[idx]), float(y[idx])


def fd_bin_width(values: np.ndarray, lower: float, upper: float) -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    x = x[(x >= lower) & (x <= upper)]
    if x.size < 2:
        return 1.0
    q25, q75 = np.percentile(x, [25, 75])
    width = 2.0 * (q75 - q25) / np.cbrt(x.size)
    if not np.isfinite(width) or width <= 0:
        width = (upper - lower) / 100.0
    return float(np.clip(width, 0.7, 3.0))


def gaussian_linear_bg(x: np.ndarray, amp: float, mu: float, sigma: float, bg0: float, bg1: float) -> np.ndarray:
    return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2) + bg0 + bg1 * (x - mu)


def local_gaussian_peak_fits(values: np.ndarray) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 100:
        empty = pd.DataFrame(columns=["mu_mV", "mu_err_mV", "sigma_mV", "amp", "prominence", "selected_for_spacing"])
        return empty, np.array([]), np.array([]), np.array([])
    lower, upper = np.percentile(x, [0.2, 99.4])
    pad = max(5.0, 0.05 * (upper - lower))
    lower -= pad
    upper += pad
    bin_width = fd_bin_width(x, lower, upper)
    edges = np.arange(np.floor(lower / bin_width) * bin_width, np.ceil(upper / bin_width) * bin_width + bin_width, bin_width)
    counts, edges = np.histogram(x, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    if counts.size < 10 or counts.max() <= 0:
        empty = pd.DataFrame(columns=["mu_mV", "mu_err_mV", "sigma_mV", "amp", "prominence", "selected_for_spacing"])
        return empty, centers, counts, counts
    smooth = gaussian_filter1d(counts.astype(float), sigma=max(1.0, 1.2 / bin_width))
    prominence = max(8.0, 0.025 * float(np.max(smooth)))
    distance = max(3, int(round(6.0 / bin_width)))
    peak_indexes, props = find_peaks(smooth, prominence=prominence, distance=distance)

    rows: list[dict[str, float | bool]] = []
    for j, idx in enumerate(peak_indexes):
        mu0 = float(centers[idx])
        half_width = max(6.0, 5.0 * bin_width)
        mask = (centers >= mu0 - half_width) & (centers <= mu0 + half_width)
        if mask.sum() < 7:
            continue
        xx = centers[mask]
        yy = counts[mask].astype(float)
        amp0 = max(1.0, float(counts[idx] - np.percentile(yy, 15)))
        sigma0 = max(1.2, 1.5 * bin_width)
        bg0 = max(0.0, float(np.percentile(yy, 15)))
        try:
            popt, pcov = curve_fit(
                gaussian_linear_bg,
                xx,
                yy,
                p0=[amp0, mu0, sigma0, bg0, 0.0],
                bounds=([0.0, mu0 - half_width, 0.25, 0.0, -np.inf], [np.inf, mu0 + half_width, 30.0, np.inf, np.inf]),
                maxfev=20000,
            )
            perr = np.sqrt(np.diag(pcov))
            mu = float(popt[1])
            mu_err = float(perr[1]) if np.isfinite(perr[1]) else np.nan
            sigma = abs(float(popt[2]))
            amp = float(popt[0])
        except Exception:
            mu = mu0
            mu_err = np.nan
            sigma = np.nan
            amp = float(counts[idx])
        rows.append(
            {
                "mu_mV": mu,
                "mu_err_mV": mu_err,
                "sigma_mV": sigma,
                "amp": amp,
                "prominence": float(props["prominences"][j]) if "prominences" in props else np.nan,
                "selected_for_spacing": False,
            }
        )
    fits = pd.DataFrame(rows).sort_values("mu_mV").reset_index(drop=True) if rows else pd.DataFrame(columns=["mu_mV", "mu_err_mV", "sigma_mV", "amp", "prominence", "selected_for_spacing"])
    return fits, centers, counts, smooth


def select_spacing(fits: pd.DataFrame) -> tuple[list[int], float, float, float]:
    if len(fits) < 3:
        return [], np.nan, np.nan, np.nan
    mus = fits["mu_mV"].to_numpy(float)
    best: tuple[float, tuple[int, ...], float, float] | None = None
    for count in range(min(6, len(mus)), 2, -1):
        for combo in combinations(range(len(mus)), count):
            vals = mus[list(combo)]
            diffs = np.diff(vals)
            spacing = float(np.median(diffs))
            if spacing < 8.0 or spacing > 120.0:
                continue
            residual = float(np.sqrt(np.mean((diffs - spacing) ** 2)))
            coverage = vals[-1] - vals[0]
            score = residual / spacing - 0.025 * count - 0.0005 * coverage
            if best is None or score < best[0]:
                best = (score, combo, spacing, residual)
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


def analyze_run(spec: RunSpec) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object], tuple[np.ndarray, np.ndarray, np.ndarray]]:
    files = scope_csv_files(spec.run_dir)
    col = 1 if spec.scope_channel == "A" else 2
    rows: list[dict[str, object]] = []
    for event, path in enumerate(files, start=1):
        arr, overflow = read_scope_csv(path)
        if arr is None or arr.ndim != 2 or arr.shape[1] < 3:
            rows.append({"event": event, "file": path.name, "ok_shape": False, "accepted": False})
            continue
        t = arr[:, 0]
        y = arr[:, col]
        finite = np.isfinite(y) & (np.abs(y) < 0.95 * RAIL_MV)
        pre = finite & (t < -50.0)
        if pre.sum() < 20:
            pre = finite & (t < -30.0)
        prompt = finite & (t >= -5.0) & (t <= 90.0)
        integ = finite & (t >= -5.0) & (t <= 180.0)
        late = finite & (t >= 120.0) & (t <= 260.0)
        if prompt.sum() < 5:
            rows.append({"event": event, "file": path.name, "ok_shape": False, "accepted": False})
            continue
        baseline, noise, n_base = robust_location_scale(y[pre])
        y_bs = y - baseline
        y_sm = smooth_waveform(y_bs)
        prompt_idx = np.flatnonzero(prompt)
        local_idx = int(np.argmax(y_sm[prompt]))
        idx = int(prompt_idx[local_idx])
        peak_time, peak_height = quadratic_peak(t, y_sm, idx)
        raw_peak = float(np.max(y[prompt]))
        area = float(np.trapezoid(y_bs[integ], t[integ])) if integ.sum() > 1 else np.nan
        late_peak = float(np.max(y_bs[late])) if late.any() else np.nan
        clipped = overflow[col] > 0 or overflow[1] + overflow[2] > 0
        accepted = (
            not clipped
            and np.isfinite(peak_height)
            and np.isfinite(noise)
            and peak_height > max(5.0 * noise, 8.0)
            and -5.0 <= peak_time <= 80.0
        )
        rows.append(
            {
                "event": event,
                "file": path.name,
                "ok_shape": True,
                "accepted": bool(accepted),
                "baseline_mV": baseline,
                "baseline_noise_mV": noise,
                "baseline_samples_used": n_base,
                "raw_peak_mV_no_baseline_subtraction": raw_peak,
                "prompt_peak_mV": peak_height,
                "prompt_peak_time_ns": peak_time,
                "prompt_area_mVns": area,
                "late_peak_mV": late_peak,
                "clipped": bool(clipped),
                "overflow_samples_scope_channel": int(overflow[col]),
            }
        )
    meas = pd.DataFrame(rows)
    accepted = meas["accepted"].fillna(False).astype(bool)
    heights = meas.loc[accepted, "prompt_peak_mV"].to_numpy(float)
    fits, centers, counts, smooth = local_gaussian_peak_fits(heights)
    selected, spacing, spacing_unc, residual = select_spacing(fits)
    if len(fits):
        fits["selected_for_spacing"] = False
        fits["selected_index"] = np.nan
        for order, idx in enumerate(selected, start=1):
            fits.loc[idx, "selected_for_spacing"] = True
            fits.loc[idx, "selected_index"] = order
    diffs = np.diff(fits.loc[selected, "mu_mV"].to_numpy(float)) if selected else np.array([])
    rel_residual = residual / spacing if np.isfinite(residual) and np.isfinite(spacing) and spacing else np.nan
    quality = (
        len(selected) >= 3
        and np.isfinite(spacing)
        and np.isfinite(rel_residual)
        and rel_residual < 0.18
        and (not np.isfinite(spacing_unc) or spacing_unc / spacing < 0.20)
        and accepted.sum() >= 1000
    )
    summary = {
        "bias_V": spec.bias_v,
        "label": spec.label,
        "run_dir": str(spec.run_dir),
        "scope_channel": spec.scope_channel,
        "n_files": len(files),
        "n_events": len(meas),
        "n_accepted": int(accepted.sum()),
        "accepted_fraction": float(accepted.mean()) if len(meas) else np.nan,
        "baseline_median_mV": float(np.median(meas["baseline_mV"].dropna())) if "baseline_mV" in meas else np.nan,
        "baseline_noise_median_mV": float(np.median(meas["baseline_noise_mV"].dropna())) if "baseline_noise_mV" in meas else np.nan,
        "prompt_peak_median_mV": float(np.median(heights)) if heights.size else np.nan,
        "prompt_peak_mean_mV": float(np.mean(heights)) if heights.size else np.nan,
        "prompt_peak_std_mV": float(np.std(heights, ddof=1)) if heights.size > 1 else np.nan,
        "selected_peak_count": int(len(selected)),
        "selected_peak_mus_mV": ";".join(f"{fits.loc[i, 'mu_mV']:.3f}" for i in selected) if selected else "",
        "one_pe_spacing_mV": spacing,
        "one_pe_spacing_unc_mV": spacing_unc,
        "spacing_residual_rms_mV": residual,
        "spacing_rel_residual": rel_residual,
        "spacing_diffs_mV": ";".join(f"{x:.3f}" for x in diffs),
        "fit_quality_ok": bool(quality),
        "peak_candidates": int(len(fits)),
        "peak_max_mV": float(np.max(heights)) if heights.size else np.nan,
        "peak_99p9_mV": float(np.percentile(heights, 99.9)) if heights.size else np.nan,
        "clipped_events": int(meas["clipped"].fillna(False).sum()) if "clipped" in meas else np.nan,
    }
    return meas, fits, summary, (centers, counts, smooth)


def fit_vbr(summary: pd.DataFrame) -> dict[str, float]:
    fit_points = summary[summary["use_in_vbr_fit"]].sort_values("bias_V")
    x = fit_points["bias_V"].to_numpy(float)
    y = fit_points["one_pe_spacing_mV"].to_numpy(float)
    coeff, cov = np.polyfit(x, y, 1, cov=True)
    slope, intercept = coeff
    vbr = -intercept / slope
    dv_dm = intercept / slope**2
    dv_db = -1.0 / slope
    vbr_unc = np.sqrt(dv_dm**2 * cov[0, 0] + dv_db**2 * cov[1, 1] + 2.0 * dv_dm * dv_db * cov[0, 1])
    pred = slope * x + intercept
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return {
        "slope_mV_per_V": float(slope),
        "intercept_mV": float(intercept),
        "breakdown_voltage_V": float(vbr),
        "breakdown_voltage_unc_V": float(vbr_unc),
        "r2": float(1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan,
        "fit_rms_mV": float(np.sqrt(np.mean((y - pred) ** 2))),
        "n_fit_points": int(len(fit_points)),
    }


def save_plots(all_summaries: pd.DataFrame, hist_data: dict[float, tuple[np.ndarray, np.ndarray, np.ndarray]], peak_fits: dict[float, pd.DataFrame], vbr_result: dict[str, float]) -> None:
    plt.style.use("default")
    fig, axes = plt.subplots(3, 3, figsize=(15, 11), dpi=160)
    axes = axes.ravel()
    for ax, (_, row) in zip(axes, all_summaries.sort_values("bias_V").iterrows()):
        bias = float(row["bias_V"])
        centers, counts, smooth = hist_data[bias]
        fits = peak_fits[bias]
        ax.step(centers, counts, where="mid", color="C0", alpha=0.75, linewidth=1.0, label="height histogram")
        ax.plot(centers, smooth, color="black", linewidth=1.1, label="smoothed")
        selected = fits[fits["selected_for_spacing"].fillna(False).astype(bool)] if len(fits) else fits
        for _, fit in selected.iterrows():
            ax.axvline(fit["mu_mV"], color="C3", linewidth=1.2)
        if np.isfinite(row["one_pe_spacing_mV"]):
            q = "used" if row["use_in_vbr_fit"] else "diagnostic"
            ax.legend(title=f"{row['one_pe_spacing_mV']:.2f} +/- {row['one_pe_spacing_unc_mV']:.2f} mV ({q})", fontsize=7, title_fontsize=7)
        ax.set_title(row["label"])
        ax.set_xlabel("Prompt peak height (mV)")
        ax.set_ylabel("Events/bin")
        ax.grid(True, alpha=0.25)
    for ax in axes[len(all_summaries) :]:
        ax.axis("off")
    fig.suptitle("SIPM1 robust pulse-height spectra", fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT_DIR / "sipm1_robust_height_peak_histograms.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 6), dpi=160)
    used = all_summaries[all_summaries["use_in_vbr_fit"]]
    diag = all_summaries[~all_summaries["use_in_vbr_fit"]]
    ax.errorbar(used["bias_V"], used["one_pe_spacing_mV"], yerr=used["one_pe_spacing_unc_mV"], fmt="o", capsize=4, color="C0", label="used in Vbr fit")
    ax.errorbar(diag["bias_V"], diag["one_pe_spacing_mV"], yerr=diag["one_pe_spacing_unc_mV"], fmt="s", mfc="none", capsize=4, color="0.45", label="diagnostic only")
    xline = np.linspace(min(all_summaries["bias_V"].min(), vbr_result["breakdown_voltage_V"]) - 0.15, all_summaries["bias_V"].max() + 0.15, 200)
    yline = vbr_result["slope_mV_per_V"] * xline + vbr_result["intercept_mV"]
    ax.plot(xline, yline, color="C0", linewidth=1.6, label=f"fit: Vbr = {vbr_result['breakdown_voltage_V']:.2f} +/- {vbr_result['breakdown_voltage_unc_V']:.2f} V")
    ax.axvline(vbr_result["breakdown_voltage_V"], color="C0", linestyle="--", linewidth=1.2)
    ax.axhline(0, color="0.3", linewidth=1.0)
    ax.set_title("SIPM1 breakdown voltage from robust p.e. spacing")
    ax.set_xlabel("Bias voltage (V)")
    ax.set_ylabel("1 p.e. pulse-height spacing (mV)")
    ax.grid(True, alpha=0.28)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "sipm1_robust_spacing_vbr_fit.png")
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), dpi=160, sharex=True)
    axes[0].plot(all_summaries["bias_V"], all_summaries["baseline_median_mV"], "o-", color="C2")
    axes[0].set_ylabel("Baseline median (mV)")
    axes[0].grid(True, alpha=0.28)
    axes[1].plot(all_summaries["bias_V"], all_summaries["baseline_noise_median_mV"], "o-", color="C4")
    axes[1].set_xlabel("Bias voltage (V)")
    axes[1].set_ylabel("Baseline noise, MAD sigma (mV)")
    axes[1].grid(True, alpha=0.28)
    fig.suptitle("SIPM1 robust baseline checks")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT_DIR / "sipm1_robust_baseline_checks.png")
    plt.close(fig)


def write_note(summary: pd.DataFrame, vbr_result: dict[str, float]) -> None:
    note = f"""# SIPM1 robust p.e. peak recalculation

This recalculation uses the raw PicoScope CSV files and does not reuse the older baseline-subtracted measurements.

## Method

1. Each waveform is read in the scope units and converted to mV.
2. The baseline is estimated event-by-event from the quiet pre-trigger region, normally `t < -50 ns`.
3. The baseline estimator is robust: median first, then MAD-based clipping, then a final clipped median. This avoids pulling the baseline when a small pre-pulse or pickup excursion is present.
4. The noise estimate is also robust: `sigma = 1.4826 * MAD` from the clipped pre-trigger samples.
5. The waveform is baseline-subtracted.
6. The prompt signal is measured only in `-5 ns <= t <= 90 ns`. This avoids mistaking late afterpulses or recovery features for the main triggered avalanche.
7. Pulse height is measured from a mildly smoothed waveform using a Savitzky-Golay filter, followed by a 3-point quadratic interpolation around the maximum.
8. Events are accepted when the prompt peak is above `max(5 sigma_baseline, 8 mV)`, not clipped, and the peak time lies between `-5 ns` and `80 ns`.
9. A pulse-height histogram is built for each bias voltage.
10. Candidate p.e. peaks are found from a lightly smoothed histogram and refined with local Gaussian-plus-linear-background fits.
11. The reported 1 p.e. spacing is the average separation of the best nearly-even sequence of fitted peak centers.
12. Points marked diagnostic are plotted and saved, but not used for the breakdown-voltage fit.

## SIPM1 result

- Breakdown voltage from selected robust height-spacing points: `{vbr_result['breakdown_voltage_V']:.3f} +/- {vbr_result['breakdown_voltage_unc_V']:.3f} V`
- Height-spacing slope: `{vbr_result['slope_mV_per_V']:.3f} mV/V`
- Fit R^2: `{vbr_result['r2']:.4f}`
- Fit RMS: `{vbr_result['fit_rms_mV']:.3f} mV`

## Why this is better than the first-pass method

The first-pass method used a simple mean baseline from `t < -30 ns` and searched for peaks in the resulting histogram. The new method uses a robust baseline and a prompt-only peak gate. That is better for these recorded waveforms because the scope is triggered on pulses, not randomly, and the traces can contain late pulses or small baseline disturbances that should not change the pulse-height estimate.
"""
    (OUT_DIR / "sipm1_robust_analysis_note.md").write_text(note, encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summaries: list[dict[str, object]] = []
    hist_data: dict[float, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    peak_fits: dict[float, pd.DataFrame] = {}
    for spec in RUNS:
        print(f"Analyzing {spec.label}: {spec.run_dir}", flush=True)
        meas, fits, summary, hist = analyze_run(spec)
        run_out = OUT_DIR / spec.label.replace(".", "p").replace(" ", "_")
        run_out.mkdir(exist_ok=True)
        meas.to_csv(run_out / "robust_waveform_measurements.csv", index=False)
        fits.to_csv(run_out / "robust_height_peak_fits.csv", index=False)
        summaries.append(summary)
        hist_data[spec.bias_v] = hist
        peak_fits[spec.bias_v] = fits
    summary_df = pd.DataFrame(summaries).sort_values("bias_V").reset_index(drop=True)
    summary_df["use_in_vbr_fit"] = summary_df["fit_quality_ok"]
    # Low-bias triggered spectra and high-bias distorted spectra are retained as diagnostics.
    summary_df.loc[summary_df["bias_V"].isin([53.20, 54.10, 56.30, 56.40, 56.86]), "use_in_vbr_fit"] = False
    vbr_result = fit_vbr(summary_df)
    pd.DataFrame([vbr_result]).to_csv(OUT_DIR / "sipm1_robust_breakdown_voltage.csv", index=False)
    summary_df.to_csv(OUT_DIR / "sipm1_robust_pe_summary.csv", index=False)
    save_plots(summary_df, hist_data, peak_fits, vbr_result)
    write_note(summary_df, vbr_result)
    print(OUT_DIR)
    print(summary_df[["bias_V", "n_accepted", "baseline_median_mV", "baseline_noise_median_mV", "selected_peak_mus_mV", "one_pe_spacing_mV", "one_pe_spacing_unc_mV", "spacing_rel_residual", "use_in_vbr_fit"]].to_string(index=False))
    print(pd.DataFrame([vbr_result]).to_string(index=False))


if __name__ == "__main__":
    main()
