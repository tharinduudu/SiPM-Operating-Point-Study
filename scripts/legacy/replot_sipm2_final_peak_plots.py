#!/usr/bin/env python3
from __future__ import annotations

from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
from scipy.signal import find_peaks


ROOT = Path("/home/muon/brDownVstudy/robust_sipm2_area_analysis_20260901")
OUT_HEIGHT_SPECTRA = ROOT / "sipm2_height_spectra_low_to_high_high_bias_distance.png"
OUT_HEIGHT_SUMMARY = ROOT / "sipm2_height_peak_summary_high_bias_distance.csv"
OUT_HEIGHT_FIT = ROOT / "sipm2_breakdown_fit_excluding_53p20_54p10.png"
OUT_HEIGHT_FIT_CSV = ROOT / "sipm2_breakdown_fit_excluding_53p20_54p10.csv"
OUT_HEIGHT_POINTS_CSV = ROOT / "sipm2_breakdown_points_excluding_53p20_54p10.csv"
OUT_AREA_FIT = ROOT / "sipm2_area_spacing_vs_bias_excluding_53p20_54p10.png"
OUT_AREA_FIT_CSV = ROOT / "sipm2_area_spacing_fit_excluding_53p20_54p10.csv"
OUT_AREA_POINTS_CSV = ROOT / "sipm2_area_spacing_points_excluding_53p20_54p10.csv"


def run_folder_name(label: str) -> str:
    return label.replace(".", "p").replace(" ", "_")


def fd_bin_width(values: np.ndarray, lower: float, upper: float) -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    x = x[(x >= lower) & (x <= upper)]
    if x.size < 2:
        return max((upper - lower) / 80.0, 1.0)
    q25, q75 = np.percentile(x, [25, 75])
    width = 2.0 * (q75 - q25) / np.cbrt(x.size)
    if not np.isfinite(width) or width <= 0:
        width = (upper - lower) / 100.0
    return float(np.clip(width, 0.7, 4.0))


def histogram_height(values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    lower, upper = np.percentile(x, [0.2, 99.4])
    pad = max(5.0, 0.05 * (upper - lower))
    lower -= pad
    upper += pad
    bin_width = fd_bin_width(x, lower, upper)
    edges = np.arange(
        np.floor(lower / bin_width) * bin_width,
        np.ceil(upper / bin_width) * bin_width + bin_width,
        bin_width,
    )
    counts, edges = np.histogram(x, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    smooth = gaussian_filter1d(counts.astype(float), sigma=max(1.0, 1.2 / bin_width))
    return centers, counts, smooth, bin_width


def gaussian_linear_bg(x: np.ndarray, amp: float, mu: float, sigma: float, bg0: float, bg1: float) -> np.ndarray:
    return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2) + bg0 + bg1 * (x - mu)


def fit_height_peaks(
    heights: np.ndarray, min_peak_distance_mV: float
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray, float]:
    centers, counts, smooth, bin_width = histogram_height(heights)
    empty = pd.DataFrame(columns=["mu_mV", "mu_err_mV", "sigma_mV", "amp", "prominence", "selected_for_spacing"])
    if counts.size < 10 or counts.max() <= 0:
        return empty, centers, counts, smooth, bin_width

    prominence = max(8.0, 0.022 * float(np.max(smooth)))
    distance = max(3, int(round(min_peak_distance_mV / bin_width)))
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
                bounds=([0.0, mu0 - half_width, 0.25, 0.0, -np.inf], [np.inf, mu0 + half_width, 35.0, np.inf, np.inf]),
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
    if not rows:
        return empty, centers, counts, smooth, bin_width
    return pd.DataFrame(rows).sort_values("mu_mV").reset_index(drop=True), centers, counts, smooth, bin_width


def select_spacing(fits: pd.DataFrame, expected_spacing_mV: float | None = None) -> tuple[list[int], float, float, float]:
    if len(fits) < 2:
        return [], np.nan, np.nan, np.nan
    mus = fits["mu_mV"].to_numpy(float)
    best: tuple[float, tuple[int, ...], float, float] | None = None
    for count in range(min(6, len(mus)), 1, -1):
        for combo in combinations(range(len(mus)), count):
            vals = mus[list(combo)]
            diffs = np.diff(vals)
            spacing = float(np.median(diffs))
            if spacing < 8.0 or spacing > 140.0:
                continue
            residual = float(np.sqrt(np.mean((diffs - spacing) ** 2)))
            coverage = vals[-1] - vals[0]
            expected_penalty = 0.0
            if expected_spacing_mV and np.isfinite(expected_spacing_mV):
                if np.any(diffs < 0.70 * expected_spacing_mV) or np.any(diffs > 1.40 * expected_spacing_mV):
                    continue
                expected_penalty = 0.5 * abs(spacing - expected_spacing_mV) / expected_spacing_mV
            score = residual / spacing + expected_penalty - 0.028 * count - 0.0005 * coverage
            if best is None or score < best[0]:
                best = (score, combo, spacing, residual)
        if best is not None and not expected_spacing_mV:
            break
    if best is None:
        return [], np.nan, np.nan, np.nan
    _, combo, _, residual = best
    selected = list(combo)
    selected_mus = mus[selected]
    diffs = np.diff(selected_mus)
    spacing = float(np.mean(diffs))
    if len(diffs) > 1:
        spacing_unc = float(np.std(diffs, ddof=1) / np.sqrt(len(diffs)))
    else:
        mu_err = fits.loc[selected, "mu_err_mV"].to_numpy(float)
        spacing_unc = float(np.sqrt(np.nansum(mu_err**2))) if np.isfinite(mu_err).all() else np.nan
    return selected, spacing, spacing_unc, residual


def linear_fit(points: pd.DataFrame, x_col: str, y_col: str, fit_kind: str, vbr_name: str) -> dict[str, float | str | int]:
    p = points[np.isfinite(points[y_col])].sort_values(x_col)
    x = p[x_col].to_numpy(float)
    y = p[y_col].to_numpy(float)
    coeff, cov = np.polyfit(x, y, 1, cov=True)
    slope, intercept = coeff
    vbr = -intercept / slope
    dv_dm = intercept / slope**2
    dv_db = -1.0 / slope
    vbr_unc = float(np.sqrt(dv_dm**2 * cov[0, 0] + dv_db**2 * cov[1, 1] + 2.0 * dv_dm * dv_db * cov[0, 1]))
    pred = slope * x + intercept
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return {
        "fit_kind": fit_kind,
        "n_fit_points": int(len(p)),
        "slope": float(slope),
        "intercept": float(intercept),
        vbr_name: float(vbr),
        f"{vbr_name}_unc": vbr_unc,
        "r2": float(1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan,
        "fit_rms": float(np.sqrt(np.mean((y - pred) ** 2))),
    }


def fit_seed_line(rows: list[dict[str, object]]) -> tuple[float, float]:
    seed = pd.DataFrame(rows)
    seed = seed[seed["bias_V"].isin([54.90, 55.59, 55.65]) & np.isfinite(seed["one_pe_spacing_mV"])]
    if len(seed) < 2:
        return np.nan, np.nan
    slope, intercept = np.polyfit(seed["bias_V"].to_numpy(float), seed["one_pe_spacing_mV"].to_numpy(float), 1)
    return float(slope), float(intercept)


def height_summary_pass(summary: pd.DataFrame, seed_slope: float = np.nan, seed_intercept: float = np.nan) -> tuple[pd.DataFrame, dict[float, tuple[np.ndarray, np.ndarray, np.ndarray]], dict[float, pd.DataFrame]]:
    rows: list[dict[str, object]] = []
    hist_data: dict[float, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    peak_fits: dict[float, pd.DataFrame] = {}
    for _, row in summary.sort_values("bias_V").iterrows():
        bias = float(row["bias_V"])
        label = str(row["label"])
        run_dir = ROOT / run_folder_name(label)
        meas = pd.read_csv(run_dir / "robust_waveform_measurements.csv")
        accepted = meas["accepted"].fillna(False).astype(bool)
        heights = meas.loc[accepted, "prompt_peak_mV"].dropna().to_numpy(float)

        expected_gap = seed_slope * bias + seed_intercept if np.isfinite(seed_slope) and np.isfinite(seed_intercept) else np.nan
        high_bias_distance = bool(bias >= 56.30 and np.isfinite(expected_gap))
        min_peak_distance = max(6.0, 0.50 * expected_gap) if high_bias_distance else 6.0
        expected_for_selection = expected_gap if high_bias_distance else None

        fits, centers, counts, smooth, bin_width = fit_height_peaks(heights, min_peak_distance)
        selected, spacing, spacing_unc, residual = select_spacing(fits, expected_for_selection)
        if len(fits):
            fits["selected_for_spacing"] = False
            fits["selected_index"] = np.nan
            for order, idx in enumerate(selected, start=1):
                fits.loc[idx, "selected_for_spacing"] = True
                fits.loc[idx, "selected_index"] = order
        selected_mus = fits.loc[selected, "mu_mV"].to_numpy(float) if selected else np.array([])
        diffs = np.diff(selected_mus) if selected else np.array([])
        rows.append(
            {
                "bias_V": bias,
                "label": label,
                "n_accepted": int(len(heights)),
                "expected_gap_from_seed_fit_mV": expected_gap,
                "min_peak_distance_mV": min_peak_distance,
                "bin_width_mV": bin_width,
                "peak_candidates": len(fits),
                "selected_peak_count": len(selected),
                "selected_peak_mus_mV": ";".join(f"{x:.3f}" for x in selected_mus),
                "one_pe_spacing_mV": spacing,
                "one_pe_spacing_unc_mV": spacing_unc,
                "spacing_residual_rms_mV": residual,
                "spacing_diffs_mV": ";".join(f"{x:.3f}" for x in diffs),
                "high_bias_distance_applied": high_bias_distance,
                "used_in_final_vbr_fit": bool(bias not in (53.20, 54.10) and np.isfinite(spacing)),
            }
        )
        fits.to_csv(run_dir / "height_peak_fits_high_bias_distance.csv", index=False)
        hist_data[bias] = (centers, counts, smooth)
        peak_fits[bias] = fits
    return pd.DataFrame(rows).sort_values("bias_V").reset_index(drop=True), hist_data, peak_fits


def plot_height_spectra(height_points: pd.DataFrame, hist_data: dict[float, tuple[np.ndarray, np.ndarray, np.ndarray]], peak_fits: dict[float, pd.DataFrame]) -> None:
    fig, axes = plt.subplots(3, 3, figsize=(16, 13), dpi=160)
    axes_flat = np.asarray(axes).ravel()
    for ax, (_, row) in zip(axes_flat, height_points.sort_values("bias_V").iterrows()):
        bias = float(row["bias_V"])
        centers, counts, smooth = hist_data[bias]
        fits = peak_fits[bias]
        ax.step(centers, counts, where="mid", color="C0", alpha=0.75, linewidth=1.0, label="height histogram")
        ax.plot(centers, smooth, color="black", linewidth=1.1, label="smoothed")
        selected = fits[fits["selected_for_spacing"].fillna(False).astype(bool)] if len(fits) else pd.DataFrame()
        for _, peak in selected.iterrows():
            ax.axvline(float(peak["mu_mV"]), color="C3", linewidth=1.2)
        status = "used" if bool(row["used_in_final_vbr_fit"]) else "diagnostic"
        if np.isfinite(row["one_pe_spacing_mV"]):
            unc = row["one_pe_spacing_unc_mV"]
            title = f"{row['one_pe_spacing_mV']:.2f}"
            if np.isfinite(unc) and unc <= 0.5 * abs(row["one_pe_spacing_mV"]):
                title += f" +/- {unc:.2f}"
            title += f" mV ({status})"
        else:
            title = status
        ax.legend(title=title, fontsize=7, title_fontsize=7, loc="upper right")
        ax.set_title(str(row["label"]))
        ax.set_xlabel("Prompt peak height (mV)")
        ax.set_ylabel("Events/bin")
        ax.grid(True, alpha=0.25)
    for ax in axes_flat[len(height_points) :]:
        ax.axis("off")
    fig.suptitle("SIPM2 robust pulse-height spectra", fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT_HEIGHT_SPECTRA)
    plt.close(fig)


def plot_height_fit(height_points: pd.DataFrame, fit: dict[str, float | str | int]) -> None:
    points = height_points[height_points["used_in_final_vbr_fit"].astype(bool)].sort_values("bias_V")
    x = points["bias_V"].to_numpy(float)
    y = points["one_pe_spacing_mV"].to_numpy(float)
    yerr = points["one_pe_spacing_unc_mV"].to_numpy(float)
    yerr[~np.isfinite(yerr) | (yerr > 0.5 * np.abs(y))] = np.nan

    fig, ax = plt.subplots(figsize=(9.2, 5.9), dpi=180)
    ax.errorbar(x, y, yerr=yerr, fmt="o", capsize=4, color="C0", markersize=6, label="height p.e. gap points")
    for _, row in points.iterrows():
        ax.annotate(f"{row['bias_V']:.2f} V", (row["bias_V"], row["one_pe_spacing_mV"]), xytext=(6, 6), textcoords="offset points", fontsize=8)
    xline = np.linspace(min(float(fit["breakdown_voltage_V"]), float(points["bias_V"].min())) - 0.35, float(points["bias_V"].max()) + 0.25, 300)
    yline = float(fit["slope"]) * xline + float(fit["intercept"])
    ax.plot(xline, yline, color="C3", linewidth=1.8, label=f"fit: Vbr = {float(fit['breakdown_voltage_V']):.2f} +/- {float(fit['breakdown_voltage_V_unc']):.2f} V")
    ax.axhline(0, color="0.35", linewidth=1.0)
    ax.axvline(float(fit["breakdown_voltage_V"]), color="C3", linestyle="--", linewidth=1.2)
    info = (
        f"slope = {float(fit['slope']):.2f} mV/V\n"
        f"R^2 = {float(fit['r2']):.3f}\n"
        f"fit RMS = {float(fit['fit_rms']):.2f} mV\n"
        f"n = {int(fit['n_fit_points'])}"
    )
    ax.text(0.03, 0.96, info, transform=ax.transAxes, va="top", ha="left", fontsize=9, bbox=dict(facecolor="white", edgecolor="0.75", alpha=0.9))
    ax.set_title("SIPM2 p.e. peak spacing vs bias")
    ax.set_xlabel("Bias voltage (V)")
    ax.set_ylabel("1 p.e. height spacing (mV)")
    ymax = float(np.nanmax(y) * 1.18)
    ax.set_ylim(-0.08 * ymax, ymax)
    ax.set_xlim(min(float(fit["breakdown_voltage_V"]), float(points["bias_V"].min())) - 0.45, float(points["bias_V"].max()) + 0.65)
    ax.grid(True, alpha=0.28)
    ax.legend(fontsize=9, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT_HEIGHT_FIT)
    plt.close(fig)


def plot_area_fit() -> dict[str, float | str | int] | None:
    path = ROOT / "sipm2_area_spacing_points.csv"
    if not path.exists():
        return None
    area = pd.read_csv(path)
    points = area[(area["bias_V"] != 54.10) & np.isfinite(area["area_spacing_mVns"])].copy()
    points.to_csv(OUT_AREA_POINTS_CSV, index=False)
    fit = linear_fit(points, "bias_V", "area_spacing_mVns", "exclude_53p20_54p10_area_unweighted", "breakdown_voltage_from_area_V")
    pd.DataFrame([fit]).to_csv(OUT_AREA_FIT_CSV, index=False)

    x = points["bias_V"].to_numpy(float)
    y = points["area_spacing_mVns"].to_numpy(float)
    yerr = points["area_spacing_unc_mVns"].to_numpy(float)
    yerr[~np.isfinite(yerr) | (yerr > 0.5 * np.abs(y))] = np.nan
    fig, ax = plt.subplots(figsize=(9.2, 5.9), dpi=180)
    ax.errorbar(x, y, yerr=yerr, fmt="o", capsize=4, color="C0", markersize=6, label="area p.e. gap points")
    for _, row in points.iterrows():
        ax.annotate(f"{row['bias_V']:.2f} V", (row["bias_V"], row["area_spacing_mVns"]), xytext=(6, 6), textcoords="offset points", fontsize=8)
    xline = np.linspace(min(float(fit["breakdown_voltage_from_area_V"]), float(points["bias_V"].min())) - 0.35, float(points["bias_V"].max()) + 0.25, 300)
    yline = float(fit["slope"]) * xline + float(fit["intercept"])
    ax.plot(
        xline,
        yline,
        color="C3",
        linewidth=1.8,
        label=f"fit: Vbr = {float(fit['breakdown_voltage_from_area_V']):.2f} +/- {float(fit['breakdown_voltage_from_area_V_unc']):.2f} V",
    )
    ax.axhline(0, color="0.35", linewidth=1.0)
    ax.axvline(float(fit["breakdown_voltage_from_area_V"]), color="C3", linestyle="--", linewidth=1.2)
    info = (
        f"slope = {float(fit['slope']):.0f} mV ns/V\n"
        f"R^2 = {float(fit['r2']):.3f}\n"
        f"fit RMS = {float(fit['fit_rms']):.0f} mV ns\n"
        f"n = {int(fit['n_fit_points'])}"
    )
    ax.text(0.03, 0.96, info, transform=ax.transAxes, va="top", ha="left", fontsize=9, bbox=dict(facecolor="white", edgecolor="0.75", alpha=0.9))
    ax.set_title("SIPM2 p.e. area spacing vs bias")
    ax.set_xlabel("Bias voltage (V)")
    ax.set_ylabel("1 p.e. area spacing (mV ns)")
    ymax = float(np.nanmax(y) * 1.18)
    ax.set_ylim(-0.08 * ymax, ymax)
    ax.set_xlim(min(float(fit["breakdown_voltage_from_area_V"]), float(points["bias_V"].min())) - 0.45, float(points["bias_V"].max()) + 0.65)
    ax.grid(True, alpha=0.28)
    ax.legend(fontsize=9, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT_AREA_FIT)
    plt.close(fig)
    return fit


def main() -> None:
    summary = pd.read_csv(ROOT / "sipm2_robust_area_summary.csv").sort_values("bias_V").reset_index(drop=True)

    first_pass, _, _ = height_summary_pass(summary)
    seed_slope, seed_intercept = fit_seed_line(first_pass.to_dict("records"))
    height_points, hist_data, peak_fits = height_summary_pass(summary, seed_slope, seed_intercept)
    height_points.to_csv(OUT_HEIGHT_SUMMARY, index=False)

    final_height_points = height_points[height_points["used_in_final_vbr_fit"].astype(bool)].copy()
    height_fit = linear_fit(final_height_points, "bias_V", "one_pe_spacing_mV", "exclude_53p20_54p10_height_unweighted", "breakdown_voltage_V")
    pd.DataFrame([height_fit]).to_csv(OUT_HEIGHT_FIT_CSV, index=False)
    final_height_points.to_csv(OUT_HEIGHT_POINTS_CSV, index=False)

    plot_height_spectra(height_points, hist_data, peak_fits)
    plot_height_fit(height_points, height_fit)
    area_fit = plot_area_fit()

    print(OUT_HEIGHT_SPECTRA)
    print(OUT_HEIGHT_FIT)
    if area_fit is not None:
        print(OUT_AREA_FIT)
    print(height_points[["bias_V", "selected_peak_mus_mV", "one_pe_spacing_mV", "one_pe_spacing_unc_mV", "used_in_final_vbr_fit"]].to_string(index=False))
    print(pd.DataFrame([height_fit]).to_string(index=False))
    if area_fit is not None:
        print(pd.DataFrame([area_fit]).to_string(index=False))


if __name__ == "__main__":
    main()
