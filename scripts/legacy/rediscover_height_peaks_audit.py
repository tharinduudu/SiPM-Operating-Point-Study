#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks


ROOT = Path("/Users/tharinduhettiarachchi/Desktop/Design/new415")
OUT_DIR = ROOT / "height_gap_linear_calibration_20260902" / "peak_rediscovery_audit_20260903"


@dataclass(frozen=True)
class RunSpec:
    sipm: str
    bias_v: float
    group: str
    path: Path


FIT_RUNS = [
    RunSpec("SIPM1", 53.90, "fit", ROOT / "brDownVstudy_53p9V_sipm1_ch0_analysis_20260827/analysis_sipm1_ch0_53p9V/pulse_measurements.csv"),
    RunSpec("SIPM1", 54.90, "fit", ROOT / "brDownVstudy_54p9V_analysis_20260827/waveform_plots_54p9V_20260827/sipm1_ch0/pulse_measurements.csv"),
    RunSpec("SIPM1", 55.59, "fit", ROOT / "brDownVstudy_55p59V_analysis_20260827/20260827_55p59V_sipm1_ch0_10000/analysis_sipm1_ch0_55p59V/pulse_measurements.csv"),
    RunSpec("SIPM1", 55.65, "fit", ROOT / "brDownVstudy_55p65V_analysis_20260827/waveform_plots_55p65V_20260827/sipm1_ch0/pulse_measurements.csv"),
    RunSpec("SIPM1", 56.30, "fit", ROOT / "brDownVstudy_56p3V_both_analysis_20260831/waveform_plots_56p3V_20260831_both20k/sipm1_ch0/pulse_measurements.csv"),
    RunSpec("SIPM1", 56.40, "fit", ROOT / "brDownVstudy_56p4V_analysis_20260827/waveform_plots_56p4V_20260827/sipm1_ch0/pulse_measurements.csv"),
    RunSpec("SIPM1", 56.86, "fit", ROOT / "brDownVstudy_56p86V_20000_analysis_20260831/analysis_sipm1_ch0_56p86V_20000/pulse_measurements.csv"),
    RunSpec("SIPM2", 54.90, "fit", ROOT / "brDownVstudy_54p9V_analysis_20260827/waveform_plots_54p9V_20260827/sipm2_ch2/pulse_measurements.csv"),
    RunSpec("SIPM2", 55.59, "fit", ROOT / "brDownVstudy_55p59V_analysis_20260827/20260827_55p59V_sipm2_ch2_10000/analysis_sipm2_ch2_55p59V/pulse_measurements.csv"),
    RunSpec("SIPM2", 55.65, "fit", ROOT / "brDownVstudy_55p65V_analysis_20260827/waveform_plots_55p65V_20260827/sipm2_ch2/pulse_measurements.csv"),
    RunSpec("SIPM2", 56.30, "fit", ROOT / "brDownVstudy_56p3V_both_analysis_20260831/waveform_plots_56p3V_20260831_both20k/sipm2_ch2/pulse_measurements.csv"),
    RunSpec("SIPM2", 56.40, "fit", ROOT / "brDownVstudy_56p4V_analysis_20260827/waveform_plots_56p4V_20260827/sipm2_ch2/pulse_measurements.csv"),
    RunSpec("SIPM2", 56.86, "fit", ROOT / "brDownVstudy_56p86V_20kSIPM1_analysis_20260831/waveform_plots_56p86V_20260831_20kSIPM1/sipm2_ch2/pulse_measurements.csv"),
]

CHECK_RUNS = [
    RunSpec("SIPM1", 56.45, "channel-switched check", ROOT / "sep2_swapped_channel_waveform_check_20260902/56p45_sipm1_ch2_20000-0003/analysis_sipm1_ch2_56p45V_20000/pulse_measurements.csv"),
    RunSpec("SIPM1", 57.867, "channel-switched check", ROOT / "sep2_swapped_channel_waveform_check_20260902/20260831-57p867_sipm1_ch2_20000/analysis_sipm1_ch2_57p867V_20000/pulse_measurements.csv"),
    RunSpec("SIPM2", 56.45, "channel-switched check", ROOT / "sep2_swapped_channel_waveform_check_20260902/56p45_sipm2_ch0_20000-0004/analysis_sipm2_ch0_56p45V_20000/pulse_measurements.csv"),
    RunSpec("SIPM2", 57.867, "channel-switched check", ROOT / "sep2_swapped_channel_waveform_check_20260902/20260831-57p867_sipm2_ch0_20000-0002/analysis_sipm2_ch0_57p867V_20000/pulse_measurements.csv"),
]

CHANNEL_STYLE = {
    "A": {"color": "C0", "label": "detector CH0 / scope A"},
    "B": {"color": "C1", "label": "detector CH2 / scope B"},
}


def active_height_values(path: Path) -> tuple[np.ndarray, str, int]:
    data = pd.read_csv(path)
    channel = str(data["active_scope_channel"].dropna().iloc[0])
    accepted = data["accepted"].fillna(False).astype(bool)
    values = pd.to_numeric(data.loc[accepted, f"{channel}_peak_mV"], errors="coerce").dropna().to_numpy(float)
    return values[np.isfinite(values)], channel, int(accepted.sum())


def histogram(values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    bin_width = 3.0
    lo = 50.0
    hi = max(450.0, float(np.nanpercentile(values, 99.7)) + 12.0)
    hi = min(520.0, np.ceil(hi / bin_width) * bin_width)
    edges = np.arange(lo, hi + bin_width, bin_width)
    counts, edges = np.histogram(values, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    smooth = gaussian_filter1d(counts.astype(float), sigma=1.15)
    return centers, counts, smooth, bin_width


def parabolic_peak(center: np.ndarray, smooth: np.ndarray, idx: int, bin_width: float) -> float:
    if idx <= 0 or idx >= len(smooth) - 1:
        return float(center[idx])
    left = float(smooth[idx - 1])
    mid = float(smooth[idx])
    right = float(smooth[idx + 1])
    denom = left - 2.0 * mid + right
    if abs(denom) < 1e-12:
        return float(center[idx])
    delta_bins = 0.5 * (left - right) / denom
    delta_bins = float(np.clip(delta_bins, -0.5, 0.5))
    return float(center[idx] + delta_bins * bin_width)


def discover_candidates(centers: np.ndarray, smooth: np.ndarray, bin_width: float) -> pd.DataFrame:
    indexes, props = find_peaks(
        smooth,
        prominence=max(7.0, 0.004 * float(np.nanmax(smooth))),
        distance=max(2, int(round(6.0 / bin_width))),
    )
    rows = []
    for idx, prominence in zip(indexes, props["prominences"]):
        if not (55.0 <= centers[idx] <= 450.0):
            continue
        rows.append(
            {
                "mu_mV": parabolic_peak(centers, smooth, int(idx), bin_width),
                "center_mV": float(centers[idx]),
                "smooth_height": float(smooth[idx]),
                "prominence": float(prominence),
                "selected": False,
                "dropped_as_doublet": False,
                "dropped_as_leading_shoulder": False,
            }
        )
    return pd.DataFrame(rows).sort_values("mu_mV").reset_index(drop=True)


def collapse_doublets(candidates: pd.DataFrame) -> pd.DataFrame:
    if len(candidates) < 3:
        return candidates
    kept = candidates.copy()
    diffs = np.diff(kept["mu_mV"].to_numpy(float))
    large_diffs = diffs[diffs > 20.0]
    if large_diffs.size == 0:
        return kept
    typical_gap = float(np.median(large_diffs))
    close_limit = 0.45 * typical_gap
    drop_indexes: set[int] = set()
    i = 0
    while i < len(kept) - 1:
        if i in drop_indexes:
            i += 1
            continue
        j = i + 1
        if float(kept.loc[j, "mu_mV"]) - float(kept.loc[i, "mu_mV"]) < close_limit:
            drop = i if float(kept.loc[i, "smooth_height"]) < float(kept.loc[j, "smooth_height"]) else j
            drop_indexes.add(drop)
            kept.loc[drop, "dropped_as_doublet"] = True
            i += 2
        else:
            i += 1
    return kept


def select_pe_sequence(candidates: pd.DataFrame) -> tuple[pd.DataFrame, float, float, float]:
    marked = collapse_doublets(candidates)
    available = marked[~marked["dropped_as_doublet"]].copy().reset_index()
    if len(available) >= 2 and float(available.loc[1, "smooth_height"]) > float(available.loc[0, "smooth_height"]):
        marked.loc[int(available.loc[0, "index"]), "dropped_as_leading_shoulder"] = True
        available = available.iloc[1:].reset_index(drop=True)
    if len(available) < 2:
        return marked, np.nan, np.nan, np.nan

    selected_original_indexes = available["index"].to_numpy(int)
    marked.loc[selected_original_indexes, "selected"] = True
    selected = marked[marked["selected"]].sort_values("mu_mV")
    positions = selected["mu_mV"].to_numpy(float)
    peak_index = np.arange(len(positions), dtype=float)
    if len(positions) == 2:
        spacing = float(positions[1] - positions[0])
        spacing_unc = float(np.sqrt(2.0) * 3.0 / np.sqrt(12.0))
        residual_rms = 0.0
    else:
        coeff, cov = np.polyfit(peak_index, positions, 1, cov=True)
        spacing = float(coeff[0])
        residuals = positions - np.polyval(coeff, peak_index)
        residual_rms = float(np.sqrt(np.mean(residuals**2)))
        slope_unc = float(np.sqrt(cov[0, 0])) if np.isfinite(cov[0, 0]) else np.nan
        bin_unc = float((3.0 / np.sqrt(12.0)) / np.sqrt(np.sum((peak_index - peak_index.mean()) ** 2)))
        spacing_unc = float(np.hypot(slope_unc, bin_unc)) if np.isfinite(slope_unc) else bin_unc
    return marked, spacing, spacing_unc, residual_rms


def fit_vbr(points: pd.DataFrame) -> dict[str, float]:
    x = points["bias_V"].to_numpy(float)
    y = points["height_gap_mV"].to_numpy(float)
    coeff, cov = np.polyfit(x, y, 1, cov=True)
    slope, intercept = coeff
    vbr = -intercept / slope
    dv_dm = intercept / slope**2
    dv_db = -1.0 / slope
    vbr_unc = float(
        np.sqrt(
            dv_dm**2 * cov[0, 0]
            + dv_db**2 * cov[1, 1]
            + 2.0 * dv_dm * dv_db * cov[0, 1]
        )
    )
    pred = slope * x + intercept
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return {
        "n_fit_points": int(len(points)),
        "slope_mV_per_V": float(slope),
        "intercept_mV": float(intercept),
        "breakdown_voltage_V": float(vbr),
        "breakdown_voltage_unc_V": vbr_unc,
        "r2": float(1.0 - ss_res / ss_tot),
        "fit_rms_mV": float(np.sqrt(np.mean((y - pred) ** 2))),
    }


def analyze_run(run: RunSpec) -> tuple[dict[str, object], pd.DataFrame, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    values, channel, n_accepted = active_height_values(run.path)
    centers, counts, smooth, bin_width = histogram(values)
    candidates = discover_candidates(centers, smooth, bin_width)
    marked, spacing, spacing_unc, residual_rms = select_pe_sequence(candidates)
    selected = marked[marked["selected"]].sort_values("mu_mV")
    ignored = marked[marked["dropped_as_leading_shoulder"]].sort_values("mu_mV")
    doublets = marked[marked["dropped_as_doublet"]].sort_values("mu_mV")
    summary = {
        "sipm": run.sipm,
        "bias_V": run.bias_v,
        "group": run.group,
        "scope_channel": channel,
        "n_accepted": n_accepted,
        "height_gap_mV": spacing,
        "height_gap_unc_mV": spacing_unc,
        "sequence_residual_rms_mV": residual_rms,
        "selected_peak_count": int(len(selected)),
        "selected_height_peak_mV": ";".join(f"{v:.3f}" for v in selected["mu_mV"].to_numpy(float)),
        "candidate_height_peak_mV": ";".join(f"{v:.3f}" for v in marked["mu_mV"].to_numpy(float)),
        "dropped_doublet_peak_mV": ";".join(f"{v:.3f}" for v in doublets["mu_mV"].to_numpy(float)),
        "dropped_leading_peak_mV": ";".join(f"{v:.3f}" for v in ignored["mu_mV"].to_numpy(float)),
        "source": str(run.path),
    }
    marked.insert(0, "sipm", run.sipm)
    marked.insert(1, "bias_V", run.bias_v)
    marked.insert(2, "group", run.group)
    marked.insert(3, "scope_channel", channel)
    return summary, marked, (centers, counts, smooth)


def plot_histograms(summaries: pd.DataFrame, hist_data: dict[tuple[str, float, str], tuple[np.ndarray, np.ndarray, np.ndarray]]) -> list[Path]:
    paths: list[Path] = []
    for sipm in ["SIPM1", "SIPM2"]:
        rows = summaries[(summaries["sipm"].eq(sipm)) & (summaries["group"].eq("fit"))].sort_values("bias_V")
        fig, axes = plt.subplots(3, 3, figsize=(15.5, 11.0), dpi=170)
        axes = axes.ravel()
        for ax in axes[len(rows) :]:
            ax.axis("off")
        for ax, (_, row) in zip(axes, rows.iterrows()):
            centers, counts, smooth = hist_data[(row["sipm"], row["bias_V"], row["group"])]
            ax.step(centers, counts, where="mid", color="C0" if sipm == "SIPM1" else "C1", alpha=0.65, linewidth=1.0)
            ax.plot(centers, smooth, color="black", linewidth=1.2, alpha=0.86)
            for peak in str(row["selected_height_peak_mV"]).split(";"):
                if peak:
                    ax.axvline(float(peak), color="C3", linewidth=1.1, alpha=0.88)
            for peak in str(row["dropped_doublet_peak_mV"]).split(";"):
                if peak:
                    ax.axvline(float(peak), color="0.45", linestyle=":", linewidth=0.9, alpha=0.78)
            for peak in str(row["dropped_leading_peak_mV"]).split(";"):
                if peak:
                    ax.axvline(float(peak), color="0.45", linestyle="--", linewidth=0.9, alpha=0.78)
            ax.set_xlim(70, 450)
            ax.set_title(f"{row['bias_V']:.2f} V, gap={row['height_gap_mV']:.1f} mV", fontsize=12)
            ax.grid(True, alpha=0.22)
            ax.tick_params(labelsize=9)
        fig.suptitle(f"{sipm}: pulse-height spectra", fontsize=17)
        fig.supxlabel("Prompt peak height (mV)", fontsize=12)
        fig.supylabel("Events/bin", fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        out = OUT_DIR / f"{sipm.lower()}_rediscovered_height_peak_histograms.png"
        fig.savefig(out)
        plt.close(fig)
        paths.append(out)
    return paths


def plot_vbr(summaries: pd.DataFrame, fit_results: pd.DataFrame) -> Path:
    plt.style.use("default")
    fig, axes = plt.subplots(1, 2, figsize=(16.2, 6.6), dpi=170, sharey=True)
    for ax, sipm in zip(axes, ["SIPM1", "SIPM2"]):
        fit_points = summaries[(summaries["sipm"].eq(sipm)) & (summaries["group"].eq("fit"))].sort_values("bias_V")
        checks = summaries[(summaries["sipm"].eq(sipm)) & (summaries["group"].eq("channel-switched check"))].sort_values("bias_V")
        result = fit_results[fit_results["sipm"].eq(sipm)].iloc[0]
        fit_channel = str(fit_points.iloc[0]["scope_channel"])
        fit_style = CHANNEL_STYLE.get(fit_channel, CHANNEL_STYLE["A"])
        ax.errorbar(
            fit_points["bias_V"],
            fit_points["height_gap_mV"],
            yerr=fit_points["height_gap_unc_mV"],
            fmt="o",
            capsize=4,
            color=fit_style["color"],
            markersize=5.5,
            linewidth=1.1,
            label=f"fit points: {fit_style['label']}",
            zorder=4,
        )
        if not checks.empty:
            for channel, channel_points in checks.groupby("scope_channel"):
                check_style = CHANNEL_STYLE.get(str(channel), CHANNEL_STYLE["A"])
                ax.errorbar(
                    channel_points["bias_V"],
                    channel_points["height_gap_mV"],
                    yerr=channel_points["height_gap_unc_mV"],
                    fmt="D",
                    mfc="white",
                    mec=check_style["color"],
                    ecolor=check_style["color"],
                    color=check_style["color"],
                    capsize=4,
                    markersize=6.0,
                    markeredgewidth=1.35,
                    linewidth=1.05,
                    label=f"channel-switched checks: {check_style['label']}",
                    zorder=5,
                )
        x_min = min(float(result["breakdown_voltage_V"]), float(fit_points["bias_V"].min())) - 0.22
        x_max = max(float(fit_points["bias_V"].max()), float(checks["bias_V"].max()) if not checks.empty else float(fit_points["bias_V"].max())) + 0.28
        x_line = np.linspace(x_min, x_max, 300)
        y_line = float(result["slope_mV_per_V"]) * x_line + float(result["intercept_mV"])
        ax.plot(
            x_line,
            y_line,
            color=fit_style["color"],
            linewidth=1.9,
            label=f"linear fit: Vbr={result['breakdown_voltage_V']:.2f} +/- {result['breakdown_voltage_unc_V']:.2f} V, R2={result['r2']:.3f}",
            zorder=2,
        )
        ax.axvline(float(result["breakdown_voltage_V"]), color=fit_style["color"], linestyle="--", linewidth=1.35)
        ax.axhline(0, color="0.30", linewidth=1.0)
        ax.set_title(sipm, fontsize=15)
        ax.set_xlabel("Bias voltage (V)", fontsize=14)
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(-5, 85)
        ax.grid(True, alpha=0.28)
        ax.tick_params(axis="both", labelsize=12)
        ax.legend(fontsize=11, loc="upper left", framealpha=0.9)
    axes[0].set_ylabel("1 p.e. peak-height spacing (mV)", fontsize=14)
    fig.suptitle("Breakdown voltage from p.e. peak-height spacing", fontsize=21)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = OUT_DIR / "rediscovered_vbr_peak_height_line_fit.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict[str, object]] = []
    candidate_tables: list[pd.DataFrame] = []
    hist_data: dict[tuple[str, float, str], tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for run in FIT_RUNS + CHECK_RUNS:
        summary, candidates, hist = analyze_run(run)
        summary_rows.append(summary)
        candidate_tables.append(candidates)
        hist_data[(run.sipm, run.bias_v, run.group)] = hist

    summaries = pd.DataFrame(summary_rows)
    candidates = pd.concat(candidate_tables, ignore_index=True)

    fit_results = []
    for sipm in ["SIPM1", "SIPM2"]:
        fit_points = summaries[(summaries["sipm"].eq(sipm)) & (summaries["group"].eq("fit"))].sort_values("bias_V")
        result = fit_vbr(fit_points)
        result["sipm"] = sipm
        fit_results.append(result)
    fit_results_df = pd.DataFrame(fit_results)

    summaries.to_csv(OUT_DIR / "rediscovered_height_peak_summary.csv", index=False)
    candidates.to_csv(OUT_DIR / "rediscovered_height_peak_candidates.csv", index=False)
    fit_results_df.to_csv(OUT_DIR / "rediscovered_vbr_fit_results.csv", index=False)
    hist_paths = plot_histograms(summaries, hist_data)
    vbr_path = plot_vbr(summaries, fit_results_df)

    print(vbr_path)
    for path in hist_paths:
        print(path)
    print(OUT_DIR / "rediscovered_height_peak_summary.csv")
    print(OUT_DIR / "rediscovered_vbr_fit_results.csv")
    print(fit_results_df.to_string(index=False))
    print(summaries[["sipm", "bias_V", "group", "scope_channel", "height_gap_mV", "height_gap_unc_mV", "selected_peak_count", "selected_height_peak_mV", "dropped_doublet_peak_mV", "dropped_leading_peak_mV"]].to_string(index=False))


if __name__ == "__main__":
    main()
