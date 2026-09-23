#!/usr/bin/env python3
from __future__ import annotations

import math
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


ROOT = Path("/Users/tharinduhettiarachchi/Desktop/Design/new415")
OUT_DIR = ROOT / "height_gap_linear_calibration_20260902"

SIPM1_ROOT = ROOT / "robust_sipm1_pe_analysis_20260901"
SIPM2_ROOT = ROOT / "robust_sipm2_area_analysis_20260901"

CHANNEL_STYLE = {
    "A": {"color": "C0", "label": "detector CH0 / scope A"},
    "B": {"color": "C1", "label": "detector CH2 / scope B"},
}


def fit_line(points: pd.DataFrame, y_col: str) -> dict[str, float]:
    x = points["bias_V"].to_numpy(float)
    y = points[y_col].to_numpy(float)
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
        "n_fit_points": int(len(points)),
        "slope_mV_per_V": float(slope),
        "intercept_mV": float(intercept),
        "breakdown_voltage_V": float(vbr),
        "breakdown_voltage_unc_V": vbr_unc,
        "r2": float(1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan,
        "fit_rms_mV": float(np.sqrt(np.mean((y - pred) ** 2))),
    }


def read_corrected_points() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sipm1_fit = pd.read_csv(SIPM1_ROOT / "sipm1_breakdown_points_excluding_53p20_54p10.csv")
    sipm1_all = pd.read_csv(SIPM1_ROOT / "sipm1_height_peak_summary_high_bias_distance.csv")
    sipm2_fit = pd.read_csv(SIPM2_ROOT / "sipm2_breakdown_points_excluding_53p20_54p10.csv")
    sipm2_all = pd.read_csv(SIPM2_ROOT / "sipm2_height_peak_summary_high_bias_distance.csv")

    fit_rows = []
    shown_rows = []
    for sipm, fit_df, all_df in [("SIPM1", sipm1_fit, sipm1_all), ("SIPM2", sipm2_fit, sipm2_all)]:
        fit_biases = set(np.round(fit_df["bias_V"].to_numpy(float), 4))
        for _, row in fit_df.iterrows():
            fit_rows.append(
                {
                    "sipm": sipm,
                    "bias_V": float(row["bias_V"]),
                    "height_gap_mV": float(row["one_pe_spacing_mV"]),
                    "height_gap_unc_mV": float(row["one_pe_spacing_unc_mV"]),
                    "selected_height_peak_mV": row["selected_peak_mus_mV"],
                }
            )
        for _, row in all_df.iterrows():
            bias = float(row["bias_V"])
            if round(bias, 4) in fit_biases:
                continue
            if not np.isfinite(row["one_pe_spacing_mV"]):
                continue
            shown_rows.append(
                {
                    "sipm": sipm,
                    "bias_V": bias,
                    "height_gap_mV": float(row["one_pe_spacing_mV"]),
                    "height_gap_unc_mV": float(row["one_pe_spacing_unc_mV"]),
                    "selected_height_peak_mV": row["selected_peak_mus_mV"],
                    "reason": "diagnostic",
                }
            )

    swapped_all = pd.read_csv(OUT_DIR / "height_gap_points.csv")
    swapped = swapped_all[swapped_all["channel_state"].eq("after channel swap")].copy()
    return pd.DataFrame(fit_rows), pd.DataFrame(shown_rows), swapped


def display_yerr(values: pd.Series, gaps: pd.Series) -> pd.Series:
    out = pd.to_numeric(values, errors="coerce").copy()
    gap_values = pd.to_numeric(gaps, errors="coerce")
    out[(~np.isfinite(out)) | (out > 0.5 * gap_values.abs()) | (out > 30.0)] = np.nan
    return out


def make_plot(
    fit_points: pd.DataFrame,
    shown_points: pd.DataFrame,
    swapped: pd.DataFrame,
    fit_results_df: pd.DataFrame,
    out: Path,
    *,
    show_diagnostic: bool = True,
    emphasize_swapped: bool = False,
) -> Path:
    colors = {"SIPM1": "C0", "SIPM2": "C1"}
    plt.style.use("default")
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6), dpi=170, sharey=True)
    for ax, sipm in zip(axes, ("SIPM1", "SIPM2")):
        used = fit_points[fit_points["sipm"].eq(sipm)].sort_values("bias_V")
        shown = shown_points[shown_points["sipm"].eq(sipm)].sort_values("bias_V")
        swap = swapped[swapped["sipm"].eq(sipm)].sort_values("bias_V")
        ax.errorbar(
            used["bias_V"],
            used["height_gap_mV"],
            yerr=display_yerr(used["height_gap_unc_mV"], used["height_gap_mV"]),
            fmt="o",
            capsize=4,
            color=colors[sipm],
            alpha=0.55 if emphasize_swapped else 1.0,
            markersize=5.2,
            zorder=3,
            label="used in fit",
        )
        if show_diagnostic and not shown.empty:
            ax.errorbar(
                shown["bias_V"],
                shown["height_gap_mV"],
                yerr=display_yerr(shown["height_gap_unc_mV"], shown["height_gap_mV"]),
                fmt="s",
                mfc="none",
                capsize=4,
                color="0.45",
                alpha=0.55 if emphasize_swapped else 1.0,
                zorder=2,
                label="shown, not fit",
            )
        if not swap.empty:
            ax.errorbar(
                swap["bias_V"],
                swap["height_gap_mV"],
                yerr=display_yerr(swap["height_gap_unc_mV"], swap["height_gap_mV"]),
                fmt="*",
                mfc="gold" if emphasize_swapped else "white",
                mec=colors[sipm],
                ecolor=colors[sipm],
                capsize=4,
                color=colors[sipm],
                markersize=13 if emphasize_swapped else 7,
                markeredgewidth=1.4 if emphasize_swapped else 1.0,
                elinewidth=1.35 if emphasize_swapped else 1.0,
                zorder=9,
                label="Sep 2 check",
            )
            if emphasize_swapped:
                for _, row in swap.iterrows():
                    ax.annotate(
                        f"{row['bias_V']:.2f} V",
                        (float(row["bias_V"]), float(row["height_gap_mV"])),
                        xytext=(7, -12),
                        textcoords="offset points",
                        fontsize=8,
                        color=colors[sipm],
                        bbox=dict(facecolor="white", edgecolor=colors[sipm], alpha=0.78, linewidth=0.7),
                        zorder=10,
                    )
        res = fit_results_df[fit_results_df["sipm"].eq(sipm)].iloc[0]
        xmax_candidates = [float(used["bias_V"].max())]
        if not swap.empty:
            xmax_candidates.append(float(swap["bias_V"].max()))
        x_min = min(float(used["bias_V"].min()), float(res["breakdown_voltage_V"])) - 0.12
        x_max = max(xmax_candidates) + 0.18
        x_line = np.linspace(x_min, x_max, 220)
        y_line = float(res["slope_mV_per_V"]) * x_line + float(res["intercept_mV"])
        ax.plot(
            x_line,
            y_line,
            color=colors[sipm],
            linewidth=1.6,
            alpha=0.68 if emphasize_swapped else 1.0,
            zorder=1,
            label=f"fit: Vbr={res['breakdown_voltage_V']:.2f} +/- {res['breakdown_voltage_unc_V']:.2f} V, R2={res['r2']:.3f}",
        )
        ax.axvline(float(res["breakdown_voltage_V"]), color=colors[sipm], linestyle="--", linewidth=1.1)
        ax.axhline(0, color="0.35", linewidth=1.0)
        ax.set_title(sipm)
        ax.set_xlabel("Bias voltage (V)")
        ax.grid(True, alpha=0.28)
        ax.legend(fontsize=8, loc="upper left")
    axes[0].set_ylabel("1 p.e. peak-height spacing (mV)")
    fig.suptitle("Breakdown voltage from peak-height p.e. spacing", fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out)
    plt.close(fig)
    return out


def make_publication_linear_plot(rows: pd.DataFrame, fit_results_df: pd.DataFrame, out: Path) -> Path:
    plt.style.use("default")
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 6.0), dpi=180, sharey=True)
    for ax, sipm in zip(axes, ("SIPM1", "SIPM2")):
        sub = rows[rows["sipm"].eq(sipm)].sort_values(["bias_V", "group"])
        used = sub[sub["group"].eq("fit")]
        switched = sub[sub["group"].eq("Sep 2 check")]
        res = fit_results_df[fit_results_df["sipm"].eq(sipm)].iloc[0]
        used_channel = str(used.iloc[0].get("scope_channel", active_channel(Path(str(used.iloc[0]["source"])), sipm))) if not used.empty else "A"
        switched_channel = str(switched.iloc[0].get("scope_channel", active_channel(Path(str(switched.iloc[0]["source"])), sipm))) if not switched.empty else "B"
        used_style = CHANNEL_STYLE.get(used_channel, CHANNEL_STYLE["A"])
        switched_style = CHANNEL_STYLE.get(switched_channel, CHANNEL_STYLE["B"])

        ax.errorbar(
            used["bias_V"],
            used["height_gap_mV"],
            yerr=display_yerr(used["height_gap_unc_mV"], used["height_gap_mV"]),
            fmt="o",
            markersize=6.8,
            capsize=4.5,
            linewidth=1.25,
            color=used_style["color"],
            label=f"fit points: {used_style['label']}",
            zorder=4,
        )
        if not switched.empty:
            ax.errorbar(
                switched["bias_V"],
                switched["height_gap_mV"],
                yerr=display_yerr(switched["height_gap_unc_mV"], switched["height_gap_mV"]),
                fmt="D",
                markersize=6.5,
                markerfacecolor="white",
                markeredgewidth=1.5,
                capsize=4.5,
                linewidth=1.25,
                color=switched_style["color"],
                ecolor=switched_style["color"],
                label=f"channel-switched checks: {switched_style['label']}",
                zorder=5,
            )

        x_min = min(float(used["bias_V"].min()), float(res["breakdown_voltage_V"])) - 0.12
        x_max = max(float(used["bias_V"].max()), float(switched["bias_V"].max()) if not switched.empty else float(used["bias_V"].max())) + 0.18
        x_line = np.linspace(x_min, x_max, 260)
        y_line = float(res["slope_mV_per_V"]) * x_line + float(res["intercept_mV"])
        ax.plot(
            x_line,
            y_line,
            color=used_style["color"],
            linewidth=2.0,
            label=(
                f"linear fit: Vbr={res['breakdown_voltage_V']:.2f} +/- {res['breakdown_voltage_unc_V']:.2f} V, "
                f"slope={res['slope_mV_per_V']:.2f} mV/V, R2={res['r2']:.3f}"
            ),
            zorder=3,
        )
        ax.axvline(float(res["breakdown_voltage_V"]), color=used_style["color"], linestyle="--", linewidth=1.4)
        ax.axhline(0, color="0.35", linewidth=1.1)
        ax.set_title(sipm, fontsize=15)
        ax.set_xlabel("Bias voltage (V)", fontsize=13)
        ax.tick_params(axis="both", labelsize=11)
        ax.grid(True, alpha=0.28)
        ax.legend(fontsize=9.5, loc="upper left", framealpha=0.88)
    axes[0].set_ylabel("1 p.e. peak-height spacing (mV)", fontsize=13)
    fig.suptitle("Breakdown voltage from p.e. peak-height spacing", fontsize=18)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out)
    plt.close(fig)
    return out


def parse_peak_string(text: object) -> list[float]:
    if pd.isna(text):
        return []
    peaks: list[float] = []
    for part in str(text).replace(",", ";").split(";"):
        part = part.strip()
        if not part:
            continue
        try:
            peaks.append(float(part))
        except ValueError:
            continue
    return peaks


def active_channel(source: Path, sipm: str) -> str:
    try:
        df = pd.read_csv(source, nrows=8)
        if "active_scope_channel" in df.columns:
            values = df["active_scope_channel"].dropna().astype(str).str.strip()
            values = values[values.isin(["A", "B"])]
            if not values.empty:
                return str(values.iloc[0])
    except Exception:
        pass
    return "A" if sipm == "SIPM1" else "B"


def height_values(source: Path, sipm: str) -> np.ndarray:
    df = pd.read_csv(source)
    accepted = df["accepted"].fillna(False).astype(bool) if "accepted" in df.columns else pd.Series(True, index=df.index)
    channel = active_channel(source, sipm)
    col = "prompt_peak_mV" if "prompt_peak_mV" in df.columns else f"{channel}_peak_mV"
    values = pd.to_numeric(df.loc[accepted, col], errors="coerce").dropna().to_numpy(float)
    return values[np.isfinite(values)]


def fd_edges(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 20:
        return np.linspace(0.0, 1.0, 20)
    lo, hi = np.percentile(x, [0.15, 99.35])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(np.nanmin(x)), float(np.nanmax(x))
    pad = max(5.0, 0.045 * (hi - lo))
    lo -= pad
    hi += pad
    core = x[(x >= lo) & (x <= hi)]
    q25, q75 = np.percentile(core, [25, 75])
    width = 2.0 * (q75 - q25) / np.cbrt(max(1, core.size))
    if not np.isfinite(width) or width <= 0:
        width = (hi - lo) / 120.0
    width = float(np.clip(width, (hi - lo) / 260.0, (hi - lo) / 45.0))
    start = math.floor(lo / width) * width
    stop = math.ceil(hi / width) * width
    return np.arange(start, stop + width, width)


def gaussian_linear_bg(x: np.ndarray, amp: float, mu: float, sigma: float, bg0: float, bg1: float) -> np.ndarray:
    return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2) + bg0 + bg1 * (x - mu)


def refine_histogram_peaks(
    centers: np.ndarray,
    counts: np.ndarray,
    peak_indexes: np.ndarray,
    width: float,
    smooth: np.ndarray | None = None,
) -> pd.DataFrame:
    rows: list[dict[str, float]] = []
    for idx in peak_indexes:
        mu0 = float(centers[idx])
        observed_height = float(smooth[idx]) if smooth is not None else float(counts[idx])
        half_width = max(9.0, 5.0 * width)
        mask = (centers >= mu0 - half_width) & (centers <= mu0 + half_width)
        if int(mask.sum()) < 7:
            rows.append(
                {
                    "mu_mV": mu0,
                    "mu_err_mV": width / np.sqrt(12.0),
                    "sigma_mV": np.nan,
                    "amp": float(counts[idx]),
                    "observed_height": observed_height,
                    "raw_count": float(counts[idx]),
                }
            )
            continue
        xx = centers[mask]
        yy = counts[mask].astype(float)
        amp0 = max(1.0, float(counts[idx] - np.percentile(yy, 15)))
        sigma0 = max(1.5, 1.5 * width)
        bg0 = max(0.0, float(np.percentile(yy, 15)))
        try:
            popt, pcov = curve_fit(
                gaussian_linear_bg,
                xx,
                yy,
                p0=[amp0, mu0, sigma0, bg0, 0.0],
                bounds=([0.0, mu0 - half_width, 0.25, 0.0, -np.inf], [np.inf, mu0 + half_width, 40.0, np.inf, np.inf]),
                maxfev=20000,
            )
            perr = np.sqrt(np.diag(pcov))
            mu = float(popt[1])
            mu_err = float(perr[1]) if np.isfinite(perr[1]) and perr[1] < 20.0 else width / np.sqrt(12.0)
            sigma = abs(float(popt[2]))
            amp = float(popt[0])
        except Exception:
            mu = mu0
            mu_err = width / np.sqrt(12.0)
            sigma = np.nan
            amp = float(counts[idx])
        rows.append(
            {
                "mu_mV": mu,
                "mu_err_mV": mu_err,
                "sigma_mV": sigma,
                "amp": amp,
                "observed_height": observed_height,
                "raw_count": float(counts[idx]),
            }
        )
    return pd.DataFrame(rows).sort_values("mu_mV").reset_index(drop=True)


def remove_first_peak_if_second_is_larger(fits: pd.DataFrame) -> tuple[pd.DataFrame, float | None]:
    if len(fits) < 2:
        return fits, None
    height_col = "observed_height" if "observed_height" in fits.columns else "amp"
    heights = pd.to_numeric(fits[height_col], errors="coerce").to_numpy(float)
    if np.isfinite(heights[0]) and np.isfinite(heights[1]) and heights[1] > heights[0]:
        ignored_mu = float(fits.iloc[0]["mu_mV"])
        return fits.iloc[1:].reset_index(drop=True), ignored_mu
    return fits, None


def select_expected_spacing(fits: pd.DataFrame, expected_gap: float) -> tuple[list[int], float, float, float]:
    if len(fits) < 2 or not np.isfinite(expected_gap):
        return [], np.nan, np.nan, np.nan
    mus = fits["mu_mV"].to_numpy(float)
    best: tuple[float, tuple[int, ...], float, float] | None = None
    for count in range(min(5, len(mus)), 1, -1):
        for combo in combinations(range(len(mus)), count):
            vals = mus[list(combo)]
            diffs = np.diff(vals)
            if diffs.size == 0:
                continue
            if np.any(diffs < 0.70 * expected_gap) or np.any(diffs > 1.35 * expected_gap):
                continue
            spacing = float(np.mean(diffs))
            residual = float(np.sqrt(np.mean((diffs - expected_gap) ** 2)))
            scatter = float(np.std(diffs, ddof=1) / np.sqrt(diffs.size)) if diffs.size > 1 else abs(spacing - expected_gap)
            coverage = float(vals[-1] - vals[0])
            score = residual / expected_gap - 0.030 * count - 0.0004 * coverage
            if best is None or score < best[0]:
                best = (score, combo, spacing, scatter)
        if best is not None and count >= 3:
            break
    if best is None:
        return [], np.nan, np.nan, np.nan
    _, combo, spacing, scatter_unc = best
    selected = list(combo)
    selected_errs = pd.to_numeric(fits.loc[selected, "mu_err_mV"], errors="coerce").to_numpy(float)
    finite_errs = selected_errs[np.isfinite(selected_errs)]
    fit_unc = float(np.sqrt(np.sum(finite_errs**2)) / max(1, len(selected) - 1)) if finite_errs.size >= 2 else np.nan
    if np.isfinite(scatter_unc) and np.isfinite(fit_unc):
        unc = float(np.hypot(scatter_unc, fit_unc))
    elif np.isfinite(scatter_unc):
        unc = scatter_unc
    elif np.isfinite(fit_unc):
        unc = fit_unc
    else:
        unc = np.nan
    residual = float(np.sqrt(np.mean((np.diff(mus[selected]) - spacing) ** 2))) if len(selected) > 2 else 0.0
    return selected, spacing, unc, residual


def recover_gap_from_histogram(source: Path, sipm: str, expected_gap: float) -> dict[str, object] | None:
    values = height_values(source, sipm)
    if values.size < 100:
        return None
    edges = fd_edges(values)
    counts, edges = np.histogram(values, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    width = float(np.median(np.diff(edges))) if len(edges) > 1 else np.nan
    smooth = gaussian_filter1d(counts.astype(float), sigma=1.15)
    prominence = max(6.0, 0.004 * float(np.nanmax(smooth)))
    distance = max(3, int(round(0.46 * expected_gap / width))) if np.isfinite(width) and width > 0 else 3
    peak_indexes, props = find_peaks(smooth, prominence=prominence, distance=distance)
    high_edge = min(float(edges[-1]), 520.0)
    good = [idx for idx in peak_indexes if 55.0 <= centers[idx] <= high_edge and smooth[idx] > 0]
    if len(good) < 2:
        return None
    fits = refine_histogram_peaks(centers, counts, np.asarray(good, dtype=int), width, smooth)
    fits, ignored_first_peak = remove_first_peak_if_second_is_larger(fits)
    selected, spacing, spacing_unc, residual = select_expected_spacing(fits, expected_gap)
    if len(selected) < 2 or not np.isfinite(spacing):
        return None
    selected_mus = fits.loc[selected, "mu_mV"].to_numpy(float)
    diffs = np.diff(selected_mus)
    return {
        "height_gap_mV": spacing,
        "height_gap_unc_mV": spacing_unc,
        "height_gap_diffs_mV": ";".join(f"{x:.2f}" for x in diffs),
        "selected_peak_count": int(len(selected)),
        "selected_height_peak_mV": ";".join(f"{x:.3f}" for x in selected_mus),
        "spacing_residual_rms_mV": residual,
        "recovered_from_histogram": True,
        "first_peak_ignored_second_larger": ignored_first_peak is not None,
        "ignored_first_peak_mV": ignored_first_peak,
    }


def recover_missing_swapped_points(swapped: pd.DataFrame, fit_results_df: pd.DataFrame) -> pd.DataFrame:
    out = swapped.copy()
    out["recovered_from_histogram"] = False
    for idx, row in out.iterrows():
        if pd.notna(row.get("height_gap_mV")) and np.isfinite(float(row["height_gap_mV"])):
            continue
        sipm = str(row["sipm"])
        source = Path(str(row["source"]))
        if not source.exists():
            continue
        res = fit_results_df[fit_results_df["sipm"].eq(sipm)]
        if res.empty:
            continue
        expected_gap = float(res.iloc[0]["slope_mV_per_V"]) * float(row["bias_V"]) + float(res.iloc[0]["intercept_mV"])
        recovered = recover_gap_from_histogram(source, sipm, expected_gap)
        if recovered is None:
            continue
        for key, value in recovered.items():
            out.at[idx, key] = value
    return out


def expected_gap_for_row(row: pd.Series, fit_results_df: pd.DataFrame) -> float:
    res = fit_results_df[fit_results_df["sipm"].eq(str(row["sipm"]))]
    if res.empty:
        return np.nan
    return float(res.iloc[0]["slope_mV_per_V"]) * float(row["bias_V"]) + float(res.iloc[0]["intercept_mV"])


def improve_histogram_rows(rows: pd.DataFrame, fit_results_df: pd.DataFrame) -> pd.DataFrame:
    improved = rows.copy()
    improved["peak_selection_method"] = "summary"
    for idx, row in improved.iterrows():
        source = Path(str(row["source"]))
        sipm = str(row["sipm"])
        expected_gap = expected_gap_for_row(row, fit_results_df)
        recovered = recover_gap_from_histogram(source, sipm, expected_gap)
        if recovered is None:
            continue

        old_gap = pd.to_numeric(pd.Series([row["height_gap_mV"]]), errors="coerce").iloc[0]
        old_peak_count = len(parse_peak_string(row["selected_height_peak_mV"]))
        rec_gap = float(recovered["height_gap_mV"])
        rec_peak_count = int(recovered["selected_peak_count"])
        old_error = abs(float(old_gap) - expected_gap) if np.isfinite(old_gap) else np.inf
        rec_error = abs(rec_gap - expected_gap)
        margin = max(2.0, 0.05 * abs(expected_gap))

        accept_recovery = (
            not np.isfinite(old_gap)
            or bool(recovered.get("first_peak_ignored_second_larger", False))
            or rec_error < old_error - 0.25 * margin
            or (rec_peak_count > old_peak_count and rec_error <= old_error + margin)
            or (rec_peak_count == old_peak_count and rec_error <= old_error + 0.5 * margin)
        )
        if not accept_recovery:
            continue
        for key, value in recovered.items():
            improved.at[idx, key] = value
        improved.at[idx, "peak_selection_method"] = "expected-spacing recovery"
    return improved


def source_rows_for_histograms(fit_points: pd.DataFrame, swapped: pd.DataFrame) -> pd.DataFrame:
    manifest_points = pd.read_csv(OUT_DIR / "height_gap_points.csv")
    rows: list[dict[str, object]] = []

    for _, row in fit_points.iterrows():
        sipm = str(row["sipm"])
        bias = float(row["bias_V"])
        candidates = manifest_points[
            manifest_points["sipm"].eq(sipm)
            & manifest_points["channel_state"].eq("original channels")
            & np.isclose(pd.to_numeric(manifest_points["bias_V"], errors="coerce"), bias, atol=0.002)
            & manifest_points["height_gap_mV"].notna()
        ].copy()
        if candidates.empty:
            continue
        candidates = candidates.sort_values(["use_in_linear_fit", "n_rows"], ascending=[False, False])
        source = Path(str(candidates.iloc[0]["source"]))
        if not source.exists():
            continue
        scope_channel = active_channel(source, sipm)
        rows.append(
            {
                "sipm": sipm,
                "bias_V": bias,
                "source": str(source),
                "scope_channel": scope_channel,
                "height_gap_mV": float(row["height_gap_mV"]),
                "height_gap_unc_mV": float(row["height_gap_unc_mV"]),
                "selected_height_peak_mV": row["selected_height_peak_mV"],
                "group": "fit",
            }
        )

    for _, row in swapped.iterrows():
        if not np.isfinite(float(row["height_gap_mV"])):
            continue
        source = Path(str(row["source"]))
        if not source.exists():
            continue
        scope_channel = active_channel(source, str(row["sipm"]))
        rows.append(
            {
                "sipm": str(row["sipm"]),
                "bias_V": float(row["bias_V"]),
                "source": str(source),
                "scope_channel": scope_channel,
                "height_gap_mV": float(row["height_gap_mV"]),
                "height_gap_unc_mV": float(row["height_gap_unc_mV"]),
                "selected_height_peak_mV": row["selected_height_peak_mV"],
                "group": "Sep 2 check",
            }
        )

    out = pd.DataFrame(rows).sort_values(["sipm", "bias_V", "group", "source"]).reset_index(drop=True)
    return out


def common_edges_for_sipm(sub: pd.DataFrame, values_by_source: dict[str, np.ndarray]) -> np.ndarray:
    lows: list[float] = []
    highs: list[float] = []
    widths: list[float] = []
    for _, row in sub.iterrows():
        source = str(row["source"])
        values = values_by_source.get(source)
        if values is None or values.size < 20:
            continue
        edges = fd_edges(values)
        if len(edges) < 2:
            continue
        lows.append(float(edges[0]))
        highs.append(float(edges[-1]))
        widths.append(float(np.median(np.diff(edges))))
    if not lows or not highs or not widths:
        return np.linspace(0.0, 1.0, 20)
    width = float(np.nanmin(widths))
    lo = math.floor(min(lows) / width) * width
    hi = math.ceil(max(highs) / width) * width
    return np.arange(lo, hi + width, width)


def write_kept_histograms(rows: pd.DataFrame) -> list[Path]:
    paths: list[Path] = []
    for sipm in ("SIPM1", "SIPM2"):
        sub = rows[rows["sipm"].eq(sipm)].sort_values(["bias_V", "group", "source"]).reset_index(drop=True)
        if sub.empty:
            continue
        values_by_source = {
            str(row["source"]): height_values(Path(str(row["source"])), str(row["sipm"]))
            for _, row in sub.iterrows()
        }
        common_edges = common_edges_for_sipm(sub, values_by_source)
        ncols = 3
        nrows = int(math.ceil(len(sub) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(16, max(3.55 * nrows, 4.6)), dpi=170)
        axes_flat = np.asarray(axes).reshape(-1)
        for ax, (_, row) in zip(axes_flat, sub.iterrows()):
            source = Path(str(row["source"]))
            values = values_by_source.get(str(source), np.array([]))
            if values.size < 50:
                ax.axis("off")
                continue
            edges = common_edges
            counts, edges = np.histogram(values, bins=edges)
            centers = 0.5 * (edges[:-1] + edges[1:])
            smooth = gaussian_filter1d(counts.astype(float), sigma=1.15)
            is_switched = row["group"] == "Sep 2 check"
            channel = str(row.get("scope_channel", active_channel(source, str(row["sipm"]))))
            style = CHANNEL_STYLE.get(channel, CHANNEL_STYLE["A"])
            color = style["color"]
            data_label = style["label"]
            if is_switched:
                data_label += ", channel switched"

            ax.step(centers, counts, where="mid", color=color, linewidth=1.05, alpha=0.78, label=data_label)
            ax.plot(centers, smooth, color="black", linewidth=1.20, alpha=0.88, label="smoothed histogram")
            for mu in parse_peak_string(row["selected_height_peak_mV"]):
                ax.axvline(float(mu), color="C3", linewidth=1.25, alpha=0.88)
            group = str(row["group"])
            gap = float(row["height_gap_mV"])
            unc = float(row["height_gap_unc_mV"]) if pd.notna(row["height_gap_unc_mV"]) else np.nan
            if np.isfinite(unc) and unc <= 0.5 * abs(gap):
                legend_title = f"{gap:.1f} +/- {unc:.1f} mV"
            else:
                legend_title = f"{gap:.1f} mV"
            ax.legend(title=legend_title, fontsize=8.8, title_fontsize=9.0, loc="upper right", framealpha=0.88)
            label = f"{row['bias_V']:.2f} V"
            if group == "Sep 2 check":
                label += ", channel switched"
            ax.set_title(label, fontsize=11.5)
            ax.set_xlabel("Prompt peak height (mV)", fontsize=10.5)
            ax.set_ylabel("Events/bin", fontsize=10.5)
            ax.set_xlim(float(edges[0]), 450.0)
            ax.set_xticks([100, 200, 300, 400, 450])
            ax.tick_params(axis="both", labelsize=9.5)
            ax.grid(True, alpha=0.24)
        for ax in axes_flat[len(sub) :]:
            ax.axis("off")
        fig.suptitle(f"{sipm}: p.e. peak-height spectra", fontsize=18)
        fig.tight_layout(rect=[0, 0, 1, 0.958])
        out = OUT_DIR / f"{sipm.lower()}_height_gap_histogram_panel_improved_publication.png"
        fig.savefig(out)
        plt.close(fig)
        paths.append(out)
    return paths


def write_recovered_zoom_plots(swapped: pd.DataFrame) -> list[Path]:
    paths: list[Path] = []
    recovered = swapped[swapped.get("recovered_from_histogram", False).fillna(False).astype(bool)].copy()
    for _, row in recovered.iterrows():
        sipm = str(row["sipm"])
        source = Path(str(row["source"]))
        values = height_values(source, sipm)
        if values.size < 50:
            continue
        edges = fd_edges(values)
        counts, edges = np.histogram(values, bins=edges)
        centers = 0.5 * (edges[:-1] + edges[1:])
        smooth = gaussian_filter1d(counts.astype(float), sigma=1.15)
        peaks = parse_peak_string(row["selected_height_peak_mV"])

        plt.style.use("default")
        fig, axes = plt.subplots(2, 1, figsize=(10.5, 7.2), dpi=170, sharex=True)
        for ax in axes:
            ax.step(centers, counts, where="mid", color="C4", linewidth=0.95, alpha=0.72, label="height histogram")
            ax.plot(centers, smooth, color="black", linewidth=1.15, alpha=0.9, label="smoothed")
            for peak in peaks:
                ax.axvline(float(peak), color="C3", linewidth=1.25, alpha=0.9)
            ax.set_xlim(float(edges[0]), float(edges[-1]))
            ax.set_ylabel("Events/bin")
            ax.grid(True, alpha=0.24)
        if len(peaks) > 1:
            ymax_zoom = max(150.0, min(float(np.nanmax(counts)) * 0.23, 1.45 * float(np.nanmax(counts[centers > peaks[0] + 0.35 * row["height_gap_mV"]]))))
            axes[1].set_ylim(0, ymax_zoom)
        axes[1].set_xlabel("Prompt peak height (mV)")
        gap = float(row["height_gap_mV"])
        unc = float(row["height_gap_unc_mV"]) if pd.notna(row["height_gap_unc_mV"]) else np.nan
        label = f"{sipm} 57.87 V recovered p.e. spacing"
        fig.suptitle(label, fontsize=15)
        axes[0].legend(
            title=f"gap = {gap:.1f} +/- {unc:.1f} mV" if np.isfinite(unc) else f"gap = {gap:.1f} mV",
            fontsize=8,
            title_fontsize=8,
            loc="upper right",
        )
        axes[1].legend(fontsize=8, loc="upper right")
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        out = OUT_DIR / f"{sipm.lower()}_57p87_recovered_height_peak_zoom.png"
        fig.savefig(out)
        plt.close(fig)
        paths.append(out)
    return paths


def write_example_zoom_histogram(rows: pd.DataFrame) -> Path | None:
    case = rows[(rows["sipm"].eq("SIPM1")) & (np.isclose(rows["bias_V"], 56.86)) & (rows["group"].eq("fit"))]
    if case.empty:
        return None
    row = case.iloc[0]
    source = Path(str(row["source"]))
    values = height_values(source, "SIPM1")
    if values.size < 50:
        return None

    edges = fd_edges(values)
    counts, edges = np.histogram(values, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    smooth = gaussian_filter1d(counts.astype(float), sigma=1.15)
    peaks = parse_peak_string(row["selected_height_peak_mV"])
    channel = str(row.get("scope_channel", active_channel(source, "SIPM1")))
    color = CHANNEL_STYLE.get(channel, CHANNEL_STYLE["A"])["color"]

    plt.style.use("default")
    fig, ax = plt.subplots(figsize=(11.5, 6.4), dpi=180)
    ax.step(centers, counts, where="mid", color=color, linewidth=1.35, alpha=0.72, label="height histogram")
    ax.plot(centers, smooth, color="black", linewidth=1.8, alpha=0.92, label="smoothed histogram")
    for i, peak in enumerate(peaks):
        ax.axvline(float(peak), color="C3", linewidth=1.65, alpha=0.9, label="selected p.e. peaks" if i == 0 else None)
    gap = float(row["height_gap_mV"])
    unc = float(row["height_gap_unc_mV"]) if pd.notna(row["height_gap_unc_mV"]) else np.nan
    legend_title = f"1 p.e. spacing = {gap:.1f} +/- {unc:.1f} mV" if np.isfinite(unc) else f"1 p.e. spacing = {gap:.1f} mV"
    ax.legend(title=legend_title, fontsize=11, title_fontsize=11.5, loc="upper right", framealpha=0.9)
    ax.set_xlim(160, 450)
    ax.set_ylim(0, 560)
    ax.set_xticks([160, 200, 240, 280, 320, 360, 400, 450])
    ax.set_xlabel("Prompt peak height (mV)", fontsize=13)
    ax.set_ylabel("Events/bin", fontsize=13)
    ax.set_title("SIPM1 at 56.86 V: p.e. peak-height spectrum", fontsize=17)
    ax.tick_params(axis="both", labelsize=11)
    ax.grid(True, alpha=0.26)
    fig.tight_layout()
    out = OUT_DIR / "sipm1_56p86_peak_height_histogram_zoom.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fit_points, shown_points, swapped = read_corrected_points()
    base_fit_results = []
    for sipm in ("SIPM1", "SIPM2"):
        sub = fit_points[fit_points["sipm"].eq(sipm)].sort_values("bias_V")
        result = fit_line(sub, "height_gap_mV")
        result["sipm"] = sipm
        base_fit_results.append(result)
    base_fit_results_df = pd.DataFrame(base_fit_results)
    swapped_no_55x = swapped[~((swapped["bias_V"] >= 55.0) & (swapped["bias_V"] < 56.0))].copy()
    swapped_no_55x = recover_missing_swapped_points(swapped_no_55x, base_fit_results_df)

    fit_points.to_csv(OUT_DIR / "clean_height_gap_fit_points.csv", index=False)
    shown_points.to_csv(OUT_DIR / "clean_height_gap_diagnostic_points.csv", index=False)
    swapped_no_55x.to_csv(OUT_DIR / "clean_height_gap_swapped_points_no_55x.csv", index=False)
    base_fit_results_df.to_csv(OUT_DIR / "clean_height_gap_linear_fit_results.csv", index=False)

    hist_rows = source_rows_for_histograms(fit_points, swapped_no_55x)
    hist_rows.to_csv(OUT_DIR / "clean_height_gap_kept_histogram_runs.csv", index=False)
    improved_hist_rows = improve_histogram_rows(hist_rows, base_fit_results_df)
    improved_hist_rows.to_csv(OUT_DIR / "clean_height_gap_improved_histogram_runs.csv", index=False)

    improved_fit_points = improved_hist_rows[improved_hist_rows["group"].eq("fit")].copy()
    improved_fit_results = []
    for sipm in ("SIPM1", "SIPM2"):
        sub = improved_fit_points[improved_fit_points["sipm"].eq(sipm)].sort_values("bias_V")
        result = fit_line(sub, "height_gap_mV")
        result["sipm"] = sipm
        improved_fit_results.append(result)
    improved_fit_results_df = pd.DataFrame(improved_fit_results)
    improved_fit_points.to_csv(OUT_DIR / "improved_height_gap_fit_points.csv", index=False)
    improved_fit_results_df.to_csv(OUT_DIR / "improved_height_gap_linear_fit_results.csv", index=False)

    linear_out = make_publication_linear_plot(
        improved_hist_rows,
        improved_fit_results_df,
        OUT_DIR / "sipm_height_gap_vs_bias_improved_peak_separation_publication.png",
    )
    hist_paths = write_kept_histograms(improved_hist_rows)
    zoom_path = write_example_zoom_histogram(improved_hist_rows)
    print(linear_out)
    for hist_path in hist_paths:
        print(hist_path)
    if zoom_path is not None:
        print(zoom_path)
    print(improved_fit_results_df.to_string(index=False))


if __name__ == "__main__":
    main()
