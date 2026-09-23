#!/usr/bin/env python3
from __future__ import annotations

import math
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
from scipy.signal import find_peaks


ROOT = Path("/Users/tharinduhettiarachchi/Desktop/Design/new415")
MANIFEST = ROOT / "waveform_all_run_histogram_panels_20260902/all_distinct_waveform_runs_manifest.csv"
OUT_DIR = ROOT / "area_gap_linear_calibration_20260902"


@dataclass(frozen=True)
class Run:
    sipm: str
    bias_v: float
    label: str
    source: Path


def gaussian_linear_bg(x: np.ndarray, amp: float, mu: float, sigma: float, bg0: float, bg1: float) -> np.ndarray:
    return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2) + bg0 + bg1 * (x - mu)


def read_manifest() -> list[Run]:
    df = pd.read_csv(MANIFEST)
    runs: list[Run] = []
    for _, row in df.iterrows():
        source = Path(str(row["source"]))
        if not source.exists():
            continue
        runs.append(
            Run(
                sipm=str(row["sipm"]).upper(),
                bias_v=float(row["bias_V"]),
                label=str(row["panel_label"]).split(" | ")[0],
                source=source,
            )
        )
    return runs


def active_channel(df: pd.DataFrame, sipm: str) -> str:
    if "active_scope_channel" in df.columns:
        values = df["active_scope_channel"].dropna().astype(str).str.strip()
        values = values[values.isin(["A", "B"])]
        if not values.empty:
            return str(values.iloc[0])
    return "A" if sipm == "SIPM1" else "B"


def saturation_count(source: Path) -> int:
    summary = source.with_name("summary.csv")
    if summary.exists():
        try:
            s = pd.read_csv(summary).iloc[0]
            if "near_rail_saturation_candidates" in s.index and pd.notna(s["near_rail_saturation_candidates"]):
                return int(float(s["near_rail_saturation_candidates"]))
        except Exception:
            pass
    try:
        clipped = pd.read_csv(source, usecols=["clipped"])["clipped"]
        return int(clipped.fillna(False).astype(bool).sum())
    except Exception:
        pass
    return 0


def area_values(run: Run) -> tuple[np.ndarray, str, int, int]:
    df = pd.read_csv(run.source)
    channel = active_channel(df, run.sipm)
    accepted = df["accepted"].fillna(False).astype(bool) if "accepted" in df.columns else pd.Series(True, index=df.index)
    col = "prompt_area_mVns" if "prompt_area_mVns" in df.columns else f"{channel}_area_mVns"
    values = pd.to_numeric(df.loc[accepted, col], errors="coerce").dropna().to_numpy(float)
    values = values[np.isfinite(values)]
    return values, channel, int(len(df)), int(accepted.sum())


def fd_edges(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 20:
        return np.linspace(0.0, 1.0, 20)
    lo, hi = np.percentile(x, [0.15, 99.35])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(np.nanmin(x)), float(np.nanmax(x))
    pad = max(300.0, 0.045 * (hi - lo))
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


def raw_peak_candidates(centers: np.ndarray, counts: np.ndarray, bias_v: float) -> np.ndarray:
    if centers.size < 10 or np.nanmax(counts) <= 0:
        return np.array([], dtype=int)
    bin_width = float(np.median(np.diff(centers)))
    smooth = gaussian_filter1d(counts.astype(float), sigma=1.35)
    prominence = max(8.0, 0.028 * float(np.max(smooth)))
    # A lower bound that grows with over-voltage keeps the high-bias doublet
    # structure from being mistaken for two different p.e. peaks.
    min_spacing_mVns = 800.0 if bias_v < 55.0 else 1150.0
    if bias_v >= 56.0:
        min_spacing_mVns = 1550.0
    distance = max(3, int(round(min_spacing_mVns / max(bin_width, 1e-9))))
    candidates, _ = find_peaks(smooth, prominence=prominence, distance=distance)
    return candidates.astype(int)


def fit_area_peak_centers(values: np.ndarray, bias_v: float) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]:
    edges = fd_edges(values)
    counts, edges = np.histogram(values, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    smooth = gaussian_filter1d(counts.astype(float), sigma=1.35)
    peak_indexes = raw_peak_candidates(centers, counts, bias_v)
    rows: list[dict[str, float]] = []
    bin_width = float(np.median(np.diff(centers))) if centers.size > 1 else 1.0
    for peak_index in peak_indexes:
        mu0 = float(centers[peak_index])
        half_width = max(450.0, 4.2 * bin_width)
        mask = (centers >= mu0 - half_width) & (centers <= mu0 + half_width)
        if int(mask.sum()) < 7:
            continue
        xx = centers[mask]
        yy = counts[mask].astype(float)
        amp0 = max(1.0, float(counts[peak_index] - np.percentile(yy, 15)))
        sigma0 = max(120.0, 1.5 * bin_width)
        bg0 = max(0.0, float(np.percentile(yy, 15)))
        try:
            popt, pcov = curve_fit(
                gaussian_linear_bg,
                xx,
                yy,
                p0=[amp0, mu0, sigma0, bg0, 0.0],
                bounds=(
                    [0.0, mu0 - half_width, 35.0, 0.0, -np.inf],
                    [np.inf, mu0 + half_width, 5000.0, np.inf, np.inf],
                ),
                maxfev=20000,
            )
            perr = np.sqrt(np.diag(pcov))
            mu = float(popt[1])
            mu_err = float(perr[1]) if np.isfinite(perr[1]) else np.nan
            sigma = abs(float(popt[2]))
            amp = float(popt[0])
            if not np.isfinite(mu_err) or mu_err > 0.8 * half_width:
                mu_err = np.nan
        except Exception:
            mu = mu0
            mu_err = np.nan
            sigma = np.nan
            amp = float(counts[peak_index])
        rows.append({"mu_mVns": mu, "mu_err_mVns": mu_err, "sigma_mVns": sigma, "amp": amp})
    fits = pd.DataFrame(rows).sort_values("mu_mVns").reset_index(drop=True)
    return fits, centers, counts, smooth


def choose_peak_sequence(fits: pd.DataFrame, bias_v: float) -> tuple[list[int], float, float, float, str]:
    if len(fits) < 2:
        return [], np.nan, np.nan, np.nan, ""
    mus = fits["mu_mVns"].to_numpy(float)
    min_spacing = 700.0 if bias_v < 55.0 else 1000.0
    if bias_v >= 56.0:
        min_spacing = 1450.0
    max_spacing = 7000.0
    best: tuple[float, tuple[int, ...], float, float] | None = None
    for count in range(min(6, len(mus)), 1, -1):
        for combo in combinations(range(len(mus)), count):
            vals = mus[list(combo)]
            diffs = np.diff(vals)
            spacing = float(np.mean(diffs))
            if spacing < min_spacing or spacing > max_spacing:
                continue
            residual = float(np.sqrt(np.mean((diffs - spacing) ** 2))) if diffs.size else 0.0
            rel = residual / spacing if spacing else np.inf
            score = rel - 0.02 * count - 0.000001 * (vals[-1] - vals[0])
            if best is None or score < best[0]:
                best = (score, combo, spacing, residual)
        if best is not None:
            break
    if best is None:
        return [], np.nan, np.nan, np.nan, ""
    _, combo, spacing, residual = best
    selected = list(combo)
    diffs = np.diff(mus[selected])
    if diffs.size > 1:
        scatter_unc = float(np.std(diffs, ddof=1) / np.sqrt(diffs.size))
    else:
        scatter_unc = np.nan
    mu_errs = fits.loc[selected, "mu_err_mVns"].to_numpy(float)
    finite_errs = mu_errs[np.isfinite(mu_errs)]
    fit_unc = np.nan
    if finite_errs.size >= 2:
        # Adjacent differences each depend on two peak centroids. Divide by the
        # number of differences because the plotted point is the mean spacing.
        fit_unc = float(np.sqrt(np.sum(finite_errs**2)) / max(1, diffs.size))
    if np.isfinite(scatter_unc) and np.isfinite(fit_unc):
        unc = float(np.hypot(scatter_unc, fit_unc))
    elif np.isfinite(scatter_unc):
        unc = scatter_unc
    elif np.isfinite(fit_unc):
        unc = fit_unc
    else:
        unc = np.nan
    return selected, float(np.mean(diffs)), unc, residual, ";".join(f"{x:.1f}" for x in diffs)


def channel_state(source: Path) -> str:
    return "after channel swap" if "sep2_swapped_channel" in str(source) else "original channels"


def fit_line(points: pd.DataFrame) -> dict[str, float]:
    x = points["bias_V"].to_numpy(float)
    y = points["area_gap_mVns"].to_numpy(float)
    if len(points) > 2:
        coeff, cov = np.polyfit(x, y, 1, cov=True)
    else:
        coeff = np.polyfit(x, y, 1)
        cov = np.full((2, 2), np.nan)
    slope, intercept = coeff
    vbr = -intercept / slope
    dv_dm = intercept / slope**2
    dv_db = -1.0 / slope
    vbr_unc = float(np.sqrt(dv_dm**2 * cov[0, 0] + dv_db**2 * cov[1, 1] + 2.0 * dv_dm * dv_db * cov[0, 1])) if len(points) > 2 else np.nan
    pred = slope * x + intercept
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return {
        "n_fit_points": int(len(points)),
        "slope_mVns_per_V": float(slope),
        "intercept_mVns": float(intercept),
        "breakdown_voltage_V": float(vbr),
        "breakdown_voltage_unc_V": vbr_unc,
        "r2": float(1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan,
        "fit_rms_mVns": float(np.sqrt(np.mean((y - pred) ** 2))),
    }


def analyze_runs() -> tuple[pd.DataFrame, dict[str, dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, pd.DataFrame, list[int]]]]]:
    rows: list[dict[str, object]] = []
    hist: dict[str, dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, pd.DataFrame, list[int]]]] = {"SIPM1": {}, "SIPM2": {}}
    for run in read_manifest():
        values, channel, n_rows, n_accepted = area_values(run)
        if values.size < 50:
            continue
        sat = saturation_count(run.source)
        fits, centers, counts, smooth = fit_area_peak_centers(values, run.bias_v)
        selected, gap, gap_unc, residual, diffs = choose_peak_sequence(fits, run.bias_v)
        hist[run.sipm][str(run.source)] = (centers, counts, smooth, fits, selected)
        use_in_fit = (
            np.isfinite(gap)
            and channel_state(run.source) == "original channels"
            and sat / max(1, n_rows) < 0.02
        )
        rows.append(
            {
                "sipm": run.sipm,
                "bias_V": run.bias_v,
                "scope_channel": channel,
                "channel_state": channel_state(run.source),
                "n_rows": n_rows,
                "n_accepted": n_accepted,
                "saturation_candidates": sat,
                "area_peak_count": int(len(fits)),
                "selected_peak_count": int(len(selected)),
                "selected_area_peak_mVns": ";".join(f"{fits.loc[i, 'mu_mVns']:.2f}" for i in selected) if selected else "",
                "area_gap_mVns": gap,
                "area_gap_unc_mVns": gap_unc,
                "area_gap_residual_rms_mVns": residual,
                "area_gap_diffs_mVns": diffs,
                "use_in_linear_fit": bool(use_in_fit),
                "source": str(run.source),
            }
        )
    return pd.DataFrame(rows).sort_values(["sipm", "bias_V", "source"]).reset_index(drop=True), hist


def write_hist_plots(points: pd.DataFrame, hist: dict[str, dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, pd.DataFrame, list[int]]]]) -> list[Path]:
    paths: list[Path] = []
    for sipm in ("SIPM1", "SIPM2"):
        sub = points[points["sipm"].eq(sipm)].sort_values(["bias_V", "source"])
        if sub.empty:
            continue
        ncols = 3
        nrows = int(math.ceil(len(sub) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(15, max(3.15 * nrows, 4.2)), dpi=160)
        axes_flat = np.asarray(axes).reshape(-1)
        for ax, (_, row) in zip(axes_flat, sub.iterrows()):
            item = hist[sipm].get(str(row["source"]))
            if item is None:
                continue
            centers, counts, smooth, fits, selected = item
            ax.step(centers, counts, where="mid", color="C0", linewidth=0.95, alpha=0.78)
            ax.plot(centers, smooth, color="black", linewidth=1.05, alpha=0.9)
            for _, peak in fits.iterrows():
                ax.axvline(float(peak["mu_mVns"]), color="0.55", linestyle=":", linewidth=0.8, alpha=0.55)
            for idx in selected:
                ax.axvline(float(fits.loc[idx, "mu_mVns"]), color="C3", linewidth=1.15, alpha=0.9)
            suffix = "used" if row["use_in_linear_fit"] else "shown"
            gap = row["area_gap_mVns"]
            title = f"{row['bias_V']:.3g} V, {row['channel_state']}"
            ax.set_title(title, fontsize=9)
            ax.text(
                0.98,
                0.95,
                f"gap {gap:.0f} mV ns ({suffix})" if np.isfinite(gap) else "gap unresolved",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=7,
                bbox=dict(facecolor="white", edgecolor="0.82", alpha=0.86),
            )
            ax.set_xlabel("Integrated area (mV ns)")
            ax.set_ylabel("Events/bin")
            ax.grid(True, alpha=0.24)
        for ax in axes_flat[len(sub) :]:
            ax.axis("off")
        fig.suptitle(f"{sipm}: area spectra used to find p.e. gaps", fontsize=15)
        fig.tight_layout(rect=[0, 0, 1, 0.965])
        out = OUT_DIR / f"{sipm.lower()}_area_gap_histogram_panel.png"
        fig.savefig(out)
        plt.close(fig)
        paths.append(out)
    return paths


def plot_linear(points: pd.DataFrame, fits: pd.DataFrame) -> Path:
    plt.style.use("default")
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6), dpi=170, sharey=True)
    colors = {"SIPM1": "C0", "SIPM2": "C1"}
    for ax, sipm in zip(axes, ("SIPM1", "SIPM2")):
        sub = points[points["sipm"].eq(sipm)].sort_values("bias_V")
        used = sub[sub["use_in_linear_fit"]]
        original_not_used = sub[(sub["channel_state"].eq("original channels")) & (~sub["use_in_linear_fit"]) & sub["area_gap_mVns"].notna()]
        swapped = sub[(sub["channel_state"].eq("after channel swap")) & sub["area_gap_mVns"].notna()]

        if not used.empty:
            ax.errorbar(
                used["bias_V"],
                used["area_gap_mVns"],
                yerr=used["area_gap_unc_mVns"].replace(0, np.nan),
                fmt="o",
                capsize=4,
                color=colors[sipm],
                label="used in fit",
            )
        if not original_not_used.empty:
            ax.errorbar(
                original_not_used["bias_V"],
                original_not_used["area_gap_mVns"],
                yerr=original_not_used["area_gap_unc_mVns"].replace(0, np.nan),
                fmt="s",
                mfc="none",
                capsize=4,
                color="0.45",
                label="shown, not fit",
            )
        if not swapped.empty:
            ax.errorbar(
                swapped["bias_V"],
                swapped["area_gap_mVns"],
                yerr=swapped["area_gap_unc_mVns"].replace(0, np.nan),
                fmt="D",
                mfc="white",
                mec=colors[sipm],
                ecolor=colors[sipm],
                capsize=4,
                color=colors[sipm],
                label="Sep 2 after channel swap",
            )

        fit = fits[fits["sipm"].eq(sipm)]
        if not fit.empty:
            res = fit.iloc[0]
            xmin = min(float(sub["bias_V"].min()), float(res["breakdown_voltage_V"])) - 0.15
            xmax = float(sub["bias_V"].max()) + 0.20
            xline = np.linspace(xmin, xmax, 240)
            yline = float(res["slope_mVns_per_V"]) * xline + float(res["intercept_mVns"])
            ax.plot(
                xline,
                yline,
                color=colors[sipm],
                linewidth=1.6,
                label=(
                    f"fit: Vbr={res['breakdown_voltage_V']:.2f} +/- {res['breakdown_voltage_unc_V']:.2f} V, "
                    f"R2={res['r2']:.3f}"
                ),
            )
            ax.axvline(float(res["breakdown_voltage_V"]), color=colors[sipm], linestyle="--", linewidth=1.1)
        ax.axhline(0, color="0.35", linewidth=1.0)
        ax.set_title(sipm)
        ax.set_xlabel("Bias voltage (V)")
        ax.grid(True, alpha=0.28)
        ax.legend(fontsize=8, loc="upper left")
    axes[0].set_ylabel("1 p.e. integrated-area spacing (mV ns)")
    fig.suptitle("Breakdown voltage from integrated-area p.e. spacing", fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = OUT_DIR / "sipm_area_gap_vs_bias_linear_fit.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    points, hist = analyze_runs()
    points.to_csv(OUT_DIR / "area_gap_points.csv", index=False)
    fit_rows: list[dict[str, object]] = []
    for sipm in ("SIPM1", "SIPM2"):
        fit_points = points[points["sipm"].eq(sipm) & points["use_in_linear_fit"] & points["area_gap_mVns"].notna()]
        if len(fit_points) >= 2:
            result = fit_line(fit_points.sort_values("bias_V"))
            result["sipm"] = sipm
            fit_rows.append(result)
    fit_df = pd.DataFrame(fit_rows)
    fit_df.to_csv(OUT_DIR / "area_gap_linear_fit_results.csv", index=False)
    plot_path = plot_linear(points, fit_df)
    hist_paths = write_hist_plots(points, hist)
    print(OUT_DIR / "area_gap_points.csv")
    print(OUT_DIR / "area_gap_linear_fit_results.csv")
    print(plot_path)
    for path in hist_paths:
        print(path)
    if not fit_df.empty:
        print(fit_df.to_string(index=False))
    print(points[["sipm", "bias_V", "channel_state", "scope_channel", "area_gap_mVns", "area_gap_unc_mVns", "use_in_linear_fit", "saturation_candidates"]].to_string(index=False))


if __name__ == "__main__":
    main()
