#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path("/Users/tharinduhettiarachchi/Desktop/Design/new415")
CALIB_DIR = ROOT / "height_gap_linear_calibration_20260902"
POINTS_CSV = CALIB_DIR / "improved_height_gap_fit_points.csv"
OUT_PATH = CALIB_DIR / "checked_vbr_peak_height_line_fit.png"
OUT_CSV = CALIB_DIR / "checked_vbr_peak_height_line_fit_results.csv"

CHANNEL_STYLE = {
    "A": {"color": "C0", "label": "detector CH0 / scope A"},
    "B": {"color": "C1", "label": "detector CH2 / scope B"},
}


def fit_line(points: pd.DataFrame) -> dict[str, float]:
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


def display_yerr(yerr: pd.Series, y: pd.Series) -> np.ndarray:
    err = pd.to_numeric(yerr, errors="coerce").to_numpy(float)
    val = pd.to_numeric(y, errors="coerce").to_numpy(float)
    err[(~np.isfinite(err)) | (err <= 0) | (err > 0.5 * np.abs(val)) | (err > 30.0)] = np.nan
    return err


def main() -> None:
    points = pd.read_csv(POINTS_CSV).sort_values(["sipm", "bias_V"])
    checks = pd.read_csv(CALIB_DIR / "clean_height_gap_improved_histogram_runs.csv")
    checks = checks[checks["group"].eq("Sep 2 check")].copy()
    resolved_points = pd.read_csv(CALIB_DIR / "height_gap_points.csv")
    for sipm, bias in [("SIPM2", 56.86)]:
        resolved = resolved_points[
            resolved_points["sipm"].eq(sipm)
            & np.isclose(resolved_points["bias_V"], bias)
            & resolved_points["channel_state"].eq("original channels")
        ].iloc[0]
        target = points["sipm"].eq(sipm) & np.isclose(points["bias_V"], bias)
        for col in [
            "height_gap_mV",
            "height_gap_unc_mV",
            "selected_height_peak_mV",
            "selected_peak_count",
        ]:
            points.loc[target, col] = resolved[col]

    keep = {
        "SIPM1": points["sipm"].eq("SIPM1"),
        "SIPM2": points["sipm"].eq("SIPM2"),
    }
    results = []
    for sipm, mask in keep.items():
        fit_points = points[mask].copy()
        result = fit_line(fit_points)
        result["sipm"] = sipm
        results.append(result)
    results_df = pd.DataFrame(results)
    results_df.to_csv(OUT_CSV, index=False)

    plt.style.use("default")
    fig, axes = plt.subplots(1, 2, figsize=(16.2, 6.6), dpi=170, sharey=True)
    for ax, sipm in zip(axes, ["SIPM1", "SIPM2"]):
        sub = points[points["sipm"].eq(sipm)].copy()
        fit_points = sub[keep[sipm].loc[sub.index]].copy()
        check_points = checks[checks["sipm"].eq(sipm)].sort_values("bias_V")
        result = results_df[results_df["sipm"].eq(sipm)].iloc[0]
        fit_channel = str(fit_points.iloc[0].get("scope_channel", "A"))
        fit_style = CHANNEL_STYLE.get(fit_channel, CHANNEL_STYLE["A"])

        ax.errorbar(
            fit_points["bias_V"],
            fit_points["height_gap_mV"],
            yerr=display_yerr(fit_points["height_gap_unc_mV"], fit_points["height_gap_mV"]),
            fmt="o",
            color=fit_style["color"],
            capsize=4,
            markersize=5.5,
            linewidth=1.1,
            label=f"fit points: {fit_style['label']}",
            zorder=4,
        )
        if not check_points.empty:
            for channel, channel_points in check_points.groupby("scope_channel"):
                check_style = CHANNEL_STYLE.get(str(channel), CHANNEL_STYLE["A"])
                ax.errorbar(
                    channel_points["bias_V"],
                    channel_points["height_gap_mV"],
                    yerr=display_yerr(channel_points["height_gap_unc_mV"], channel_points["height_gap_mV"]),
                    fmt="D",
                    mfc="white",
                    mec=check_style["color"],
                    ecolor=check_style["color"],
                    color=check_style["color"],
                    capsize=4,
                    markersize=6.2,
                    markeredgewidth=1.4,
                    linewidth=1.1,
                    label=f"channel-switched checks: {check_style['label']}",
                    zorder=5,
                )

        x_min = min(float(result["breakdown_voltage_V"]), float(fit_points["bias_V"].min())) - 0.22
        x_max_values = [float(fit_points["bias_V"].max())]
        if not check_points.empty:
            x_max_values.append(float(check_points["bias_V"].max()))
        x_max = max(x_max_values) + 0.28
        x_line = np.linspace(x_min, x_max, 300)
        y_line = float(result["slope_mV_per_V"]) * x_line + float(result["intercept_mV"])
        ax.plot(
            x_line,
            y_line,
            color=fit_style["color"],
            linewidth=1.9,
            label=(
                f"linear fit: Vbr={result['breakdown_voltage_V']:.2f} +/- "
                f"{result['breakdown_voltage_unc_V']:.2f} V, "
                f"R2={result['r2']:.3f}"
            ),
            zorder=2,
        )
        ax.axvline(float(result["breakdown_voltage_V"]), color=fit_style["color"], linestyle="--", linewidth=1.35)
        ax.axhline(0, color="0.30", linewidth=1.0)
        ax.set_title(sipm, fontsize=15)
        ax.set_xlabel("Bias voltage (V)", fontsize=14)
        ax.grid(True, alpha=0.28)
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(-5, 85)
        ax.tick_params(axis="both", labelsize=12)
        ax.legend(fontsize=11, loc="upper left", framealpha=0.9)

    axes[0].set_ylabel("1 p.e. peak-height spacing (mV)", fontsize=14)
    fig.suptitle("Breakdown voltage from p.e. peak-height spacing", fontsize=21)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(OUT_PATH)
    plt.close(fig)
    print(OUT_PATH)
    print(OUT_CSV)
    print(results_df.to_string(index=False))


if __name__ == "__main__":
    main()
