#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageChops
from scipy.ndimage import gaussian_filter1d


ROOT = Path("/Users/tharinduhettiarachchi/Desktop/Design/new415")
ANALYSIS_DIR = ROOT / "brDownVstudy_56p86V_20000_analysis_20260831" / "analysis_sipm1_ch0_56p86V_20000"
CALIB_DIR = ROOT / "height_gap_linear_calibration_20260902"
OUT_PATH = CALIB_DIR / "sipm1_56p86_three_panel_hist_waveform_area.png"


def fd_edges(values: np.ndarray, min_bins: int = 45, max_bins: int = 180) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 2:
        return np.linspace(0, 1, 21)
    lo, hi = np.percentile(x, [0.15, 99.35])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(np.nanmin(x)), float(np.nanmax(x))
    pad = 0.05 * (hi - lo)
    lo -= pad
    hi += pad
    q25, q75 = np.percentile(x[(x >= lo) & (x <= hi)], [25, 75])
    width = 2.0 * (q75 - q25) / np.cbrt(x.size)
    if not np.isfinite(width) or width <= 0:
        width = (hi - lo) / 80.0
    bins = int(np.clip(np.ceil((hi - lo) / width), min_bins, max_bins))
    return np.linspace(lo, hi, bins + 1)


def parse_peaks(text: object) -> list[float]:
    if pd.isna(text):
        return []
    peaks: list[float] = []
    for part in str(text).replace(",", ";").split(";"):
        try:
            peaks.append(float(part.strip()))
        except ValueError:
            pass
    return peaks


def read_cropped_plot(path: Path, pad: int = 18) -> np.ndarray:
    img = Image.open(path).convert("RGB")
    white = Image.new("RGB", img.size, (255, 255, 255))
    bbox = ImageChops.difference(img, white).getbbox()
    if bbox is None:
        return np.asarray(img)
    left, top, right, bottom = bbox
    crop = (
        max(0, left - pad),
        max(0, top - pad),
        min(img.width, right + pad),
        min(img.height, bottom + pad),
    )
    return np.asarray(img.crop(crop))


def main() -> None:
    CALIB_DIR.mkdir(parents=True, exist_ok=True)
    meas = pd.read_csv(ANALYSIS_DIR / "pulse_measurements.csv")
    summary = pd.read_csv(ANALYSIS_DIR / "summary.csv").iloc[0]
    rows = pd.read_csv(CALIB_DIR / "clean_height_gap_improved_histogram_runs.csv")
    fit_row = rows[(rows["sipm"].eq("SIPM1")) & (np.isclose(rows["bias_V"], 56.86)) & (rows["group"].eq("fit"))].iloc[0]

    active = str(summary["active_scope_channel"])
    accepted = meas["accepted"].fillna(False).astype(bool)
    heights_mV = pd.to_numeric(meas.loc[accepted, f"{active}_peak_mV"], errors="coerce").dropna().to_numpy(float)
    areas_mVns = pd.to_numeric(meas.loc[accepted, f"{active}_area_mVns"], errors="coerce").dropna().to_numpy(float)
    area_vns = areas_mVns / 1000.0
    height_v = heights_mV / 1000.0
    peaks = parse_peaks(fit_row["selected_height_peak_mV"])

    plt.style.use("default")
    fig, axes = plt.subplots(
        3,
        1,
        figsize=(13.2, 14.6),
        dpi=180,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [1.05, 1.75, 1.10], "hspace": 0.10},
    )
    fig.suptitle("SIPM1 at 56.86 V", fontsize=20)

    ax = axes[0]
    edges = fd_edges(heights_mV)
    counts, edges = np.histogram(heights_mV, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    smooth = gaussian_filter1d(counts.astype(float), sigma=1.15)
    ax.step(centers, counts, where="mid", color="C0", linewidth=1.2, alpha=0.78, label="height histogram")
    ax.plot(centers, smooth, color="black", linewidth=1.5, alpha=0.9, label="smoothed histogram")
    for i, peak in enumerate(peaks):
        ax.axvline(peak, color="C3", linewidth=1.25, alpha=0.9, label="selected p.e. peaks" if i == 0 else None)
    gap = float(fit_row["height_gap_mV"])
    unc = float(fit_row["height_gap_unc_mV"])
    ax.set_xlim(100, 450)
    ax.set_xlabel("Prompt peak height (mV)", fontsize=13)
    ax.set_ylabel("Events/bin", fontsize=13)
    ax.set_title("Peak-height histogram", fontsize=15)
    ax.legend(title=f"1 p.e. spacing = {gap:.1f} +/- {unc:.1f} mV", fontsize=10.5, title_fontsize=11)
    ax.grid(True, alpha=0.25)

    ax = axes[1]
    waveform_img = read_cropped_plot(ANALYSIS_DIR / "all_waveforms_overlay_zoom.png")
    ax.imshow(waveform_img)
    ax.set_axis_off()
    ax.set_title("Baseline-subtracted waveforms", fontsize=15, pad=8)

    ax = axes[2]
    ax.scatter(area_vns, height_v, s=8, alpha=0.18, color="C0", edgecolors="none", rasterized=True)
    finite = np.isfinite(area_vns) & np.isfinite(height_v)
    if finite.sum() > 10:
        coeff = np.polyfit(area_vns[finite], height_v[finite], 1)
        xline = np.linspace(np.nanpercentile(area_vns[finite], 0.2), np.nanpercentile(area_vns[finite], 99.5), 200)
        ax.plot(xline, coeff[0] * xline + coeff[1], color="black", linewidth=1.4, label="linear trend")
    ax.set_xlabel("Integrated area (V ns)", fontsize=13)
    ax.set_ylabel("Peak height (V)", fontsize=13)
    ax.set_title("Peak height vs integrated area", fontsize=15)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=10.5, loc="upper left")

    for ax in (axes[0], axes[2]):
        ax.tick_params(axis="both", labelsize=11)

    fig.savefig(OUT_PATH)
    plt.close(fig)
    print(OUT_PATH)


if __name__ == "__main__":
    main()
