#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import math
import re
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
OUT_DIR = ROOT / "waveform_all_run_histogram_panels_20260902"


@dataclass(frozen=True)
class RunMeasurement:
    sipm: str
    bias_v: float
    run_name: str
    source: Path
    digest: str
    n_rows: int
    n_accepted: int
    area: np.ndarray
    height: np.ndarray


def file_digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def infer_sipm(path: Path) -> str | None:
    parts = [part.lower() for part in path.parts]
    local_parts = list(reversed(parts[-5:]))
    for part in local_parts:
        if "sipm2_ch2" in part or "analysis_sipm2" in part:
            return "SIPM2"
        if "sipm1_ch0" in part or "analysis_sipm1" in part:
            return "SIPM1"

    text = str(path).lower()
    if re.search(r"(^|[/_-])ch2([/_-]|$)", text):
        return "SIPM2"
    if re.search(r"(^|[/_-])ch0([/_-]|$)", text):
        return "SIPM1"
    if "sipm2" in text:
        return "SIPM2"
    if "sipm1" in text:
        return "SIPM1"
    return None


def infer_bias(path: Path) -> float | None:
    text = str(path)
    # Match forms used in this project: 55p59V, 53p2v, 56.86, 54p1.
    dot = re.search(r"(?<!\d)(5\d)\.(\d{1,3})(?!\d)", text)
    if dot:
        return float(f"{dot.group(1)}.{dot.group(2)}")
    pform = re.search(r"(?i)(?<!\d)(5\d)p(\d{1,3})\s*v?", text)
    if pform:
        whole = pform.group(1)
        frac = pform.group(2)
        return float(f"{whole}.{frac}")
    return None


def display_run_name(path: Path, bias_v: float, digest: str, duplicate_index: int) -> str:
    parts = path.parts
    candidates = [part for part in parts if "brDownVstudy_" in part or re.search(r"202608", part)]
    if candidates:
        name = candidates[-1]
    else:
        name = path.parent.name
    name = name.replace("brDownVstudy_", "").replace("_analysis_20260831", "").replace("_analysis_20260827", "")
    name = name.replace("waveform_plots_", "").replace("analysis_", "")
    label = f"{bias_v:.2f} V"
    if duplicate_index > 1:
        label += f" run {duplicate_index}"
    return f"{label}\n{name}\n{digest[:7]}"


def active_channel(df: pd.DataFrame, sipm: str) -> str:
    if "active_scope_channel" in df.columns:
        values = df["active_scope_channel"].dropna().astype(str)
        if not values.empty and values.iloc[0].strip() in {"A", "B"}:
            return values.iloc[0].strip()
    return "A" if sipm == "SIPM1" else "B"


def extract_values(path: Path, sipm: str) -> tuple[np.ndarray, np.ndarray, int, int]:
    df = pd.read_csv(path)
    channel = active_channel(df, sipm)
    accepted = df["accepted"].fillna(False).astype(bool) if "accepted" in df.columns else np.ones(len(df), dtype=bool)

    area_col = f"{channel}_area_mVns"
    height_col = f"{channel}_peak_mV"
    if "prompt_area_mVns" in df.columns:
        area_col = "prompt_area_mVns"
    if "prompt_peak_mV" in df.columns:
        height_col = "prompt_peak_mV"

    area = pd.to_numeric(df.loc[accepted, area_col], errors="coerce").dropna().to_numpy(float)
    height = pd.to_numeric(df.loc[accepted, height_col], errors="coerce").dropna().to_numpy(float)
    area = area[np.isfinite(area)]
    height = height[np.isfinite(height)]
    return area, height, int(len(df)), int(np.sum(accepted))


def collect_runs() -> list[RunMeasurement]:
    rows: list[RunMeasurement] = []
    seen: set[str] = set()
    for path in sorted(ROOT.glob("**/pulse_measurements.csv")):
        if "unit_mismatch_old" in str(path):
            continue
        sipm = infer_sipm(path)
        bias = infer_bias(path)
        if sipm is None or bias is None:
            continue
        digest = file_digest(path)
        if digest in seen:
            continue
        seen.add(digest)
        area, height, n_rows, n_acc = extract_values(path, sipm)
        if len(area) < 50 or len(height) < 50:
            continue
        rows.append(
            RunMeasurement(
                sipm=sipm,
                bias_v=bias,
                run_name="",
                source=path,
                digest=digest,
                n_rows=n_rows,
                n_accepted=n_acc,
                area=area,
                height=height,
            )
        )

    output: list[RunMeasurement] = []
    for sipm in ("SIPM1", "SIPM2"):
        group = sorted([r for r in rows if r.sipm == sipm], key=lambda r: (r.bias_v, str(r.source)))
        counts_by_bias: dict[float, int] = {}
        for run in group:
            key = round(run.bias_v, 3)
            counts_by_bias[key] = counts_by_bias.get(key, 0) + 1
            name = display_run_name(run.source, run.bias_v, run.digest, counts_by_bias[key])
            output.append(
                RunMeasurement(
                    sipm=run.sipm,
                    bias_v=run.bias_v,
                    run_name=name,
                    source=run.source,
                    digest=run.digest,
                    n_rows=run.n_rows,
                    n_accepted=run.n_accepted,
                    area=run.area,
                    height=run.height,
                )
            )
    return output


def fd_edges(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    lo, hi = np.percentile(x, [0.2, 99.4])
    pad = max((hi - lo) * 0.06, 1.0)
    lo -= pad
    hi += pad
    if hi <= lo:
        hi = lo + 1.0
    q25, q75 = np.percentile(x[(x >= lo) & (x <= hi)], [25, 75])
    width = 2.0 * (q75 - q25) / np.cbrt(max(1, len(x)))
    if not np.isfinite(width) or width <= 0:
        width = (hi - lo) / 100.0
    width = float(np.clip(width, (hi - lo) / 260.0, (hi - lo) / 35.0))
    start = math.floor(lo / width) * width
    stop = math.ceil(hi / width) * width
    return np.arange(start, stop + width, width)


def mark_candidate_peaks(ax: plt.Axes, centers: np.ndarray, counts: np.ndarray, quantity: str) -> None:
    if len(centers) < 10 or np.max(counts) <= 0:
        return
    bin_width = float(np.median(np.diff(centers)))
    smooth = gaussian_filter1d(counts.astype(float), sigma=1.4)
    if quantity == "area":
        min_spacing = max(3, int(round(700.0 / max(bin_width, 1e-9))))
    else:
        min_spacing = max(3, int(round(8.0 / max(bin_width, 1e-9))))
    peaks, _ = find_peaks(smooth, prominence=max(6.0, 0.025 * np.max(smooth)), distance=min_spacing)
    for peak in peaks:
        ax.axvline(centers[peak], color="C3", linewidth=0.9, alpha=0.75)
    ax.plot(centers, smooth, color="black", linewidth=1.0, alpha=0.9)


def plot_panel(runs: list[RunMeasurement], sipm: str, quantity: str) -> Path:
    selected = [r for r in runs if r.sipm == sipm]
    selected.sort(key=lambda r: (r.bias_v, r.run_name))
    n = len(selected)
    ncols = 3
    nrows = int(math.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, max(3.2 * nrows, 4.5)), dpi=160)
    axes_flat = np.asarray(axes).reshape(-1)

    for ax, run in zip(axes_flat, selected):
        values = run.area if quantity == "area" else run.height
        edges = fd_edges(values)
        counts, edges = np.histogram(values, bins=edges)
        centers = 0.5 * (edges[:-1] + edges[1:])
        ax.step(centers, counts, where="mid", color="C0", linewidth=1.0, alpha=0.82)
        mark_candidate_peaks(ax, centers, counts, quantity)
        ax.set_title(run.run_name, fontsize=9)
        ax.set_ylabel("Events/bin")
        ax.grid(True, alpha=0.25)
        ax.text(
            0.98,
            0.95,
            f"accepted {run.n_accepted}/{run.n_rows}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=7,
            bbox=dict(facecolor="white", edgecolor="0.8", alpha=0.85),
        )

    for ax in axes_flat[n:]:
        ax.axis("off")

    if quantity == "area":
        xlabel = "Prompt integrated area (mV ns)"
        suffix = "area"
    else:
        xlabel = "Prompt peak height (mV)"
        suffix = "height"
    for ax in axes_flat[:n]:
        ax.set_xlabel(xlabel)

    fig.suptitle(f"{sipm} all distinct waveform-run {suffix} histograms", fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    out = OUT_DIR / f"{sipm.lower()}_all_distinct_run_{suffix}_histograms.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def write_manifest(runs: list[RunMeasurement]) -> Path:
    rows = []
    for run in sorted(runs, key=lambda r: (r.sipm, r.bias_v, r.run_name)):
        rows.append(
            {
                "sipm": run.sipm,
                "bias_V": run.bias_v,
                "panel_label": run.run_name.replace("\n", " | "),
                "n_rows": run.n_rows,
                "n_accepted": run.n_accepted,
                "sha256": run.digest,
                "source": str(run.source),
            }
        )
    manifest = OUT_DIR / "all_distinct_waveform_runs_manifest.csv"
    pd.DataFrame(rows).to_csv(manifest, index=False)
    return manifest


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    runs = collect_runs()
    manifest = write_manifest(runs)
    outputs = [
        plot_panel(runs, "SIPM1", "area"),
        plot_panel(runs, "SIPM2", "area"),
    ]
    # Height panels are saved too for checking the peak-finding behavior, but the area
    # panels remain the primary output for the area-gap calibration.
    outputs.extend(
        [
            plot_panel(runs, "SIPM1", "height"),
            plot_panel(runs, "SIPM2", "height"),
        ]
    )
    print(manifest)
    for out in outputs:
        print(out)


if __name__ == "__main__":
    main()
