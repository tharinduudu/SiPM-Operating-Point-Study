#!/usr/bin/env python3
"""Produce pedestal-constrained SiPM gain fits and expert diagnostic plots."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import brentq, least_squares
from scipy.stats import norm


TEMPERATURE_COEFFICIENT_V_PER_C = 0.054
REFERENCE_TEMPERATURE_C = 20.0
COLORS = {"SIPM1": "#1f77b4", "SIPM2": "#ff7f0e"}


@dataclass
class CombFit:
    gap: float
    gap_unc: float
    offset: float
    sigma_cell: float
    reduced_deviance: float
    pedestal_mean: float
    pedestal_sigma: float
    pe_indexes: np.ndarray
    means: np.ndarray
    sigmas: np.ndarray
    amplitudes: np.ndarray
    centers: np.ndarray
    counts: np.ndarray
    model: np.ndarray
    components: np.ndarray
    visible_pe: int
    success: bool


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 10.5,
            "axes.titlesize": 12,
            "axes.labelsize": 10.5,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "axes.grid": True,
            "grid.alpha": 0.22,
            "figure.dpi": 170,
            "savefig.dpi": 220,
            "savefig.bbox": "tight",
        }
    )


def poisson_deviance_residual(observed: np.ndarray, expected: np.ndarray) -> np.ndarray:
    expected = np.maximum(expected, 1e-9)
    term = expected - observed
    positive = observed > 0
    term[positive] += observed[positive] * np.log(observed[positive] / expected[positive])
    return np.sign(observed - expected) * np.sqrt(2.0 * np.maximum(term, 0.0))


def comb_values(
    x: np.ndarray,
    parameters: np.ndarray,
    pedestal_mean: float,
    pedestal_sigma: float,
    pe_indexes: np.ndarray,
    x_min: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    gap, offset, sigma_cell, background, tail_amplitude, tail_scale = parameters[:6]
    amplitudes = parameters[6:]
    means = pedestal_mean + offset + pe_indexes * gap
    sigmas = np.sqrt(pedestal_sigma**2 + pe_indexes * sigma_cell**2)
    components = np.asarray(
        [amplitude * np.exp(-0.5 * np.square((x - mean) / sigma)) for amplitude, mean, sigma in zip(amplitudes, means, sigmas)]
    )
    tail = tail_amplitude * np.exp(-(x - x_min) / tail_scale)
    model = background + tail + np.sum(components, axis=0)
    return model, components, means, sigmas


def fit_charge_comb(values: np.ndarray, gap_initial: float, pedestal: np.ndarray) -> CombFit:
    pedestal_mean = float(np.mean(pedestal))
    pedestal_sigma = float(np.std(pedestal, ddof=1))
    low, high = np.percentile(values, [0.15, 99.65])
    low = min(low, pedestal_mean + 0.45 * gap_initial)
    bin_width = max(10.0, min(gap_initial / 45.0, pedestal_sigma / 5.0))
    bins = max(70, int(np.ceil((high - low) / bin_width)))
    counts, edges = np.histogram(values, bins=bins, range=(low, high))
    centers = 0.5 * (edges[:-1] + edges[1:])
    maximum_pe = int(np.clip(np.ceil((high - pedestal_mean) / gap_initial) + 1, 4, 7))
    pe_indexes = np.arange(0, maximum_pe + 1, dtype=float)

    amplitudes = []
    for pe_index in pe_indexes:
        location = pedestal_mean + pe_index * gap_initial
        nearest = int(np.clip(np.searchsorted(centers, location), 0, len(centers) - 1))
        left = max(0, nearest - 2)
        right = min(len(counts), nearest + 3)
        amplitudes.append(max(float(np.max(counts[left:right])) if right > left else 0.0, 0.5))
    parameters0 = np.asarray(
        [
            gap_initial,
            0.0,
            max(0.10 * gap_initial, 25.0),
            max(0.2, float(np.percentile(counts, 10))),
            max(0.2, float(np.percentile(counts, 25))),
            max(gap_initial, 300.0),
            *amplitudes,
        ],
        dtype=float,
    )
    maximum_count = max(float(np.max(counts)), 1.0)
    lower = np.asarray(
        [
            0.72 * gap_initial,
            -0.18 * gap_initial,
            max(10.0, 0.015 * gap_initial),
            0.0,
            0.0,
            0.25 * gap_initial,
            *([0.0] * len(pe_indexes)),
        ]
    )
    upper = np.asarray(
        [
            1.28 * gap_initial,
            0.18 * gap_initial,
            0.45 * gap_initial,
            maximum_count,
            3.0 * maximum_count,
            12.0 * gap_initial,
            *([3.0 * maximum_count] * len(pe_indexes)),
        ]
    )

    def residual(parameters: np.ndarray) -> np.ndarray:
        model, _, _, _ = comb_values(
            centers, parameters, pedestal_mean, pedestal_sigma, pe_indexes, centers[0]
        )
        return poisson_deviance_residual(counts.astype(float), model)

    result = least_squares(
        residual,
        parameters0,
        bounds=(lower, upper),
        x_scale="jac",
        max_nfev=30000,
    )
    model, components, means, sigmas = comb_values(
        centers, result.x, pedestal_mean, pedestal_sigma, pe_indexes, centers[0]
    )
    degrees_of_freedom = max(len(counts) - len(result.x), 1)
    reduced_deviance = float(np.sum(np.square(residual(result.x))) / degrees_of_freedom)
    try:
        covariance = np.linalg.pinv(result.jac.T @ result.jac, rcond=1e-12)
        covariance *= max(reduced_deviance, 1.0)
        gap_unc = float(np.sqrt(max(covariance[0, 0], 0.0)))
    except np.linalg.LinAlgError:
        gap_unc = np.nan
    fitted_amplitudes = result.x[6:]
    amplitude_threshold = max(0.04 * float(np.max(fitted_amplitudes)), 1.0)
    visible_candidates = pe_indexes[(pe_indexes >= 1) & (fitted_amplitudes >= amplitude_threshold)]
    visible_pe = int(visible_candidates[0]) if len(visible_candidates) else -1
    success = bool(
        result.success
        and np.isfinite(gap_unc)
        and gap_unc > 0
        and gap_unc < 0.12 * result.x[0]
        and len(visible_candidates) >= 2
        and 0.35 <= reduced_deviance <= 20.0
    )
    return CombFit(
        gap=float(result.x[0]),
        gap_unc=gap_unc,
        offset=float(result.x[1]),
        sigma_cell=float(result.x[2]),
        reduced_deviance=reduced_deviance,
        pedestal_mean=pedestal_mean,
        pedestal_sigma=pedestal_sigma,
        pe_indexes=pe_indexes,
        means=means,
        sigmas=sigmas,
        amplitudes=fitted_amplitudes,
        centers=centers,
        counts=counts.astype(float),
        model=model,
        components=components,
        visible_pe=visible_pe,
        success=success,
    )


def linear_fit(x: np.ndarray, y: np.ndarray, yerr: np.ndarray, minimum_error: float) -> dict[str, object]:
    yerr = np.maximum(np.asarray(yerr, float), minimum_error)

    def solve(scatter: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
        errors = np.sqrt(np.square(yerr) + scatter**2)
        weights = 1.0 / np.square(errors)
        x_reference = float(np.average(x, weights=weights))
        centered = x - x_reference
        design = np.column_stack([centered, np.ones(len(x))])
        covariance_centered = np.linalg.inv(design.T @ (weights[:, None] * design))
        parameters_centered = covariance_centered @ (design.T @ (weights * y))
        slope, center_intercept = parameters_centered
        intercept = center_intercept - slope * x_reference
        transform = np.asarray([[1.0, 0.0], [-x_reference, 1.0]])
        covariance = transform @ covariance_centered @ transform.T
        residuals = y - (slope * x + intercept)
        reduced_chi_squared = float(np.sum(np.square(residuals / errors)) / max(len(x) - 2, 1))
        return np.asarray([slope, intercept]), covariance, residuals, reduced_chi_squared

    _, _, _, initial_reduced = solve(0.0)
    intrinsic_scatter = 0.0
    if initial_reduced > 1.0 and len(x) > 2:
        def objective(scatter: float) -> float:
            return solve(scatter)[3] - 1.0

        high = max(float(np.std(y)), minimum_error)
        while objective(high) > 0 and high < 10.0 * max(float(np.ptp(y)), 1.0):
            high *= 2.0
        intrinsic_scatter = float(brentq(objective, 0.0, high))
    parameters, covariance, residuals, reduced_chi_squared = solve(intrinsic_scatter)
    slope, intercept = parameters
    vbr = -intercept / slope
    derivative = np.asarray([intercept / slope**2, -1.0 / slope])
    vbr_unc = float(np.sqrt(max(derivative @ covariance @ derivative, 0.0)))
    predicted = slope * x + intercept
    ss_total = float(np.sum(np.square(y - np.mean(y))))
    r_squared = 1.0 - float(np.sum(np.square(residuals))) / ss_total if ss_total > 0 else np.nan
    return {
        "slope": float(slope),
        "slope_unc": float(np.sqrt(max(covariance[0, 0], 0.0))),
        "intercept": float(intercept),
        "vbr": float(vbr),
        "vbr_unc": vbr_unc,
        "intrinsic_scatter": intrinsic_scatter,
        "reduced_chi_squared": reduced_chi_squared,
        "r_squared": r_squared,
        "residuals": residuals,
        "predicted": predicted,
        "effective_errors": np.sqrt(np.square(yerr) + intrinsic_scatter**2),
    }


def leave_one_out(x: np.ndarray, y: np.ndarray, yerr: np.ndarray, minimum_error: float) -> np.ndarray:
    values = []
    for omitted in range(len(x)):
        keep = np.arange(len(x)) != omitted
        values.append(linear_fit(x[keep], y[keep], yerr[keep], minimum_error)["vbr"])
    return np.asarray(values, float)


def point_records(study: Path, sipm: str, analysis_name: str) -> list[dict[str, object]]:
    records = []
    for summary_path in sorted((study / analysis_name).glob(f"bias_*/{sipm}/point_summary.json")):
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        table = pd.read_csv(summary_path.parent / "pulse_measurements.csv")
        records.append(
            {
                "tag": summary_path.parents[1].name,
                "summary": summary,
                "measurements": table,
                "directory": summary_path.parent,
            }
        )
    return records


def save_figure(figure: plt.Figure, path: Path) -> None:
    figure.savefig(path)
    plt.close(figure)


def plot_pedestals(study: Path, output: Path) -> dict[str, dict[str, float]]:
    figure, axes = plt.subplots(2, 2, figsize=(12.0, 8.2))
    summaries: dict[str, dict[str, float]] = {}
    for column, sipm in enumerate(("SIPM1", "SIPM2")):
        color = COLORS[sipm]
        table = pd.read_csv(study / "pedestal_aligned" / sipm / "pedestal_measurements.csv")
        area = table["pulse_area_mV_ns"].to_numpy(float)
        noise = table["baseline_noise_mV"].to_numpy(float)
        area_mean = float(np.mean(area))
        area_sigma = float(np.std(area, ddof=1))
        summaries[sipm] = {
            "area_mean_mV_ns": area_mean,
            "area_sigma_mV_ns": area_sigma,
            "baseline_noise_median_mV": float(np.median(noise)),
        }
        low, high = np.percentile(area, [0.2, 99.8])
        bins = np.linspace(low, high, 90)
        counts, edges, _ = axes[0, column].hist(area, bins=bins, color=color, alpha=0.72)
        centers = 0.5 * (edges[:-1] + edges[1:])
        expected = len(area) * np.mean(np.diff(edges)) * norm.pdf(centers, area_mean, area_sigma)
        axes[0, column].plot(centers, expected, color="black", linewidth=1.7, label="Gaussian reference")
        axes[0, column].axvline(area_mean, color="black", linestyle="--", linewidth=1.1)
        axes[0, column].set_title(f"{sipm}: HV-off charge pedestal")
        axes[0, column].set_xlabel("Integrated baseline charge (mV ns)")
        axes[0, column].set_ylabel("Events per bin")
        axes[0, column].legend(title=f"mean = {area_mean:.1f}\nsigma = {area_sigma:.1f} mV ns")
        axes[1, column].hist(noise, bins="fd", color=color, alpha=0.72)
        axes[1, column].axvline(np.median(noise), color="black", linestyle="--", linewidth=1.1)
        axes[1, column].set_title(f"{sipm}: event baseline noise")
        axes[1, column].set_xlabel("Pretrigger RMS noise (mV)")
        axes[1, column].set_ylabel("Events per bin")
    figure.suptitle("Matched random-trigger pedestal characterization", fontsize=15)
    figure.tight_layout()
    save_figure(figure, output / "pedestal_characterization.png")
    return summaries


def fit_and_plot_charge_spectra(
    study: Path, output: Path, sipm: str, pedestal_values: np.ndarray
) -> tuple[pd.DataFrame, dict[str, CombFit]]:
    records = point_records(study, sipm, "analysis_aligned")
    figure, axes = plt.subplots(2, 3, figsize=(15.2, 8.8))
    rows = []
    fits: dict[str, CombFit] = {}
    for axis, record in zip(axes.flat, records):
        summary = record["summary"]
        table = record["measurements"]
        accepted = table.loc[table["accepted"].fillna(False), "pulse_area_mV_ns"].dropna().to_numpy(float)
        gap_initial = float(summary["one_pe_area_spacing_mV_ns"])
        if not np.isfinite(gap_initial):
            height_gap = float(summary.get("one_pe_spacing_mV", np.nan))
            gap_initial = max(700.0, 50.0 * height_gap) if np.isfinite(height_gap) else 900.0
        fit = fit_charge_comb(accepted, gap_initial, pedestal_values)
        fits[str(record["tag"])] = fit
        accepted_fit = bool(fit.success and summary.get("area_quality_pass", False))
        rows.append(
            {
                "sipm": sipm,
                "tag": record["tag"],
                "bias_v": float(summary["bias_v"]),
                "events": int(len(table)),
                "accepted_events": int(table["accepted"].fillna(False).sum()),
                "overflow_events": int(table.get("overflow", pd.Series(False, index=table.index)).fillna(False).sum()),
                "gap_mV_ns": fit.gap,
                "gap_unc_mV_ns": fit.gap_unc,
                "first_visible_pe": fit.visible_pe,
                "comb_offset_mV_ns": fit.offset,
                "reduced_deviance": fit.reduced_deviance,
                "fit_success": accepted_fit,
            }
        )
        axis.step(fit.centers, fit.counts, where="mid", color=COLORS[sipm], linewidth=1.0, label="data")
        axis.plot(fit.centers, fit.model, color="black", linewidth=1.6, label="constrained fit")
        if accepted_fit:
            for pe_index, component, mean, amplitude in zip(
                fit.pe_indexes, fit.components, fit.means, fit.amplitudes
            ):
                if pe_index >= 1 and amplitude >= 0.04 * max(fit.amplitudes):
                    axis.plot(fit.centers, component, color="0.55", linewidth=0.8, alpha=0.7)
                    axis.axvline(mean, color="0.45", linestyle=":", linewidth=0.8)
                    axis.text(mean, 0.94 * axis.get_ylim()[1], f"{int(pe_index)} p.e.", rotation=90, va="top", ha="right", fontsize=8)
        axis.set_title(f"{float(summary['bias_v']):.3f} V")
        axis.set_xlabel("Peak-aligned pulse charge (mV ns)")
        axis.set_ylabel("Events per bin")
        legend_title = (
            f"G = {fit.gap:.0f} +/- {fit.gap_unc:.0f} mV ns\nD/ndf = {fit.reduced_deviance:.2f}"
            if accepted_fit
            else "excluded: one resolved population"
        )
        axis.legend(title=legend_title)
    figure.suptitle(f"{sipm}: pedestal-constrained multi-photoelectron charge spectra", fontsize=15)
    figure.tight_layout()
    save_figure(figure, output / f"{sipm.lower()}_charge_spectra.png")
    return pd.DataFrame(rows), fits


def plot_height_spectra(study: Path, output: Path, sipm: str) -> None:
    records = point_records(study, sipm, "analysis_aligned")
    figure, axes = plt.subplots(2, 3, figsize=(15.2, 8.8))
    for axis, record in zip(axes.flat, records):
        summary = record["summary"]
        table = record["measurements"]
        values = table.loc[table["accepted"].fillna(False), "peak_height_mV"].dropna().to_numpy(float)
        low, high = np.percentile(values, [0.2, 99.7])
        axis.hist(values[(values >= low) & (values <= high)], bins="fd", histtype="step", color=COLORS[sipm], linewidth=1.1)
        for index, center in enumerate(summary.get("selected_peak_centers_mV", []), start=1):
            axis.axvline(center, color="black", linestyle=":" if index > 1 else "--", linewidth=1.0)
        spacing = summary.get("one_pe_spacing_mV")
        uncertainty = summary.get("one_pe_spacing_unc_mV")
        axis.set_title(f"{float(summary['bias_v']):.3f} V")
        axis.set_xlabel("Quadratic-interpolated pulse height (mV)")
        axis.set_ylabel("Events per bin")
        if isinstance(spacing, (int, float)) and np.isfinite(spacing):
            axis.legend([f"spacing = {spacing:.2f} +/- {uncertainty:.2f} mV"])
        else:
            axis.legend(["one p.e. population only"])
    figure.suptitle(f"{sipm}: pulse-height spectra and resolved p.e. populations", fontsize=15)
    figure.tight_layout()
    save_figure(figure, output / f"{sipm.lower()}_height_spectra.png")


def temperature_corrected_bias(manifest: dict[str, object]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    biases = np.asarray([point["effective_bias_v"] for point in manifest["points"]], float)
    temperatures = np.asarray([point["environment"]["temperature_C"] for point in manifest["points"]], float)
    corrected = biases - TEMPERATURE_COEFFICIENT_V_PER_C * (temperatures - REFERENCE_TEMPERATURE_C)
    return biases, temperatures, corrected


def metric_table_from_summaries(study: Path, sipm: str, analysis_name: str, metric: str) -> pd.DataFrame:
    rows = []
    for record in point_records(study, sipm, analysis_name):
        summary = record["summary"]
        if metric == "height":
            value = summary.get("one_pe_spacing_mV")
            uncertainty = summary.get("one_pe_spacing_unc_mV")
            passed = summary.get("quality_pass")
        else:
            value = summary.get("one_pe_area_spacing_mV_ns")
            uncertainty = summary.get("one_pe_area_spacing_unc_mV_ns")
            passed = summary.get("area_quality_pass")
        if passed and isinstance(value, (int, float)) and np.isfinite(value):
            rows.append({"bias_v": float(summary["bias_v"]), "value": float(value), "uncertainty": float(uncertainty)})
    return pd.DataFrame(rows)


def fit_metric(table: pd.DataFrame, corrected_bias_map: dict[float, float], minimum_error: float) -> dict[str, object]:
    x = np.asarray([corrected_bias_map[round(value, 6)] for value in table["bias_v"]], float)
    y = table["value"].to_numpy(float)
    yerr = table["uncertainty"].to_numpy(float)
    result = linear_fit(x, y, yerr, minimum_error)
    result["x"] = x
    result["y"] = y
    result["yerr"] = yerr
    result["leave_one_out_vbr"] = leave_one_out(x, y, yerr, minimum_error)
    return result


def plot_breakdown_fits(output: Path, fits: dict[str, dict[str, dict[str, object]]], metric: str) -> None:
    unit = "mV ns" if metric == "charge" else "mV"
    label = "1 p.e. charge spacing" if metric == "charge" else "1 p.e. peak-height spacing"
    figure, axes = plt.subplots(2, 2, figsize=(12.8, 8.8), gridspec_kw={"height_ratios": [3.0, 1.0]}, sharex="col")
    for column, sipm in enumerate(("SIPM1", "SIPM2")):
        result = fits[sipm][metric]
        x = result["x"]
        y = result["y"]
        error = result["effective_errors"]
        fit_x = np.linspace(result["vbr"] - 0.15, max(x) + 0.18, 300)
        axes[0, column].errorbar(x, y, yerr=error, fmt="o", color=COLORS[sipm], capsize=4, markersize=6.5, label="accepted bias points")
        axes[0, column].plot(fit_x, result["slope"] * fit_x + result["intercept"], color=COLORS[sipm], linewidth=2.0, label=(f"Vbr(20 C) = {result['vbr']:.3f} +/- {result['vbr_unc']:.3f} V\n" f"slope = {result['slope']:.1f} +/- {result['slope_unc']:.1f} {unit}/V\n" f"R2 = {result['r_squared']:.4f}"))
        axes[0, column].axhline(0, color="0.35", linewidth=1.0)
        axes[0, column].axvline(result["vbr"], color=COLORS[sipm], linestyle="--", linewidth=1.1)
        axes[0, column].set_title(sipm)
        axes[0, column].set_ylabel(f"{label} ({unit})")
        axes[0, column].legend()
        axes[1, column].axhline(0, color="0.35", linewidth=1.0)
        axes[1, column].errorbar(x, result["residuals"], yerr=error, fmt="o", color=COLORS[sipm], capsize=3)
        axes[1, column].set_xlabel("Bias normalized to 20 C (V)")
        axes[1, column].set_ylabel(f"Residual ({unit})")
    figure.suptitle(f"Breakdown voltage from {label}", fontsize=15)
    figure.tight_layout()
    save_figure(figure, output / f"{metric}_gain_breakdown_fits.png")


def plot_quality_control(study: Path, output: Path, manifest: dict[str, object], comb_tables: dict[str, pd.DataFrame]) -> None:
    biases, temperatures, corrected = temperature_corrected_bias(manifest)
    figure, axes = plt.subplots(2, 2, figsize=(12.5, 8.5))
    for sipm in ("SIPM1", "SIPM2"):
        color = COLORS[sipm]
        table = comb_tables[sipm]
        axes[0, 0].plot(table["bias_v"], 100.0 * table["accepted_events"] / table["events"], "o-", color=color, label=sipm)
        axes[0, 1].plot(table["bias_v"], 100.0 * table["overflow_events"] / table["events"], "o-", color=color, label=sipm)
        noises = []
        for record in point_records(study, sipm, "analysis_aligned"):
            accepted = record["measurements"]["accepted"].fillna(False)
            noises.append(record["measurements"].loc[accepted, "baseline_noise_mV"].median())
        axes[1, 0].plot(table["bias_v"], noises, "o-", color=color, label=sipm)
    axes[0, 0].set_title("Waveform acceptance")
    axes[0, 0].set_ylabel("Accepted events (%)")
    axes[0, 0].set_xlabel("Effective bias (V)")
    axes[0, 1].set_title("ADC overflow rejection")
    axes[0, 1].set_ylabel("Overflow events (%)")
    axes[0, 1].set_xlabel("Effective bias (V)")
    axes[1, 0].set_title("Pretrigger baseline noise")
    axes[1, 0].set_ylabel("Median RMS noise (mV)")
    axes[1, 0].set_xlabel("Effective bias (V)")
    axes[1, 1].plot(biases, temperatures, "o-", color="#2ca02c", label="measured temperature")
    axes[1, 1].set_title("Temperature during sequential scan")
    axes[1, 1].set_ylabel("Temperature (C)")
    axes[1, 1].set_xlabel("Effective bias (V)")
    secondary = axes[1, 1].twinx()
    secondary.plot(biases, 1000.0 * (biases - corrected), "s--", color="#9467bd", label="bias correction")
    secondary.set_ylabel("20 C normalization correction (mV)")
    for axis in axes.flat:
        handles, labels = axis.get_legend_handles_labels()
        if handles:
            axis.legend()
    figure.suptitle("Acquisition and environmental quality control", fontsize=15)
    figure.tight_layout()
    save_figure(figure, output / "acquisition_quality_control.png")


def plot_observable_correlations(study: Path, output: Path) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(12.5, 5.2))
    for axis, sipm in zip(axes, ("SIPM1", "SIPM2")):
        path = study / "analysis_aligned" / "bias_55p500" / sipm / "pulse_measurements.csv"
        table = pd.read_csv(path)
        table = table[table["accepted"].fillna(False)]
        artist = axis.hexbin(table["peak_height_mV"], table["pulse_area_mV_ns"], C=table["peak_time_ns"], reduce_C_function=np.mean, gridsize=65, mincnt=1, cmap="viridis")
        correlation = table[["peak_height_mV", "pulse_area_mV_ns"]].corr().iloc[0, 1]
        axis.set_title(f"{sipm} at 55.500 V: r = {correlation:.3f}")
        axis.set_xlabel("Pulse height (mV)")
        axis.set_ylabel("Peak-aligned charge (mV ns)")
        colorbar = figure.colorbar(artist, ax=axis)
        colorbar.set_label("Mean pulse-peak time (ns)")
    figure.suptitle("Pulse-height, charge and trigger-time relationship", fontsize=15)
    figure.tight_layout()
    save_figure(figure, output / "pulse_observable_correlations.png")


def plot_digitizer_diagnostics(study: Path, output: Path) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(12.5, 5.0))
    for axis, sipm in zip(axes, ("SIPM1", "SIPM2")):
        path = study / "analysis_aligned" / "bias_55p500" / sipm / "pulse_measurements.csv"
        table = pd.read_csv(path)
        table = table[table["accepted"].fillna(False)]
        sampled = table["sampled_peak_height_mV"].to_numpy(float)
        interpolated = table["peak_height_mV"].to_numpy(float)
        low, high = np.percentile(interpolated, [0.2, 99.5])
        bins = np.linspace(low, high, 180)
        axis.hist(sampled, bins=bins, histtype="step", linewidth=1.0, color="0.55", label="largest ADC sample")
        axis.hist(interpolated, bins=bins, histtype="step", linewidth=1.4, color=COLORS[sipm], label="three-sample quadratic maximum")
        axis.set_title(f"{sipm} at 55.500 V")
        axis.set_xlabel("Pulse height (mV)")
        axis.set_ylabel("Events per bin")
        axis.legend()
    figure.suptitle("Digitizer-level comb and sub-sample peak interpolation", fontsize=15)
    figure.tight_layout()
    save_figure(figure, output / "digitizer_interpolation_diagnostic.png")


def plot_gate_systematic(study: Path, output: Path, corrected_bias_map: dict[float, float]) -> dict[str, dict[str, float]]:
    definitions = [
        ("Fixed -5 to 100 ns", "analysis_refined"),
        ("Fixed -5 to 250 ns", "analysis_expert"),
        ("Peak -20 to +200 ns", "analysis_aligned"),
    ]
    figure, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), sharey=False)
    results: dict[str, dict[str, float]] = {"SIPM1": {}, "SIPM2": {}}
    markers = ["o", "s", "D"]
    for axis, sipm in zip(axes, ("SIPM1", "SIPM2")):
        for (label, directory), marker in zip(definitions, markers):
            table = metric_table_from_summaries(study, sipm, directory, "area")
            fit = fit_metric(table, corrected_bias_map, 10.0)
            results[sipm][label] = float(fit["vbr"])
            axis.errorbar(fit["x"], fit["y"], yerr=fit["effective_errors"], fmt=marker, capsize=3, label=f"{label}: Vbr={fit['vbr']:.3f} V")
        axis.set_title(sipm)
        axis.set_xlabel("Bias normalized to 20 C (V)")
        axis.set_ylabel("1 p.e. charge spacing (mV ns)")
        axis.legend()
    figure.suptitle("Integration-window sensitivity", fontsize=15)
    figure.tight_layout()
    save_figure(figure, output / "integration_gate_systematic.png")
    return results


def plot_leave_one_out(output: Path, fits: dict[str, dict[str, dict[str, object]]]) -> None:
    figure, axis = plt.subplots(figsize=(9.8, 5.4))
    positions = np.arange(2)
    width = 0.34
    for offset, metric, marker in ((-width / 2, "charge", "o"), (width / 2, "height", "s")):
        for index, sipm in enumerate(("SIPM1", "SIPM2")):
            result = fits[sipm][metric]
            values = result["leave_one_out_vbr"]
            x = np.full(len(values), positions[index] + offset)
            axis.scatter(x, values, color=COLORS[sipm], marker=marker, alpha=0.7, label=f"{metric} omissions" if index == 0 else None)
            axis.errorbar(positions[index] + offset, result["vbr"], yerr=result["vbr_unc"], fmt=marker, color="black", capsize=4)
    axis.set_xticks(positions, ["SIPM1", "SIPM2"])
    axis.set_ylabel("Fitted Vbr at 20 C (V)")
    axis.set_title("Leave-one-bias-point-out stability")
    axis.legend()
    save_figure(figure, output / "leave_one_out_stability.png")


def plot_vbr_summary(output: Path, summary_rows: list[dict[str, object]]) -> None:
    table = pd.DataFrame(summary_rows)
    figure, axis = plt.subplots(figsize=(10.5, 5.6))
    y_positions = np.arange(len(table))[::-1]
    for position, row in zip(y_positions, table.itertuples(index=False)):
        axis.errorbar(row.vbr_V, position, xerr=row.uncertainty_V, fmt=row.marker, color=row.color, capsize=4, markersize=7)
    axis.set_yticks(y_positions, table["label"])
    axis.set_xlabel("Breakdown voltage at 20 C (V)")
    axis.set_title("Breakdown-voltage estimates and internal uncertainty")
    axis.grid(True, axis="x", alpha=0.25)
    axis.grid(False, axis="y")
    save_figure(figure, output / "breakdown_voltage_summary.png")


def plot_analysis_flow(output: Path) -> None:
    figure, axis = plt.subplots(figsize=(13.0, 3.4))
    axis.axis("off")
    labels = [
        "Raw waveform\n(1 ns sampling, 10:1 probe)",
        "Sigma-clipped pretrigger mean\nsubtract event baseline",
        "Peak-aligned charge gate\n-20 ns to +200 ns",
        "Pedestal-constrained comb\nQn = Q0 + nG",
        "Linear G vs bias fit\nintercept gives Vbr",
    ]
    x_positions = np.linspace(0.08, 0.92, len(labels))
    for index, (x, label) in enumerate(zip(x_positions, labels)):
        axis.text(x, 0.52, label, ha="center", va="center", transform=axis.transAxes, bbox={"boxstyle": "round,pad=0.45", "facecolor": "white", "edgecolor": "0.25"})
        if index < len(labels) - 1:
            axis.annotate("", xy=(x_positions[index + 1] - 0.09, 0.52), xytext=(x + 0.09, 0.52), xycoords=axis.transAxes, arrowprops={"arrowstyle": "->", "linewidth": 1.7, "color": "0.25"})
    axis.set_title("Waveform-to-breakdown-voltage analysis chain", fontsize=15)
    save_figure(figure, output / "analysis_flow.png")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    study = args.study.resolve()
    output = (args.output or study / "expert_results").resolve()
    plots = output / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    configure_style()

    manifest = json.loads((study / "scan_manifest.json").read_text(encoding="utf-8"))
    biases, temperatures, corrected_biases = temperature_corrected_bias(manifest)
    corrected_bias_map = {round(float(bias), 6): float(corrected) for bias, corrected in zip(biases, corrected_biases)}

    pedestal_summaries = plot_pedestals(study, plots)
    comb_tables: dict[str, pd.DataFrame] = {}
    for sipm in ("SIPM1", "SIPM2"):
        pedestal = pd.read_csv(study / "pedestal_aligned" / sipm / "pedestal_measurements.csv")["pulse_area_mV_ns"].to_numpy(float)
        comb_table, _ = fit_and_plot_charge_spectra(study, plots, sipm, pedestal)
        comb_tables[sipm] = comb_table
        comb_table.to_csv(output / f"{sipm.lower()}_charge_comb_fits.csv", index=False)
        plot_height_spectra(study, plots, sipm)

    fits: dict[str, dict[str, dict[str, object]]] = {}
    for sipm in ("SIPM1", "SIPM2"):
        fits[sipm] = {}
        charge_table = comb_tables[sipm].loc[comb_tables[sipm]["fit_success"]].copy()
        charge_metric = pd.DataFrame({"bias_v": charge_table["bias_v"], "value": charge_table["gap_mV_ns"], "uncertainty": charge_table["gap_unc_mV_ns"]})
        fits[sipm]["charge"] = fit_metric(charge_metric, corrected_bias_map, 10.0)
        height_table = metric_table_from_summaries(study, sipm, "analysis_aligned", "height")
        fits[sipm]["height"] = fit_metric(height_table, corrected_bias_map, 0.25)

    plot_breakdown_fits(plots, fits, "charge")
    plot_breakdown_fits(plots, fits, "height")
    gate_results = plot_gate_systematic(study, plots, corrected_bias_map)
    plot_quality_control(study, plots, manifest, comb_tables)
    plot_observable_correlations(study, plots)
    plot_digitizer_diagnostics(study, plots)
    plot_leave_one_out(plots, fits)
    plot_analysis_flow(plots)

    summary_rows = []
    final_results: dict[str, dict[str, object]] = {}
    for sipm in ("SIPM1", "SIPM2"):
        charge = fits[sipm]["charge"]
        height = fits[sipm]["height"]
        method_systematic = abs(float(charge["vbr"]) - float(height["vbr"]))
        gate_systematic = abs(gate_results[sipm]["Peak -20 to +200 ns"] - gate_results[sipm]["Fixed -5 to 250 ns"])
        internal_total = float(np.sqrt(charge["vbr_unc"] ** 2 + method_systematic**2 + gate_systematic**2))
        final_results[sipm] = {
            "recommended_method": "pedestal-constrained peak-aligned charge spacing",
            "reference_temperature_C": REFERENCE_TEMPERATURE_C,
            "temperature_coefficient_V_per_C": TEMPERATURE_COEFFICIENT_V_PER_C,
            "breakdown_voltage_V": charge["vbr"],
            "fit_uncertainty_V": charge["vbr_unc"],
            "height_crosscheck_V": height["vbr"],
            "height_charge_difference_V": method_systematic,
            "long_gate_systematic_V": gate_systematic,
            "internal_total_uncertainty_V": internal_total,
            "charge_fit_slope_mV_ns_per_V": charge["slope"],
            "height_fit_slope_mV_per_V": height["slope"],
            "charge_r_squared": charge["r_squared"],
            "height_r_squared": height["r_squared"],
            "leave_one_out_charge_min_V": float(np.min(charge["leave_one_out_vbr"])),
            "leave_one_out_charge_max_V": float(np.max(charge["leave_one_out_vbr"])),
        }
        summary_rows.extend(
            [
                {"label": f"{sipm} charge (primary)", "vbr_V": charge["vbr"], "uncertainty_V": charge["vbr_unc"], "marker": "o", "color": COLORS[sipm]},
                {"label": f"{sipm} height (cross-check)", "vbr_V": height["vbr"], "uncertainty_V": height["vbr_unc"], "marker": "s", "color": COLORS[sipm]},
                {"label": f"{sipm} recommended (internal total)", "vbr_V": charge["vbr"], "uncertainty_V": internal_total, "marker": "D", "color": "black"},
            ]
        )
    plot_vbr_summary(plots, summary_rows)

    serializable_fits = {}
    for sipm, metric_fits in fits.items():
        serializable_fits[sipm] = {}
        for metric, result in metric_fits.items():
            serializable_fits[sipm][metric] = {
                key: value.tolist() if isinstance(value, np.ndarray) else value
                for key, value in result.items()
                if key not in {"predicted", "residuals", "effective_errors", "x", "y", "yerr"}
            }
    results = {
        "scan_root": study.name,
        "events_per_bias_per_sipm": 10000,
        "effective_bias_points_V": biases.tolist(),
        "temperatures_C": temperatures.tolist(),
        "temperature_normalized_biases_V": corrected_biases.tolist(),
        "pedestal": pedestal_summaries,
        "fits": serializable_fits,
        "gate_systematic": gate_results,
        "recommended": final_results,
    }
    (output / "expert_results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    pd.DataFrame(
        [
            {
                "sipm": sipm,
                **values,
            }
            for sipm, values in final_results.items()
        ]
    ).to_csv(output / "breakdown_voltage_results.csv", index=False)
    print(json.dumps(final_results, indent=2))


if __name__ == "__main__":
    main()
