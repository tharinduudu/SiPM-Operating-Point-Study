#!/usr/bin/env python3
"""Assess the low-overvoltage p.e. resolution without forcing a gain fit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import least_squares

from expert_vbr_analysis import fit_charge_comb, poisson_deviance_residual


REFERENCE = {
    "SIPM1": {
        "slope_mV_ns_per_V": 531.6362382578809,
        "slope_unc_mV_ns_per_V": 5.1613853581116,
        "vbr_V": 50.63024796943835,
        "vbr_unc_V": 0.04814993121301643,
        "intrinsic_scatter_mV_ns": 7.927988606545367,
        "color": "#1f77b4",
    },
    "SIPM2": {
        "slope_mV_ns_per_V": 550.7615796322087,
        "slope_unc_mV_ns_per_V": 11.674113644620045,
        "vbr_V": 51.95649336083294,
        "vbr_unc_V": 0.08380240217701126,
        "intrinsic_scatter_mV_ns": 13.271161932676616,
        "color": "#ff7f0e",
    },
}
TEMPERATURE_COEFFICIENT_V_PER_C = 0.054
REFERENCE_TEMPERATURE_C = 20.0


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "legend.fontsize": 9,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "axes.grid": True,
            "grid.alpha": 0.22,
            "figure.dpi": 170,
            "savefig.dpi": 220,
            "savefig.bbox": "tight",
        }
    )


def point_rows(study: Path, sipm: str) -> list[dict[str, object]]:
    manifest = json.loads((study / "scan_manifest.json").read_text(encoding="utf-8"))
    rows = []
    for point in manifest["points"]:
        tag = str(point["tag"])
        directory = study / "analysis_aligned" / tag / sipm
        summary_path = directory / "point_summary.json"
        measurement_path = directory / "pulse_measurements.csv"
        if summary_path.exists():
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        else:
            summary = {
                "bias_v": float(point["effective_bias_v"]),
                "quality_pass": False,
                "area_quality_pass": False,
                "selected_peak_centers_mV": [],
                "one_pe_spacing_mV": np.nan,
                "one_pe_spacing_unc_mV": np.nan,
                "one_pe_area_spacing_mV_ns": np.nan,
                "one_pe_area_spacing_unc_mV_ns": np.nan,
            }
        measurements = pd.read_csv(measurement_path)
        rows.append(
            {
                "tag": tag,
                "summary": summary,
                "measurements": measurements,
            }
        )
    return rows


def corrected_biases(manifest: dict[str, object]) -> dict[str, tuple[float, float]]:
    result = {}
    for point in manifest["points"]:
        measured = float(point["effective_bias_v"])
        temperature = float(point["environment"]["temperature_C"])
        corrected = measured - TEMPERATURE_COEFFICIENT_V_PER_C * (
            temperature - REFERENCE_TEMPERATURE_C
        )
        result[str(point["tag"])] = (measured, corrected)
    return result


def null_model_delta_bic(fit: object) -> float:
    """Compare the p.e. comb with one Gaussian plus smooth background."""
    x = fit.centers
    counts = fit.counts
    x_min = float(x[0])
    span = max(float(x[-1] - x[0]), 1.0)
    peak_index = int(np.argmax(counts))
    maximum = max(float(np.max(counts)), 1.0)
    parameters0 = np.asarray(
        [
            float(x[peak_index]),
            max(float(fit.pedestal_sigma), span / 12.0, 10.0),
            maximum,
            max(0.1, float(np.percentile(counts, 10))),
            max(0.1, float(np.percentile(counts, 25))),
            max(span / 3.0, 50.0),
        ]
    )
    lower = np.asarray([x_min, 5.0, 0.0, 0.0, 0.0, 10.0])
    upper = np.asarray([float(x[-1]), span, 4.0 * maximum, maximum, 4.0 * maximum, 4.0 * span])

    def residual(parameters: np.ndarray) -> np.ndarray:
        mean, sigma, amplitude, background, tail_amplitude, tail_scale = parameters
        model = (
            background
            + tail_amplitude * np.exp(-(x - x_min) / tail_scale)
            + amplitude * np.exp(-0.5 * np.square((x - mean) / sigma))
        )
        return poisson_deviance_residual(counts, model)

    result = least_squares(
        residual,
        parameters0,
        bounds=(lower, upper),
        x_scale="jac",
        max_nfev=30000,
    )
    null_deviance = float(np.sum(np.square(residual(result.x))))
    comb_parameters = 6 + len(fit.pe_indexes)
    comb_degrees_of_freedom = max(len(counts) - comb_parameters, 1)
    comb_deviance = float(fit.reduced_deviance * comb_degrees_of_freedom)
    sample_size = max(len(counts), 2)
    null_bic = null_deviance + 6.0 * np.log(sample_size)
    comb_bic = comb_deviance + comb_parameters * np.log(sample_size)
    return float(null_bic - comb_bic)


def fit_points(
    study: Path,
    sipm: str,
    pedestal: np.ndarray,
    bias_map: dict[str, tuple[float, float]],
) -> tuple[pd.DataFrame, dict[str, object]]:
    reference = REFERENCE[sipm]
    records = point_rows(study, sipm)
    fits: dict[str, object] = {}
    output_rows = []
    for record in records:
        tag = str(record["tag"])
        summary = record["summary"]
        table = record["measurements"]
        measured_bias, corrected_bias = bias_map[tag]
        expected_gap = float(reference["slope_mV_ns_per_V"]) * (
            corrected_bias - float(reference["vbr_V"])
        )
        values = table.loc[
            table["accepted"].fillna(False), "pulse_area_mV_ns"
        ].dropna().to_numpy(float)
        fit = None
        reasons = []
        delta_bic = np.nan
        multistart_relative_spread = np.nan
        if expected_gap <= 0:
            reasons.append("no positive resolvable spacing expected")
        elif len(values) < 200:
            reasons.append("fewer than 200 events passed the five-sigma pulse cut")
        else:
            fit = fit_charge_comb(values, max(expected_gap, 180.0), pedestal)
            fits[tag] = fit
            trial_fits = [
                fit_charge_comb(values, max(scale * expected_gap, 180.0), pedestal)
                for scale in (0.90, 1.00, 1.10)
            ]
            stable_gaps = np.asarray(
                [trial.gap for trial in trial_fits if trial.success], dtype=float
            )
            if len(stable_gaps) >= 2:
                multistart_relative_spread = float(
                    (np.max(stable_gaps) - np.min(stable_gaps)) / np.median(stable_gaps)
                )
            delta_bic = null_model_delta_bic(fit)
            amplitude_threshold = max(0.04 * float(np.max(fit.amplitudes)), 1.0)
            visible = np.flatnonzero(
                (fit.pe_indexes >= 1) & (fit.amplitudes >= amplitude_threshold)
            )
            if len(visible) >= 2:
                first, second = visible[:2]
                resolution = float(
                    (fit.means[second] - fit.means[first])
                    / np.hypot(fit.sigmas[first], fit.sigmas[second])
                )
            else:
                resolution = np.nan
            prior_agreement = bool(abs(fit.gap - expected_gap) <= 0.22 * expected_gap)
            away_from_bound = bool(0.75 * expected_gap < fit.gap < 1.25 * expected_gap)
            if not fit.success:
                reasons.append("constrained comb fit failed")
            if not np.isfinite(delta_bic) or delta_bic < 10.0:
                reasons.append("multiple populations are not preferred by BIC")
            if not np.isfinite(multistart_relative_spread) or multistart_relative_spread > 0.12:
                reasons.append("comb spacing is not stable across initial values")
            if not prior_agreement:
                reasons.append("spacing is inconsistent with the frozen higher-bias line")
            if not away_from_bound:
                reasons.append("comb spacing reached a fit boundary")
            if not np.isfinite(resolution) or resolution < 1.25:
                reasons.append("adjacent populations are not sufficiently separated")
        validated = fit is not None and not reasons
        output_rows.append(
            {
                "sipm": sipm,
                "tag": tag,
                "commanded_bias_v": float(summary["bias_v"]),
                "measured_bias_v": measured_bias,
                "temperature_corrected_bias_v": corrected_bias,
                "temperature_c": REFERENCE_TEMPERATURE_C
                + (measured_bias - corrected_bias) / TEMPERATURE_COEFFICIENT_V_PER_C,
                "accepted_events": int(table["accepted"].fillna(False).sum()),
                "total_events": int(len(table)),
                "expected_gap_mV_ns": expected_gap,
                "fitted_gap_mV_ns": float(fit.gap) if fit is not None else np.nan,
                "fitted_gap_unc_mV_ns": float(fit.gap_unc) if fit is not None else np.nan,
                "reference_residual_mV_ns": (
                    float(fit.gap - expected_gap) if fit is not None else np.nan
                ),
                "reference_residual_pct": (
                    float(100.0 * (fit.gap - expected_gap) / expected_gap)
                    if fit is not None and expected_gap > 0
                    else np.nan
                ),
                "reduced_deviance": float(fit.reduced_deviance) if fit is not None else np.nan,
                "delta_bic_comb_over_single": delta_bic,
                "multistart_relative_spread": multistart_relative_spread,
                "validated": validated,
                "reason": "resolved" if validated else "; ".join(reasons),
            }
        )
    return pd.DataFrame(output_rows), fits


def plot_charge_spectra(
    study: Path,
    output: Path,
    sipm: str,
    table: pd.DataFrame,
    fits: dict[str, object],
) -> None:
    records = point_rows(study, sipm)
    all_values = []
    by_tag = {}
    for record in records:
        values = record["measurements"].loc[
            record["measurements"]["accepted"].fillna(False), "pulse_area_mV_ns"
        ].dropna().to_numpy(float)
        by_tag[str(record["tag"])] = values
        all_values.append(values)
    combined = np.concatenate(all_values)
    x_low = min(0.0, float(np.percentile(combined, 0.1)))
    x_high = float(np.percentile(combined, 99.7))
    figure, axes = plt.subplots(2, 4, figsize=(17.0, 8.4), sharex=True)
    color = str(REFERENCE[sipm]["color"])
    for axis, record in zip(axes.flat, records):
        tag = str(record["tag"])
        values = by_tag[tag]
        row = table.loc[table["tag"] == tag].iloc[0]
        fit = fits.get(tag)
        if fit is None:
            bins = "fd" if len(values) >= 2 else 50
            counts, edges = np.histogram(values, bins=bins, range=(x_low, x_high))
            centers = 0.5 * (edges[:-1] + edges[1:])
            axis.step(centers, counts, where="mid", color=color, linewidth=1.0)
        else:
            axis.step(fit.centers, fit.counts, where="mid", color=color, linewidth=1.0, label="measured bins")
            axis.plot(
                fit.centers,
                fit.model,
                color="black" if bool(row["validated"]) else "0.45",
                linestyle="-" if bool(row["validated"]) else "--",
                linewidth=1.5,
                label="comb fit" if bool(row["validated"]) else "diagnostic fit",
            )
            if bool(row["validated"]):
                threshold = max(0.04 * float(np.max(fit.amplitudes)), 1.0)
                for pe_index, mean, amplitude in zip(fit.pe_indexes, fit.means, fit.amplitudes):
                    if pe_index >= 1 and amplitude >= threshold:
                        axis.axvline(mean, color="0.35", linestyle=":", linewidth=0.9)
        axis.set_yscale("log")
        axis.set_ylim(bottom=0.8)
        axis.set_xlim(x_low, x_high)
        axis.set_title(f"{float(row['measured_bias_v']):.3f} V")
        if bool(row["validated"]):
            label = f"resolved: {row['fitted_gap_mV_ns']:.0f} +/- {row['fitted_gap_unc_mV_ns']:.0f} mV ns"
        else:
            label = "p.e. spacing not resolved"
        axis.legend([label], loc="upper right")
    for axis in axes[-1, :]:
        axis.set_xlabel("Pulse area (mV ns)")
    for axis in axes[:, 0]:
        axis.set_ylabel("Events per bin")
    figure.suptitle(f"{sipm}: lower-bias pulse-area spectra", fontsize=16)
    figure.tight_layout()
    figure.savefig(output / f"{sipm.lower()}_lower_bias_charge_spectra.png")
    plt.close(figure)


def plot_height_spectra(study: Path, output: Path, sipm: str) -> None:
    records = point_rows(study, sipm)
    values_by_tag = {}
    all_values = []
    for record in records:
        values = record["measurements"].loc[
            record["measurements"]["accepted"].fillna(False), "peak_height_mV"
        ].dropna().to_numpy(float)
        values_by_tag[str(record["tag"])] = values
        all_values.append(values)
    combined = np.concatenate(all_values)
    x_low, x_high = np.percentile(combined, [0.1, 99.7])
    figure, axes = plt.subplots(2, 4, figsize=(17.0, 8.4), sharex=True)
    color = str(REFERENCE[sipm]["color"])
    for axis, record in zip(axes.flat, records):
        tag = str(record["tag"])
        summary = record["summary"]
        values = values_by_tag[tag]
        bins = "fd" if len(values) >= 2 else 50
        counts, edges = np.histogram(values, bins=bins, range=(x_low, x_high))
        centers = 0.5 * (edges[:-1] + edges[1:])
        width = float(np.mean(np.diff(centers)))
        smooth = gaussian_filter1d(counts.astype(float), sigma=max(1.0, 4.0 / width))
        axis.step(centers, counts, where="mid", color=color, linewidth=1.0, label="measured bins")
        axis.plot(centers, smooth, color="black", linewidth=1.4, label="4 mV Gaussian smoothing")
        if bool(summary.get("quality_pass", False)):
            for center in summary.get("selected_peak_centers_mV", []):
                axis.axvline(center, color="0.35", linestyle=":", linewidth=0.9)
            status = f"resolved: {summary['one_pe_spacing_mV']:.1f} +/- {summary['one_pe_spacing_unc_mV']:.1f} mV"
        else:
            status = "p.e. spacing not resolved"
        axis.set_yscale("log")
        axis.set_ylim(bottom=0.8)
        axis.set_xlim(x_low, x_high)
        axis.set_title(f"{float(summary['bias_v']):.3f} V")
        axis.legend([status], loc="upper right")
    for axis in axes[-1, :]:
        axis.set_xlabel("Pulse height (mV)")
    for axis in axes[:, 0]:
        axis.set_ylabel("Events per bin")
    figure.suptitle(f"{sipm}: lower-bias pulse-height spectra", fontsize=16)
    figure.tight_layout()
    figure.savefig(output / f"{sipm.lower()}_lower_bias_height_spectra.png")
    plt.close(figure)


def plot_resolution_summary(output: Path, tables: dict[str, pd.DataFrame]) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(13.4, 5.8), sharex=True, sharey=True)
    for axis, sipm in zip(axes, ("SIPM1", "SIPM2")):
        reference = REFERENCE[sipm]
        table = tables[sipm]
        color = str(reference["color"])
        fit_x = np.linspace(float(reference["vbr_V"]), 53.55, 300)
        fit_y = float(reference["slope_mV_ns_per_V"]) * (fit_x - float(reference["vbr_V"]))
        axis.plot(fit_x, fit_y, color=color, linewidth=2.0, label="higher-bias reference fit")
        valid = table["validated"]
        if valid.any():
            axis.errorbar(
                table.loc[valid, "temperature_corrected_bias_v"],
                table.loc[valid, "fitted_gap_mV_ns"],
                yerr=table.loc[valid, "fitted_gap_unc_mV_ns"],
                fmt="o",
                color=color,
                markeredgecolor="black",
                capsize=4,
                label="statistically resolved points",
            )
        rejected = ~valid
        if rejected.any():
            axis.scatter(
                table.loc[rejected, "temperature_corrected_bias_v"],
                np.zeros(int(rejected.sum())),
                marker="x",
                s=65,
                color="0.25",
                label="spacing not resolved",
                zorder=4,
            )
        axis.axhline(0, color="0.3", linewidth=1.0)
        axis.axvline(float(reference["vbr_V"]), color=color, linestyle="--", linewidth=1.2)
        axis.set_title(
            f"{sipm}\nreference Vbr = {reference['vbr_V']:.3f} +/- {reference['vbr_unc_V']:.3f} V"
        )
        axis.set_xlabel("Bias voltage corrected to 20 C (V)")
        axis.legend(loc="upper left")
    axes[0].set_ylabel("One-photoelectron charge spacing (mV ns)")
    axes[0].set_xlim(50.4, 53.55)
    axes[0].set_ylim(-80, None)
    figure.suptitle("Lower-bias resolution test against the frozen gain lines", fontsize=15)
    figure.tight_layout()
    figure.savefig(output / "lower_bias_resolution_summary.png")
    plt.close(figure)


def plot_reference_residuals(output: Path, tables: dict[str, pd.DataFrame]) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(13.4, 5.5), sharex=True, sharey=True)
    for axis, sipm in zip(axes, ("SIPM1", "SIPM2")):
        table = tables[sipm]
        color = str(REFERENCE[sipm]["color"])
        fitted = np.isfinite(table["fitted_gap_mV_ns"])
        valid = fitted & table["validated"]
        diagnostic = fitted & ~table["validated"]
        if valid.any():
            errors = (
                100.0
                * table.loc[valid, "fitted_gap_unc_mV_ns"]
                / table.loc[valid, "expected_gap_mV_ns"]
            )
            axis.errorbar(
                table.loc[valid, "temperature_corrected_bias_v"],
                table.loc[valid, "reference_residual_pct"],
                yerr=errors,
                fmt="o",
                color=color,
                markeredgecolor="black",
                capsize=4,
                label="statistically resolved spacing",
            )
        if diagnostic.any():
            axis.scatter(
                table.loc[diagnostic, "temperature_corrected_bias_v"],
                table.loc[diagnostic, "reference_residual_pct"],
                marker="x",
                s=60,
                color="0.35",
                label="diagnostic fit not accepted",
            )
        axis.axhline(0, color="black", linewidth=1.1)
        axis.axhspan(-5, 5, color="0.75", alpha=0.22, label="within 5%")
        axis.set_title(sipm)
        axis.set_xlabel("Bias voltage corrected to 20 C (V)")
        axis.legend(loc="best")
    axes[0].set_ylabel("Difference from higher-bias gain line (%)")
    axes[0].set_xlim(51.9, 53.5)
    figure.suptitle("Low-bias p.e. spacing relative to the frozen higher-bias fit", fontsize=15)
    figure.tight_layout()
    figure.savefig(output / "lower_bias_reference_residuals.png")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    configure_style()
    args.output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((args.study / "scan_manifest.json").read_text(encoding="utf-8"))
    bias_map = corrected_biases(manifest)
    tables = {}
    result = {
        "reference_source": "independent 2026-09-21 higher-bias repeat",
        "reference": REFERENCE,
        "validation_rule": (
            "A lower-bias charge spacing is accepted only when the pedestal-constrained comb fit succeeds, "
            "the comb is strongly preferred over a one-Gaussian smooth-background model (delta BIC >= 10), "
            "the fitted spacing is stable across three starting values, adjacent populations are separated, "
            "and the result is not on a fit boundary or inconsistent with the frozen reference line."
        ),
        "sipms": {},
    }
    for sipm in ("SIPM1", "SIPM2"):
        pedestal_parts = []
        for label in ("pedestal_pre_aligned", "pedestal_post_aligned"):
            pedestal_table = pd.read_csv(
                args.study / label / sipm / "pedestal_measurements.csv"
            )
            pedestal_parts.append(pedestal_table["pulse_area_mV_ns"].dropna().to_numpy(float))
        pedestal = np.concatenate(pedestal_parts)
        table, fits = fit_points(args.study, sipm, pedestal, bias_map)
        tables[sipm] = table
        table.to_csv(args.output / f"{sipm.lower()}_lower_bias_results.csv", index=False)
        plot_charge_spectra(args.study, args.output, sipm, table, fits)
        plot_height_spectra(args.study, args.output, sipm)
        result["sipms"][sipm] = {
            "resolved_points": int(table["validated"].sum()),
            "tested_points": int(len(table)),
            "lowest_resolved_bias_v": (
                float(table.loc[table["validated"], "measured_bias_v"].min())
                if table["validated"].any()
                else None
            ),
            "rows": table.to_dict(orient="records"),
        }
    plot_resolution_summary(args.output, tables)
    plot_reference_residuals(args.output, tables)
    (args.output / "lower_bias_results.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
