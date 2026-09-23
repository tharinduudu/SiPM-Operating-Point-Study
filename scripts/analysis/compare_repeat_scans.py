#!/usr/bin/env python3
"""Compare two independently analyzed SiPM breakdown-voltage scans."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest, chi2


COLORS = {"SIPM1": "#1f77b4", "SIPM2": "#ff7f0e"}
MARKERS = {"original": "o", "repeat": "s"}


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 10.5,
            "axes.titlesize": 12,
            "axes.labelsize": 10.5,
            "legend.fontsize": 9,
            "axes.grid": True,
            "grid.alpha": 0.22,
            "figure.dpi": 170,
            "savefig.dpi": 220,
            "savefig.bbox": "tight",
        }
    )


def load_run(path: Path) -> tuple[dict, dict[str, pd.DataFrame]]:
    result = json.loads((path / "expert_results.json").read_text(encoding="utf-8"))
    tables = {
        sipm: pd.read_csv(path / f"{sipm.lower()}_charge_comb_fits.csv")
        for sipm in ("SIPM1", "SIPM2")
    }
    return result, tables


def matched_table(original: pd.DataFrame, repeat: pd.DataFrame) -> pd.DataFrame:
    left = original.copy()
    right = repeat.copy()
    left["bias_key"] = left["bias_v"].round(3)
    right["bias_key"] = right["bias_v"].round(3)
    merged = left.merge(right, on="bias_key", suffixes=("_original", "_repeat"), validate="one_to_one")
    merged["delta_gap_mV_ns"] = merged["gap_mV_ns_repeat"] - merged["gap_mV_ns_original"]
    merged["combined_gap_unc_mV_ns"] = np.hypot(
        merged["gap_unc_mV_ns_original"], merged["gap_unc_mV_ns_repeat"]
    )
    merged["standardized_delta"] = merged["delta_gap_mV_ns"] / merged["combined_gap_unc_mV_ns"]
    merged["mean_gap_mV_ns"] = 0.5 * (merged["gap_mV_ns_original"] + merged["gap_mV_ns_repeat"])
    merged["fractional_delta_pct"] = 100.0 * merged["delta_gap_mV_ns"] / merged["gap_mV_ns_original"]
    return merged


def weighted_mean(values: np.ndarray, errors: np.ndarray) -> tuple[float, float]:
    weights = 1.0 / np.square(errors)
    return float(np.sum(weights * values) / np.sum(weights)), float(np.sqrt(1.0 / np.sum(weights)))


def run_summary(
    sipm: str,
    original_result: dict,
    repeat_result: dict,
    matched: pd.DataFrame,
) -> dict[str, float | int]:
    valid = matched[matched["fit_success_original"].astype(bool) & matched["fit_success_repeat"].astype(bool)].copy()
    differences = valid["delta_gap_mV_ns"].to_numpy(float)
    errors = valid["combined_gap_unc_mV_ns"].to_numpy(float)
    statistic = float(np.sum(np.square(differences / errors)))
    dof = int(len(valid))
    probability = float(chi2.sf(statistic, dof)) if dof else np.nan
    mean_delta, mean_delta_unc = weighted_mean(differences, errors)
    old = original_result["recommended"][sipm]
    new = repeat_result["recommended"][sipm]
    vbr_delta = float(new["breakdown_voltage_V"] - old["breakdown_voltage_V"])
    vbr_combined_fit_unc = float(np.hypot(new["fit_uncertainty_V"], old["fit_uncertainty_V"]))
    vbr_combined_internal_unc = float(
        np.hypot(new["internal_total_uncertainty_V"], old["internal_total_uncertainty_V"])
    )
    combined_vbr, combined_vbr_fit_unc = weighted_mean(
        np.asarray([old["breakdown_voltage_V"], new["breakdown_voltage_V"]], float),
        np.asarray([old["fit_uncertainty_V"], new["fit_uncertainty_V"]], float),
    )
    positive_differences = int(np.sum(differences > 0))
    negative_differences = int(np.sum(differences < 0))
    sign_test_p = float(binomtest(positive_differences, positive_differences + negative_differences, 0.5).pvalue)
    height_delta = float(new["height_crosscheck_V"] - old["height_crosscheck_V"])
    return {
        "matched_points": dof,
        "gap_difference_chi_squared": statistic,
        "gap_difference_degrees_of_freedom": dof,
        "gap_difference_p_value": probability,
        "weighted_mean_gap_difference_mV_ns": mean_delta,
        "weighted_mean_gap_difference_uncertainty_mV_ns": mean_delta_unc,
        "weighted_mean_gap_difference_z": mean_delta / mean_delta_unc,
        "positive_gap_differences": positive_differences,
        "negative_gap_differences": negative_differences,
        "exploratory_two_sided_sign_test_p_value": sign_test_p,
        "rms_fractional_gap_difference_pct": float(np.sqrt(np.mean(np.square(valid["fractional_delta_pct"])))),
        "maximum_absolute_standardized_gap_difference": float(np.max(np.abs(valid["standardized_delta"]))),
        "original_vbr_V": float(old["breakdown_voltage_V"]),
        "repeat_vbr_V": float(new["breakdown_voltage_V"]),
        "vbr_difference_V": vbr_delta,
        "vbr_combined_fit_uncertainty_V": vbr_combined_fit_unc,
        "vbr_difference_fit_z": vbr_delta / vbr_combined_fit_unc,
        "vbr_combined_internal_uncertainty_V": vbr_combined_internal_unc,
        "vbr_difference_internal_z": vbr_delta / vbr_combined_internal_unc,
        "fit_weighted_combined_vbr_V": combined_vbr,
        "fit_weighted_combined_vbr_uncertainty_V": combined_vbr_fit_unc,
        "height_crosscheck_difference_V": height_delta,
        "original_charge_slope_mV_ns_per_V": float(old["charge_fit_slope_mV_ns_per_V"]),
        "repeat_charge_slope_mV_ns_per_V": float(new["charge_fit_slope_mV_ns_per_V"]),
        "charge_slope_change_pct": 100.0
        * (float(new["charge_fit_slope_mV_ns_per_V"]) / float(old["charge_fit_slope_mV_ns_per_V"]) - 1.0),
    }


def plot_gap_overlay(output: Path, results: dict[str, dict], tables: dict[str, dict[str, pd.DataFrame]]) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(12.8, 5.2))
    for axis, sipm in zip(axes, ("SIPM1", "SIPM2")):
        for label in ("original", "repeat"):
            table = tables[label][sipm]
            valid = table[table["fit_success"].astype(bool)]
            axis.errorbar(
                valid["bias_v"],
                valid["gap_mV_ns"],
                yerr=valid["gap_unc_mV_ns"],
                fmt=MARKERS[label],
                markersize=5.5,
                capsize=3,
                color=COLORS[sipm],
                alpha=1.0 if label == "repeat" else 0.50,
                markerfacecolor=COLORS[sipm] if label == "repeat" else "white",
                label=f"{label} spectra",
            )
            fit = results[label]["fits"][sipm]["charge"]
            x = np.linspace(valid["bias_v"].min() - 0.1, valid["bias_v"].max() + 0.1, 200)
            y = fit["slope"] * (x - fit["vbr"])
            axis.plot(x, y, color=COLORS[sipm], alpha=1.0 if label == "repeat" else 0.45, linestyle="-" if label == "repeat" else "--", label=f"{label} fit: Vbr={fit['vbr']:.3f} V")
        axis.set_title(sipm)
        axis.set_xlabel("Effective bias (V)")
        axis.set_ylabel("Photoelectron charge spacing (mV ns)")
        axis.legend()
    figure.suptitle("Independent scan repeat: charge-spacing response")
    figure.tight_layout()
    figure.savefig(output / "repeatability_charge_gap_overlay.png")
    plt.close(figure)


def plot_differences(output: Path, matched: dict[str, pd.DataFrame]) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(12.8, 8.4))
    for column, sipm in enumerate(("SIPM1", "SIPM2")):
        table = matched[sipm]
        valid = table[table["fit_success_original"].astype(bool) & table["fit_success_repeat"].astype(bool)]
        top = axes[0, column]
        top.errorbar(valid["bias_key"], valid["delta_gap_mV_ns"], yerr=valid["combined_gap_unc_mV_ns"], fmt="o", color=COLORS[sipm], capsize=3)
        top.axhline(0.0, color="black", linewidth=1)
        top.set_title(f"{sipm}: repeat - original")
        top.set_xlabel("Effective bias (V)")
        top.set_ylabel("Charge-spacing difference (mV ns)")
        bottom = axes[1, column]
        bottom.axhspan(-2.0, 2.0, color="#d9ead3", alpha=0.6, label="+/-2 combined sigma")
        bottom.axhline(0.0, color="black", linewidth=1)
        bottom.plot(valid["bias_key"], valid["standardized_delta"], "o-", color=COLORS[sipm])
        bottom.set_xlabel("Effective bias (V)")
        bottom.set_ylabel("Difference / combined fit uncertainty")
        bottom.legend()
    figure.suptitle("Point-by-point repeatability residuals")
    figure.tight_layout()
    figure.savefig(output / "repeatability_point_differences.png")
    plt.close(figure)


def plot_identity_bland_altman(output: Path, matched: dict[str, pd.DataFrame]) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(12.8, 8.2))
    for column, sipm in enumerate(("SIPM1", "SIPM2")):
        table = matched[sipm]
        valid = table[table["fit_success_original"].astype(bool) & table["fit_success_repeat"].astype(bool)]
        x = valid["gap_mV_ns_original"].to_numpy(float)
        y = valid["gap_mV_ns_repeat"].to_numpy(float)
        low = min(x.min(), y.min())
        high = max(x.max(), y.max())
        identity = axes[0, column]
        identity.plot([low, high], [low, high], color="black", linewidth=1, label="identity")
        identity.errorbar(x, y, xerr=valid["gap_unc_mV_ns_original"], yerr=valid["gap_unc_mV_ns_repeat"], fmt="o", color=COLORS[sipm], capsize=3)
        identity.set_title(f"{sipm}: repeat versus original")
        identity.set_xlabel("Original spacing (mV ns)")
        identity.set_ylabel("Repeat spacing (mV ns)")
        identity.legend()
        difference = y - x
        mean = 0.5 * (x + y)
        mean_difference = float(np.mean(difference))
        sd_difference = float(np.std(difference, ddof=1)) if len(difference) > 1 else 0.0
        bland = axes[1, column]
        bland.scatter(mean, difference, color=COLORS[sipm])
        bland.axhline(mean_difference, color="black", label=f"mean={mean_difference:.1f}")
        bland.axhline(mean_difference + 1.96 * sd_difference, color="0.4", linestyle="--", label="95% limits")
        bland.axhline(mean_difference - 1.96 * sd_difference, color="0.4", linestyle="--")
        bland.set_xlabel("Mean spacing of the two scans (mV ns)")
        bland.set_ylabel("Repeat - original (mV ns)")
        bland.legend()
    figure.suptitle("Identity and Bland-Altman repeatability views")
    figure.tight_layout()
    figure.savefig(output / "repeatability_identity_bland_altman.png")
    plt.close(figure)


def plot_vbr_summary(output: Path, results: dict[str, dict]) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(12.3, 5.1), sharey=False)
    for axis, sipm in zip(axes, ("SIPM1", "SIPM2")):
        positions = np.arange(2)
        for metric, marker, offset in (("charge", "o", -0.09), ("height", "s", 0.09)):
            values = [results[label]["fits"][sipm][metric]["vbr"] for label in ("original", "repeat")]
            errors = [results[label]["fits"][sipm][metric]["vbr_unc"] for label in ("original", "repeat")]
            axis.errorbar(positions + offset, values, yerr=errors, fmt=marker, color=COLORS[sipm], capsize=4, label=metric)
        axis.set_xticks(positions, ["Original", "Repeat"])
        axis.set_ylabel("Breakdown voltage at 20 C (V)")
        axis.set_title(sipm)
        axis.legend()
    figure.suptitle("Independent breakdown-voltage estimates")
    figure.tight_layout()
    figure.savefig(output / "repeatability_breakdown_voltage.png")
    plt.close(figure)


def plot_environment(output: Path, original: dict, repeat: dict) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(12.3, 4.8))
    for label, result, marker in (("Original", original, "o"), ("Repeat", repeat, "s")):
        axes[0].plot(result["effective_bias_points_V"], result["temperatures_C"], marker=marker, label=label)
        correction_mv = 1000.0 * (np.asarray(result["effective_bias_points_V"]) - np.asarray(result["temperature_normalized_biases_V"]))
        axes[1].plot(result["effective_bias_points_V"], correction_mv, marker=marker, label=label)
    axes[0].set_xlabel("Effective bias (V)")
    axes[0].set_ylabel("Recorded temperature (C)")
    axes[0].set_title("Temperature during each scan")
    axes[1].set_xlabel("Effective bias (V)")
    axes[1].set_ylabel("20 C bias correction (mV)")
    axes[1].set_title("Applied temperature normalization")
    for axis in axes:
        axis.legend()
    figure.tight_layout()
    figure.savefig(output / "repeatability_environment.png")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", type=Path, required=True, help="original expert_results directory")
    parser.add_argument("--repeat", type=Path, required=True, help="repeat expert_results directory")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    configure_style()

    original_result, original_tables = load_run(args.original.resolve())
    repeat_result, repeat_tables = load_run(args.repeat.resolve())
    results = {"original": original_result, "repeat": repeat_result}
    tables = {"original": original_tables, "repeat": repeat_tables}
    matched = {sipm: matched_table(original_tables[sipm], repeat_tables[sipm]) for sipm in ("SIPM1", "SIPM2")}

    summaries = {
        sipm: run_summary(sipm, original_result, repeat_result, matched[sipm])
        for sipm in ("SIPM1", "SIPM2")
    }
    for sipm, table in matched.items():
        table.to_csv(output / f"{sipm.lower()}_point_comparison.csv", index=False)
    (output / "repeatability_results.json").write_text(json.dumps(summaries, indent=2) + "\n", encoding="utf-8")
    pd.DataFrame([{"sipm": sipm, **values} for sipm, values in summaries.items()]).to_csv(
        output / "repeatability_summary.csv", index=False
    )

    plot_gap_overlay(output, results, tables)
    plot_differences(output, matched)
    plot_identity_bland_altman(output, matched)
    plot_vbr_summary(output, results)
    plot_environment(output, original_result, repeat_result)
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
