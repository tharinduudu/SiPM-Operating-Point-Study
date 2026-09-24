# Figure and result catalog

This is a guide to the evidence in the repository. The purpose is to make it possible to follow the study without opening every directory.

## First manual waveforms

| Figure | What it shows |
|---|---|
| `figures/historical/manual/2026-08-27_sipm1_area_waveform_peak_panel.png` | Pulse-area histogram, waveform family, and height/area behavior from the first campaign |
| `figures/historical/manual/2026-08-27_refined_height_spacing_comparison.png` | Early local-Gaussian p.e. spacing comparison |
| `figures/historical/manual/2026-08-27_initial_height_spacing_vbr.png` | First provisional breakdown fit from peak-height gaps |
| `figures/historical/manual/2026-08-27_sipm1_saturation_check.png` | Full waveform overlay and height/area saturation diagnostic |
| `figures/historical/manual/2026-08-27_sipm2_saturation_check.png` | Same check for the second SiPM |

The full compact output, including CSV peak centers and event measurements, is in `data/processed/manual_campaigns/2026-08-27_55p59V_waveform_analysis/`. The first all-voltage line-fit tables are in `data/processed/manual_campaigns/2026-08-27_initial_breakdown_analysis/`.

## Robust manual reanalysis

| Figure | What it shows |
|---|---|
| `figures/historical/method_evolution/2026-09-01_sipm1_height_spectra.png` | SIPM1 peak-height spectra ordered from low to high voltage |
| `figures/historical/method_evolution/2026-09-01_sipm1_area_spectra.png` | SIPM1 pulse-area spectra and expected gap lines |
| `figures/historical/method_evolution/2026-09-01_sipm2_height_spectra.png` | SIPM2 height spectra with larger high-bias peak distance |
| `figures/historical/method_evolution/2026-09-01_sipm2_area_spectra.png` | SIPM2 area spectra across the available voltages |

The corresponding tables and notes are in `data/processed/manual_campaigns/2026-09-01_*`.

## High-statistics checks

The `figures/historical/high_statistics/` directory contains the 50,000-event height, area, and saturation panels for both original SiPMs at 55.59 V. These plots show why event count alone does not decide whether every shoulder is a new p.e. population.

## Automated-method checks

| Figure | Question answered |
|---|---|
| `analysis_flow.png` | What is the order of the final analysis? |
| `pedestal_characterization.png` | Is the zero-p.e. response stable and narrow enough? |
| `pulse_observable_correlations.png` | Do height and area remain correlated? |
| `acquisition_quality_control.png` | Did acceptance, noise, or pulse timing change with bias? |
| `sipm1_charge_spectra.png`, `sipm2_charge_spectra.png` | Are the fitted p.e. populations visible in the raw bins? |
| `charge_gain_breakdown_fits.png` | What Vbr follows from the primary area gap? |
| `height_gain_breakdown_fits.png` | Does peak height give a compatible intercept? |
| `integration_gate_systematic.png` | How much does the area gate move the result? |
| `leave_one_out_stability.png` | Does one voltage point control the intercept? |
| `digitizer_interpolation_diagnostic.png` | Is peak-height quantization affecting the result? |

These files are in `figures/historical/method_evolution/`.

## Final original-pair results

The accepted automated original-pair tables are in `data/processed/original_pair/`. The central results are the charge-spacing fits, repeated-scan comparisons, lower-bias resolvability tests, and interleaved checks. `docs/05_results.md` explains which numerical result is recommended and which values are retained only as cross-checks.

## SiPM 1 (△) and SiPM 2 (★) results

The final plots are in `figures/results/` and compact data in `data/processed/new_sipm_pair/`:

- CH3 trigger scan and the 25 mV correction;
- first clean two-device spectra and Vbr fit;
- independent repeat spectra and Vbr fit;
- run-to-run reproducibility comparison;
- weighted SiPM 1 (△) and SiPM 2 (★) result.

Rejected and diagnostic scans are still listed in `data/campaign_catalog.csv` and the lab log. They are not silently deleted from the history.

## Reading a result safely

Before quoting a number, check all four items:

1. Is the physical SiPM identity known, including swaps?
2. Was the trigger population valid for that channel?
3. Was the p.e. spacing resolved by the stated acceptance rule?
4. Does the uncertainty include only internal analysis, or also absolute bias calibration?

The final answer is only as reliable as those four statements.
