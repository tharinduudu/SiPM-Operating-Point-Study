# Data guide

## What is stored in Git

The [processed](processed) directory contains compact CSV and JSON outputs used in the figures and result tables. These files are small enough to review and version normally.

The repository does not contain the multi-gigabyte waveform campaigns. Copying selected raw files into Git would make the history difficult to use and would separate those files from their manifests. Raw data remain on the Ubuntu acquisition computer under:

```text
/home/muon/brDownVstudy/experiments/
```

The original top-level locations are listed in [raw_data_index.csv](raw_data_index.csv). The more detailed [campaign catalog](campaign_catalog.csv) records the devices, channel mapping, voltage coverage, event scale, purpose, and whether the campaign was accepted, diagnostic, or rejected.

## Processed directories

### `processed/original_pair`

- first automated expert results;
- repeated-scan breakdown values;
- interleaved validation points;
- lower-bias resolution results for SIPM1 and SIPM2.

### `processed/new_triangle_star_pair`

- CH3 trigger scan and corrected bias-scan summaries;
- first clean two-device Vbr result;
- repeated clean result;
- two-run repeatability summary.

### `processed/manual_campaigns`

- first 55.59 V waveform panels, event tables, peak fits, and saturation checks;
- first provisional all-voltage height-gap breakdown fit;
- robust 31 August to 2 September height and area reanalysis;
- 50,000-event high-statistics height, area, and saturation checks.

These files document how the algorithm developed. They are not substituted for the later automated result.

## Units

- voltage: V unless a column explicitly states mV;
- time: ns;
- pulse-area proxy: mV ns at the corrected probe-tip scale;
- temperature: C;
- uncertainties: one standard uncertainty unless stated otherwise.

## Reproducibility note

Every result must be tied to a run manifest and a physical-device mapping. The same label can appear on different PCB or scope channels after a swap. Never infer physical identity from a filename alone when the lab log records a mapping correction.
