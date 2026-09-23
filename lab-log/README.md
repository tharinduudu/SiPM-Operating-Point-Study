# Laboratory log

This log records what changed between experiments and why. It is not a polished replacement for the raw run manifests. The manifests preserve instrument settings; these notes preserve the reasoning that connected one run to the next.

The first entry covers several manual sub-campaigns because the method was changing quickly. A run-level summary is also available in [`data/campaign_catalog.csv`](../data/campaign_catalog.csv).

| Period | Main work | Record |
|---|---|---|
| 2026-08-27 to 2026-09-16 | Manual waveform captures, added voltages, channel swaps, and high-statistics repeats | [Manual waveform period](2026-08-27_to_09-16_manual_waveforms.md) |
| 2026-09-17 to 2026-09-18 | Automated PicoScope acquisition and first full breakdown scans | [Automation and validation](2026-09-17_to_18_automation.md) |
| 2026-09-19 to 2026-09-21 | Repeatability, low-bias limits, and provisional operating points | [Repeatability and operating point](2026-09-19_to_21_repeatability.md) |
| 2026-09-22 | Scintillator zero-event pilot and first Triangle/Star tests | [Scintillator and new SiPMs](2026-09-22_scintillator_and_new_sipms.md) |
| 2026-09-23 | CH3 trigger diagnosis, corrected scans, and sanity-check repeat | [Trigger correction and final repeat](2026-09-23_trigger_fix_and_sanity_check.md) |
| 2026-09-23 | Small-tile reference trigger, four-channel bias setup, and production acquisition | [Small-tile trigger calibration](2026-09-23_small_tile_trigger_calibration.md) |

## How to add a new entry

Each entry should state:

1. the physical devices and their labels;
2. PCB and scope channel mapping;
3. effective-bias range and temperature;
4. trigger, probe, sampling, and event count;
5. what was changed from the previous run;
6. what was observed before fitting;
7. numerical result and uncertainty;
8. decision: accepted, diagnostic only, or rejected;
9. final hardware state;
10. raw-data and processed-output paths.

Unexpected results should be written down before changing the setup. Failed runs are part of the experimental record.
