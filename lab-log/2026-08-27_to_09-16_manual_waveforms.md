# Manual waveform period: 2026-08-27 to 2026-09-16

## Purpose

The first aim was to see whether the amplified dark pulses contained recognizable photoelectron structure and whether that structure changed with SiPM bias. The PicoScope software was used to save waveform sets at several voltages. Histograms of peak height and integrated area were then compared.

## 27 August recordings

SIPM1 was on PCB CH0 and SIPM2 on PCB CH2. The first series used 10,000 events per recorded device and point.

| Bias | SIPM1 | SIPM2 |
|---:|---|---|
| 53.90 V | recorded | not recorded |
| 54.90 V | recorded | recorded |
| 55.59 V | recorded | recorded |
| 55.65 V | recorded | recorded |
| 56.40 V | recorded | recorded |
| 56.86 V | recorded | recorded |

The early script used the pre-trigger mean, found histogram maxima after Gaussian smoothing, fitted each maximum locally, and selected a nearly even set of peak centers. This produced a first p.e. gap estimate, but it was sensitive to baseline disturbances and to small sub-peaks inside a broad population.

## What was learned

- Distinct pulse populations could be seen over part of the bias range.
- Pulse area and peak height both showed increasing separation with bias.
- Pulse area was less sensitive to changes in waveform shape.
- At higher pulse amplitudes, saturation and non-Gaussian shoulders had to be checked rather than assumed absent.
- At low voltage, changing the scope trigger changed which events entered the histogram.
- A channel swap can invalidate device labels unless the physical SiPM identity is recorded separately from the electronic channel.

The manual data were useful for learning the waveform shape and building the analysis. They were not treated as the final calibration because acquisition settings and metadata were not controlled uniformly across every run.

## Development during this period

Additional 20,000-event voltages were recorded from 2026-08-31 to 2026-09-02, including 53.20, 54.10, and 56.30 V. Channel-switched checks were made near 55.27, 56.45, and 57.867 V. The same physical devices were moved between readout channels to check whether observed differences followed the SiPM or the electronics.

The analysis was restarted from raw CSV rather than reusing the first baseline-subtracted tables. The baseline became a median/MAD clipped estimate. The maximum had to lie in the prompt region, exceed `max(5 x baseline noise, 8 mV)`, and remain away from clipping. This was the first analysis in the study that explicitly separated physical pulse acceptance from finding peaks in a histogram.

On 2026-09-09, 20,000- and 50,000-event captures were made at the important comparison voltages, including 55.59 V. The aim was to inspect small populations and determine whether fitted p.e. spacing was stable with more statistics.

The high-statistics spectra showed that the simple Gaussian-comb description does not reproduce every tail and shoulder. It can still provide a stable period between populations, but the fit covariance by itself is not a complete uncertainty.

The detailed method and figure locations are now recorded in `docs/08_complete_experiment_history.md`, `docs/09_waveform_processing.md`, `docs/10_peak_finding_and_pe_gap.md`, and `docs/12_figure_and_result_catalog.md`.

## Decision

The next scan needed to be automated. Every voltage point had to use the same probe correction, sample interval, waveform window, trigger definition, event count, and file metadata. The raw manual campaigns were retained for comparison rather than mixed into the later automated result.

## Raw records

- `/home/muon/brDownVstudy/experiments/2026-08-27_initial_bias_waveform_series/`
- `/home/muon/brDownVstudy/experiments/2026-08-31_to_09-02_extended_bias_and_channel_swap/`
- `/home/muon/brDownVstudy/experiments/2026-09-09_high_statistics_repeats/`
