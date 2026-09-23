# Manual waveform period: 2026-08-27 to 2026-09-16

## Purpose

The first aim was to see whether the amplified dark pulses contained recognizable photoelectron structure and whether that structure changed with SiPM bias. The PicoScope software was used to save waveform sets at several voltages. Histograms of peak height and integrated area were then compared.

## What was learned

- Distinct pulse populations could be seen over part of the bias range.
- Pulse area and peak height both showed increasing separation with bias.
- Pulse area was less sensitive to changes in waveform shape.
- At higher pulse amplitudes, saturation and non-Gaussian shoulders had to be checked rather than assumed absent.
- At low voltage, changing the scope trigger changed which events entered the histogram.
- A channel swap can invalidate device labels unless the physical SiPM identity is recorded separately from the electronic channel.

The manual data were useful for learning the waveform shape and building the analysis. They were not treated as the final calibration because acquisition settings and metadata were not controlled uniformly across every run.

## Development during this period

Additional voltages were recorded from 2026-08-31 to 2026-09-02. The same physical devices were also moved between readout channels to check whether observed differences followed the SiPM or the electronics. On 2026-09-09, 20,000- and 50,000-event captures were made to inspect small populations and determine whether fitted p.e. spacing was stable with more statistics.

The high-statistics spectra showed that the simple Gaussian-comb description does not reproduce every tail and shoulder. It can still provide a stable period between populations, but the fit covariance by itself is not a complete uncertainty.

## Decision

The next scan needed to be automated. Every voltage point had to use the same probe correction, sample interval, waveform window, trigger definition, event count, and file metadata. The raw manual campaigns were retained for comparison rather than mixed into the later automated result.

## Raw records

- `/home/muon/brDownVstudy/experiments/2026-08-27_initial_bias_waveform_series/`
- `/home/muon/brDownVstudy/experiments/2026-08-31_to_09-02_extended_bias_and_channel_swap/`
- `/home/muon/brDownVstudy/experiments/2026-09-09_high_statistics_repeats/`
