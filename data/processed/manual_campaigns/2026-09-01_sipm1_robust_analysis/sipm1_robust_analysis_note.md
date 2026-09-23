# SIPM1 robust p.e. peak recalculation

This recalculation uses the raw PicoScope CSV files and does not reuse the older baseline-subtracted measurements.

## Method

1. Each waveform is read in the scope units and converted to mV.
2. The baseline is estimated event-by-event from the quiet pre-trigger region, normally `t < -50 ns`.
3. The baseline estimator is robust: median first, then MAD-based clipping, then a final clipped median. This avoids pulling the baseline when a small pre-pulse or pickup excursion is present.
4. The noise estimate is also robust: `sigma = 1.4826 * MAD` from the clipped pre-trigger samples.
5. The waveform is baseline-subtracted.
6. The prompt signal is measured only in `-5 ns <= t <= 90 ns`. This avoids mistaking late afterpulses or recovery features for the main triggered avalanche.
7. Pulse height is measured from a mildly smoothed waveform using a Savitzky-Golay filter, followed by a 3-point quadratic interpolation around the maximum.
8. Events are accepted when the prompt peak is above `max(5 sigma_baseline, 8 mV)`, not clipped, and the peak time lies between `-5 ns` and `80 ns`.
9. A pulse-height histogram is built for each bias voltage.
10. Candidate p.e. peaks are found from a lightly smoothed histogram and refined with local Gaussian-plus-linear-background fits.
11. The reported 1 p.e. spacing is the average separation of the best nearly-even sequence of fitted peak centers.
12. Points marked diagnostic are plotted and saved, but not used for the breakdown-voltage fit.

## SIPM1 result

- Breakdown voltage from selected robust height-spacing points: `51.448 +/- 0.217 V`
- Height-spacing slope: `12.612 mV/V`
- Fit R^2: `0.9929`
- Fit RMS: `0.751 mV`

## Why this is better than the first-pass method

The first-pass method used a simple mean baseline from `t < -30 ns` and searched for peaks in the resulting histogram. The new method uses a robust baseline and a prompt-only peak gate. That is better for these recorded waveforms because the scope is triggered on pulses, not randomly, and the traces can contain late pulses or small baseline disturbances that should not change the pulse-height estimate.
