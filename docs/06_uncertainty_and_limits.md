# Uncertainty and limits

## What is included in the quoted error

The reported internal Vbr uncertainty includes:

- uncertainty of the p.e.-spacing fit at each bias point;
- slope-intercept covariance in the linear extrapolation;
- additional point scatter where needed;
- run-to-run weighting for the combined Triangle/Star values.

These terms answer: **if the same calibrated voltages were delivered, how precisely does this waveform method locate the intercept?**

## What is not included

### Absolute bias calibration

The MAX1932 command is write-only. The effective bias is calculated from calibration relations for the high side and channel DAC. It is not continuously measured at the SiPM terminals. A traceable high-impedance measurement is required at every scan point before the internal Vbr error can be called an absolute voltage uncertainty.

### Individual temperature coefficient

The 54 mV/C correction is a provisional device-family value. Individual devices can differ. A controlled temperature scan should fit both `Vbr(T)` and the operating-point behavior.

### Trigger selection

The recorded dark-pulse spectrum is conditional on the PicoScope trigger. It is not an unbiased sample of all avalanches below the trigger. This matters most near breakdown, where 1 p.e. pulses are small and neighboring populations overlap.

### Full-chain response

The measured p.e. spacing includes the analogue gain. It is a good relative gain proxy while the readout chain remains unchanged. It is not an absolute SiPM gain in electrons until the electronics transfer function is calibrated.

### PDE and DCR

Equal overvoltage usually brings nominally similar SiPMs closer to equal gain and PDE, but it does not guarantee equality. PDE also depends on microcell properties and wavelength. DCR depends strongly on device defects and temperature. Optical crosstalk and afterpulsing can differ as well.

The present self-triggered acquisition does not measure physical DCR because its event throughput and dead time are not calibrated. The count of saved waveforms per second is an acquisition rate, not automatically an avalanche rate.

## Statistical cautions found during this study

1. A clean-looking line can still be biased when trigger acceptance is poor.
2. A high R2 does not validate a small set of selected points.
3. A converged multi-peak fit can be unstable when peaks overlap.
4. A very small covariance error can understate run-to-run variation.
5. Repeating a scan does not remove a common voltage-calibration error.
6. A scope threshold that fixes acquisition is not necessarily the best detector comparator threshold.

## Present confidence statement

The evidence supports an internal, repeatable breakdown-voltage comparison for the Triangle and Star SiPMs. It does not yet support a traceable absolute Vbr at the 0.1 V level, nor does it prove matched photon sensitivity. The central operating settings are suitable for the next controlled comparison, provided their provisional status is retained.
