# Trigger correction and sanity check: 2026-09-23

## CH3 trigger scan

The Triangle SiPM was held at an effective bias of 56.401 V on PCB CH3. Five thousand events were requested at each PicoScope trigger threshold from 10 to 45 mV.

Acceptance was 1.06% at 10 mV and 1.62% at 15 mV. It rose to 86.44% at 20 mV and 99.82% at 25 mV, remaining near 99.9% through 45 mV. The recovered p.e. spacing was stable across the high-acceptance region.

**Decision:** use 25 mV for the waveform calibration. This was the lowest tested point on the stable plateau. It is not a change to the detector comparator threshold.

## Corrected CH3 scan

Seven effective-bias points from 54.799 to 56.999 V were recorded with 10,000 events per point. Combined acceptance was 99.806%. All seven area spectra passed the quality requirements.

Result:

```text
Vbr = 50.524 +/- 0.114 V
slope = 519.25 +/- 11.93 mV ns/V
R2 = 0.99892
```

This agreed with the earlier clean-channel Triangle measurement within 0.35 combined standard deviations.

## First matched two-device scan

The Triangle and Star devices were measured with the same 25 mV trigger, bias list, sampling, and event count.

| SiPM | Vbr | Acceptance |
|---|---:|---:|
| Triangle | 50.534 +/- 0.113 V | 99.794% |
| Star | 50.667 +/- 0.097 V | 99.930% |

Their difference was 0.90 sigma and was not significant.

## Independent sanity-check repeat

The full scan was performed again. The mean temperature was 20.403 C, with a measured range from 20.365 to 20.462 C.

| SiPM | Repeat Vbr | Shift from first clean run |
|---|---:|---:|
| Triangle | 50.486 +/- 0.121 V | -0.048 +/- 0.166 V |
| Star | 50.488 +/- 0.092 V | -0.179 +/- 0.134 V |

Neither shift was significant. The weighted two-run results were 50.512 +/- 0.083 V for Triangle and 50.574 +/- 0.067 V for Star.

## Interpretation

The corrected scan shows that the earlier unstable Triangle result was caused by the measurement path and trigger condition, not demonstrated SiPM physics. The CH3 artifact is still present below the chosen acquisition threshold and should be investigated electrically.

The two new SiPMs can provisionally be set near 53.51 and 53.57 V for 3.0 V overvoltage at about 20.4 C. These remain predicted effective biases until the high-voltage path is measured independently.

## Final state

High voltage was turned off after the scans. Temperature compensation remained stopped.

## Raw records

- `/home/muon/brDownVstudy/experiments/2026-09-23_new_sipms_triangle_star/`
