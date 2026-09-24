# EPIC-tile trigger calibration: 2026-09-23

## Physical arrangement

The two GSU gLOWCOST tiles were placed between the two EPIC reference tiles. The channel labels for this experiment are:

| Position | SiPM label | PCB channel | PicoScope channel |
|---|---|---:|---|
| Top EPIC tile | not individually calibrated | CH0 | C |
| Top GSU gLOWCOST tile | SiPM 1 (Triangle) | CH2 | A |
| Bottom GSU gLOWCOST tile | SiPM 2 (Star) | CH3 | B |
| Bottom EPIC tile | not individually calibrated | CH1 | D |

The purpose is to use the EPIC tiles to select a particle crossing the stack, then measure the probability and signal distribution seen by the two GSU gLOWCOST tiles.

## Bias preparation

The first attempt set Triangle and Star to their measured `Vbr + 3 V` values, but CH0 and CH1 were not yet assigned controlled biases. That run was rejected.

A second attempt used the Hamamatsu typical absolute Vbr of 53 V for the EPIC-tile SiPMs on CH0 and CH1. At 20.82 C this placed both channels near 55.77 V. This was also rejected as the operating-point choice: CH0 saturated the scope and its response was not on the same board-referenced voltage scale as the measured Triangle and Star intercepts.

The corrected provisional setting uses the mean measured Vbr of Triangle and Star, 50.543 V at 20.4 C, as a proxy for each EPIC-tile SiPM. This is not presented as a measurement of their Vbr. At 20.94 C the applied effective biases were:

| PCB channel | Effective bias |
|---:|---:|
| CH0, top EPIC tile | 53.572 V |
| CH1, bottom EPIC tile | 53.572 V |
| CH2, SiPM 1 (Triangle) | 53.540 V |
| CH3, SiPM 2 (Star) | 53.602 V |

The shared MAX1932 setting was `0xF8`. The low-side DAC settings were `0x215`, `0x215`, `0x220`, and `0x20B` for CH0 through CH3. Temperature compensation remained stopped during the run.

## Trigger development

The first C-D analogue-AND pilot used 700 mV on both reference channels. All 100 records arrived at the 5 s auto-trigger interval, so this condition was too strict.

Single-channel tests showed that D pulses were much smaller than C pulses. A D-triggered pilot at the overbiased EPIC-tile setting produced real, simultaneous pulses, while an HV-off control produced only timeout records. This confirmed that the live waveforms were detector signals rather than a persistent scope artifact.

After correcting the EPIC-tile biases, a 50 mV trigger on D collected mostly low-amplitude dark pulses. A 5,000-event run was used to define a cleaner reference region. The useful starting cuts were:

- C positive excursion at least 3,000 mV;
- D positive excursion at least 120 mV;
- absolute C-D peak-time difference no more than 40 ns.

Six events passed these strict cuts in 81.49 s, corresponding to 4.4 events/min. Four of the six had Triangle and Star excursions above 200 mV. This sample is too small for an efficiency result, but its rate is compatible with the few-per-minute FPGA CH0-CH1 coincidence scale.

## Production acquisition

The production run was started at 18:11 EDT:

```text
/home/muon/brDownVstudy/experiments/2026-09-23_triangle_star_scintillator_calibration/
```

The run directory is identified by its `20260923_181800` start-time suffix on the acquisition computer.

Settings:

| Item | Value |
|---|---:|
| Requested events | 1,000 |
| Hardware trigger | Scope D, 150 mV at probe tip |
| Auto-trigger | disabled |
| Scope input ranges | A 200 mV, B 200 mV, C 2,000 mV, D 100 mV |
| Probe attenuation | 10:1 |
| Sampling interval | 4 ns actual |
| Time window | 800 ns before, 400 ns after trigger |

The A and B ranges were chosen to retain sensitivity to small signals. Saturation of large muon pulses is acceptable for the first pass because the main quantity is hit or miss probability, not the full large-pulse amplitude.

An initial production attempt used a 500 mV scope-input range on C. It was stopped after about 20 events because earlier pilots showed that C could overflow this range. The restarted run uses 2,000 mV on C. The analysis was also corrected to reject an event only when a **reference** channel overflows; saturation on A or B is retained as an unambiguous detector hit.

The observed starting rate was about 3 events/min, giving an estimated duration of 5 to 6 hours. The acquisition configuration and the live bias-state JSON were copied into the run directory.

## Status and limitations

This is a full-chain scintillator response measurement, not an intrinsic PDE measurement. The two EPIC-tile SiPM breakdown voltages are still unknown, and their provisional bias must eventually be replaced by individual measurements or a demonstrated reference-efficiency plateau. The production result must not be quoted until the reference cuts, accidental contribution, waveform saturation, and GSU gLOWCOST-tile zero-event threshold have been checked.
