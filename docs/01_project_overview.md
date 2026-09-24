# Project overview

## Why this study was started

The gLOWCOST detector uses plastic scintillators, wavelength-shifting fibres, and SiPMs to detect charged particles. A charged particle crossing the scintillator produces light. Part of this light reaches the SiPM and produces a short current pulse. The readout amplifies the pulse and compares it with a threshold before the FPGA processes it. The processed counts are sent to the Raspberry Pi that hosts the detector readout, where the data are prepared and transferred to the Georgia State University servers.

The SiPM must be biased above its breakdown voltage. The difference

```text
Vover = Vbias - Vbr
```

is called the overvoltage. Increasing overvoltage increases the avalanche gain and usually raises the photon-detection efficiency, but it also raises dark counts, optical crosstalk, and afterpulsing. The best operating point is therefore not simply the largest safe voltage.

The SiPMs also do not have exactly the same breakdown voltage. Applying one common bias can place two devices at different overvoltages. The first goal of this study was to measure the breakdown voltage of each SiPM using the existing detector electronics. The larger goal is to operate several SiPMs at comparable sensitivity and then verify that condition with light or particle data. A pulsed-LED zero-event measurement is the planned controlled calibration after the particle-trigger method has been validated.

## Questions being tested

1. Can the detector readout and PicoScope resolve photoelectron populations from dark pulses?
2. Does the spacing between those populations change linearly with effective bias?
3. Is the extrapolated breakdown voltage reproducible when a scan is repeated?
4. Does the result follow the physical SiPM when channels are exchanged?
5. Which trigger settings record the pulses without selecting a misleading subset?
6. After the breakdown voltages are known, how should equal response be tested?

## How the work developed

The first waveform recordings were made manually. They showed visible groups of pulse heights, but the file naming, trigger settings, and channel changes were not controlled well enough for a final comparison. The acquisition was then automated so that each voltage point used the same sampling, event count, bias settling time, and metadata record.

The first SiPM pair produced repeatable p.e.-area spacing curves. A lower-bias scan was then attempted to move closer to breakdown. This was useful, but it also showed an important limit: once neighboring p.e. populations overlap and the trigger accepts only the larger pulses, a fitted spacing can be biased even when a fit converges.

A scintillator coincidence setup was prepared next. This is a different measurement. It can measure the detection efficiency of the complete scintillator-fibre-SiPM chain, but it does not directly give the intrinsic SiPM photon-detection efficiency.

Two new bare devices, identified as SiPM 1 (△) and SiPM 2 (★), were later connected to PCB channels 3 and 2. Early scans gave inconsistent behavior on channel 3. Swapping scope inputs showed that the problem stayed with the PCB path. A trigger scan found a sharp acceptance change between 15 and 25 mV. The full bias scan was repeated at 25 mV, followed by another complete repeat as a sanity check.

## Main conclusions at this stage

- Pulse-area spacing gives a repeatable internal measurement of breakdown voltage.
- The SiPM 1 (△) and SiPM 2 (★) devices have compatible breakdown voltages within the present uncertainty.
- PCB CH3 has a low-amplitude prompt artifact. A 25 mV PicoScope trigger avoids it during this calibration, but the electronics path still deserves inspection.
- Very low-bias points are not automatically better points. Trigger efficiency and peak overlap make them less reliable in the present setup.
- The reported Vbr errors do not yet include an independently measured absolute bias-voltage systematic.
- A pulsed LED zero-event measurement is still required before claiming equal PDE. A scintillator coincidence can test the full detector efficiency while that light source is being prepared.

## Device groups must remain separate

Two different physical pairs appear in this repository:

| Group | Labels | Purpose |
|---|---|---|
| Original pair | SIPM1 and SIPM2 | Method development, repeated scans, lower-bias study, and initial operating-point run |
| New pair | SiPM 1 (△) and SiPM 2 (★) | Channel diagnosis and final two-run comparison on 2026-09-23 |

Their values must not be combined. The labels identify physical devices, not permanent PCB channels.
