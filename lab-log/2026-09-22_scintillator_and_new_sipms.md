# Scintillator pilot and new SiPMs: 2026-09-22

## Scintillator zero-event pilot

Two small reference tiles were placed above and below the standard detector tiles. Their coincidence selected particles crossing the central region. The two standard-tile waveforms were then inspected for signal and zero events.

This test was prepared as a full-chain efficiency measurement. It includes the scintillator, wavelength-shifting fibre, optical coupling, SiPM, amplifier, and analysis threshold. It should not be described as an intrinsic PDE measurement.

The pilot established the acquisition geometry and the corrected zero-fraction calculation. More reference coincidences are required before comparing operating points with useful precision.

## New SiPM pair

Two light-blocked SiPMs were connected and physically marked Triangle and Star. The first comparisons included channel swaps intended to separate a device effect from a readout-path effect.

The Star measurement on the clean channel was reproducible. The Triangle result on PCB CH3 was not. Only about 2 to 4% of CH3 acquisitions survived the pulse selection, and the apparent Vbr moved between repeats. A high R2 from three remaining points was rejected because the distant intercept was poorly constrained.

The important observation was that poor acceptance followed the PCB CH3 measurement path rather than the physical SiPM. Median CH3 acquisitions contained an approximately 8 to 9 mV prompt feature at trigger time, while only a small subset contained the later SiPM-like pulse family.

## Decision

Do not average the invalid CH3 estimates. First determine a scope trigger that suppresses the prompt artifact while retaining the resolved p.e. spacing, then repeat the complete scan for both physical devices with matched settings.

## Raw records

- `/home/muon/brDownVstudy/experiments/2026-09-22_scintillator_zero_event_test/`
- The first Triangle/Star runs are retained inside the 2026-09-23 campaign after directory reorganization.
