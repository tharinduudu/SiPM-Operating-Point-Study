# Scripts

## Acquisition

The [acquisition](acquisition) scripts control the PicoScope capture, communicate with the Raspberry Pi bias controller, and save run metadata. Review host names, SSH settings, channel maps, and hardware addresses before use.

## Analysis

The [analysis](analysis) scripts perform baseline subtraction, waveform measurements, p.e.-spectrum fitting, breakdown-voltage fitting, repeatability checks, and figure generation.

The scripts were copied from the working study directory at the time this repository was assembled. Processed outputs are retained so that later code changes do not silently rewrite the recorded result.

## Safety

Bias-control code can command detector high voltage. A scan must leave HV off on normal completion and on failure. Because the MAX1932 path has no readback, verify the physical output before connecting or moving SiPMs.
