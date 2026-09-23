# Next measurements

## 1. Calibrate the delivered voltage

The first priority is to measure the high side and channel low side with a high-impedance calibrated meter at several commanded values. The effective SiPM bias is their difference. The calibration should be repeated under the same load and wiring used in the waveform scan.

Required output:

- measured voltage versus MAX1932 code;
- measured channel low-side voltage versus DAC code;
- residuals from the calibration model;
- uncertainty of the effective-bias difference;
- check for hysteresis by approaching points from both directions.

## 2. Measure Vbr as a function of temperature

Place the detector and temperature sensor in a controlled enclosure. At each stable temperature, allow thermal equilibrium and repeat a shortened p.e.-spacing scan. Fit

```text
Vbr(T) = Vbr(Tref) + kT (T - Tref)
```

for each physical SiPM. The fitted `kT` should replace the provisional 54 mV/C value in the compensation program.

Temperature steps should be held long enough for the SiPM, board, and sensor to agree. A stable air temperature alone is not proof that the SiPM junction has equilibrated.

## 3. Pulsed-LED zero-event measurement

A pulsed blue LED near 470 nm is suitable for testing the wavelength-shifting-fibre readout region. The LED should be driven by a stable pulser and split or diffused so that both SiPMs receive the same weak light field.

For each light pulse, integrate the waveform in a fixed signal gate. Measure the pedestal fraction with the LED on and with random/dark gates. The mean detected primary avalanches can be estimated from the corrected zero fraction:

```text
mu = -ln(P0)
```

The light level should remain low enough that the pedestal remains populated. Repeating the measurement at several overvoltages gives a relative PDE curve. The operating points can then be adjusted until the two devices have the required matched response.

## 4. Scintillator efficiency check

Use the small reference tiles to select through-going particles. Record enough coincidences to compare the standard tiles at their provisional `Vbr + 3 V` settings. The result should be reported as full-chain detection efficiency, not intrinsic PDE.

The following should be varied deliberately:

- discriminator or software pulse threshold;
- SiPM overvoltage;
- trigger geometry;
- tile/fibre position;
- temperature.

## 5. Investigate PCB CH3

The prompt artifact should be measured directly at several points along the CH3 analogue path. Compare CH3 with a clean channel using the same SiPM and probe:

- connector input;
- amplifier input and output;
- local power rails;
- ground reference;
- neighboring digital activity.

The 25 mV acquisition trigger is a useful workaround for the present calibration. It is not the final explanation.

## 6. Choose the detector operating point

The final operating point should be chosen from several measurements together:

1. stable p.e. gain and resolved calibration peaks;
2. acceptable dark-count, crosstalk, and afterpulse behavior;
3. high and stable scintillator detection efficiency;
4. matched response between channels;
5. sufficient margin below amplifier and digitizer saturation;
6. stable performance over the required temperature range.

The target does not have to be exactly 3.0 V overvoltage. That value is a controlled starting point for comparison.
