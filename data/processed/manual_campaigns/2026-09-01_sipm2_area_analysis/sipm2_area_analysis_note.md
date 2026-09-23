# SIPM2 area-based p.e. spacing

This repeats the SIPM1 area method for SIPM2 using detector channel CH2 recorded on PicoScope Channel B.

Available SIPM2 waveform runs were found at 54.10, 54.90, 55.59, 55.65, 56.30, 56.40, and 56.86 V. No 53.20 V or 53.90 V SIPM2 waveform folder was found under `/home/muon`.

The area is the prompt baseline-subtracted waveform integral:

```text
area = integral[V(t) - baseline] dt
```

The unit is `mV ns`.

All available SIPM2 area-gap points were used in the line fit:

```text
area spacing slope = 646.308 mV ns/V
area-fit Vbr = 52.035 +/- 0.342 V
R^2 = 0.9603
fit RMS = 115.6 mV ns
n = 7
```
