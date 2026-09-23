# 55.59 V waveform analysis

Raw waveform CSV files were kept on the Ubuntu PC. This folder contains extracted measurements, plots, and summary files only.

## sipm1 on ch0
- Run: `20260827_55p59V_sipm1_ch0_10000`
- Files: 10000; accepted after basic cuts: 9949
- Active scope channel: Channel A (median A peak 172.01 mV, median B peak 7.10 mV)
- Selected p.e.-like peak positions: 161.882;215.992;270.782 mV
- Estimated 1 p.e. pulse-height spacing: 54.45 +/- 0.55 mV
- Median pulse height: 172.30 mV
- Median pulse area: 8552.8 mV ns
- Median baseline RMS: 1.98 mV

## sipm2 on ch2
- Run: `20260827_55p59V_sipm2_ch2_10000`
- Files: 10000; accepted after basic cuts: 9871
- Active scope channel: Channel B (median A peak 5.82 mV, median B peak 171.16 mV)
- Selected p.e.-like peak positions: 170.025;212.112;253.921 mV
- Estimated 1 p.e. pulse-height spacing: 41.95 +/- 0.38 mV
- Median pulse height: 171.36 mV
- Median pulse area: 8515.7 mV ns
- Median baseline RMS: 2.39 mV

## Comparison
- 1 p.e. spacing ratio sipm2/sipm1: 0.770
- At the same nominal 55.59 V, this ratio suggests the second SiPM/readout channel has a smaller voltage step per p.e. than the first one.
- The most likely physics/electronics interpretation is different effective overvoltage, different SiPM breakdown voltage, or channel gain differences. A multi-bias scan is needed to separate those possibilities.

## Interpretation note
The peak spacing, not the absolute position of the first visible peak, is the useful 1 p.e. calibration. The absolute peak position depends on trigger threshold and on whether the acquisition preferentially recorded multi-p.e. pulses.