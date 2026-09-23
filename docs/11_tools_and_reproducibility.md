# Tools and reproducibility

## Detector hardware

| Item | Use in this study |
|---|---|
| Hamamatsu S13360-2050VE MPPC | SiPM under calibration |
| gLOWCOST four-channel readout PCB | Bias routing and non-inverting analogue amplification |
| Raspberry Pi detector controller | Programs MAX1932 and per-channel DAC values; records temperature |
| MAX1932 high-voltage supply | Common high-side SiPM bias |
| Per-channel DAC trim | Sets each channel's effective SiPM bias relative to the high side |
| PicoScope 3406B | Multi-channel waveform digitizer |
| 10:1 passive probes | Scope connections to amplified analogue outputs |
| BME280 sensor | Local temperature record |
| Light-tight SiPM connection / scintillator tiles | Dark-pulse and cosmic-scintillation tests |

The readout amplifier gain and bandwidth shape the measured height and area. The p.e. gap values in this repository describe the complete SiPM plus readout plus scope chain. They are not the bare SiPM charge unless that chain is calibrated.

## Computers and control

- Ubuntu waveform PC: `/home/muon/brDownVstudy`, used for PicoScope/PicoSDK acquisition and primary storage.
- Raspberry Pi: controlled over SSH during automated bias scans.
- Analysis and repository copy: macOS workstation.

The current hardware addresses belong in local configuration, not in published results. Passwords and private credentials are not stored in this repository.

## Software used on the acquisition PC

The recorded environment was:

| Software | Version |
|---|---:|
| Ubuntu/Linux kernel | 6.8 series |
| Python | 3.10.12 |
| NumPy | 2.2.6 |
| pandas | 2.3.3 |
| SciPy | 1.15.3 |
| Matplotlib | 3.10.8 |

PicoScope software was used for the early manual recordings. PicoSDK was used for automated acquisition. The Python analysis uses:

- `numpy` for arrays, integration, and linear algebra;
- `pandas` for run tables and event measurements;
- `scipy.ndimage.gaussian_filter1d` for diagnostic histogram smoothing;
- `scipy.signal.find_peaks` for historical candidate discovery;
- `scipy.optimize.curve_fit` for historical local Gaussian peak fits;
- `scipy.optimize.least_squares` for simultaneous spectrum models;
- `scipy.optimize.brentq` for the added-scatter solution;
- `matplotlib` for all scientific plots.

## Code layout

| Directory | Purpose |
|---|---|
| `scripts/acquisition/` | PicoSDK capture, Raspberry Pi bias commands, automated scans |
| `scripts/analysis/` | Current waveform processing, spectrum fit, Vbr fit, checks, and reporting |
| `scripts/legacy/` | Historical height/area peak-finding scripts used before the final comb method |
| `config/` | Reproducible scan settings and channel maps |
| `data/processed/` | Compact tables and selected per-run analysis outputs |
| `figures/` | Publication and historical diagnostic figures |
| `lab-log/` | Chronological observations and decisions |

The legacy scripts contain some original absolute paths. They are evidence of the analysis history and may require path edits before rerunning. The current scripts take command-line paths or configuration files.

## Reproducing one point

1. Verify the physical SiPM-to-PCB-to-scope mapping.
2. Record temperature and actual measured high-side voltage.
3. Load the scan configuration and confirm the 10:1 probe correction.
4. Set the bias, wait for settling, then acquire the requested events.
5. Run `scripts/analysis/analyze_point.py` on that voltage directory.
6. Inspect acceptance, baseline noise, waveform overlay, saturation, and both spectra.
7. Accept a p.e. gap only if the sequence/model tests pass.

The analysis entry points have `--help` output. A typical local form is:

```bash
python3 scripts/analysis/analyze_point.py --help
```

## Reproducing a Vbr result

1. Analyze every voltage independently.
2. Freeze the accepted-point list before looking at the intercept.
3. Fit area spacing versus temperature-normalized effective bias.
4. Repeat with peak height as a cross-check.
5. Run integration-gate and leave-one-out tests.
6. Compare an independent repeated scan.
7. Quote internal and absolute-systematic uncertainties separately.

## Raw-data policy

The waveform archive is much larger than a normal Git repository. Raw files remain on the acquisition PC under `/home/muon/brDownVstudy/experiments/`. `data/campaign_catalog.csv` gives the campaign path, device mapping, event scale, and analysis status. Selected compact outputs are committed so the plots and decisions can be audited without copying tens of gigabytes.

## Minimum metadata for future scans

Every new run should record:

- date, time, operator, and purpose;
- physical SiPM label and serial number if available;
- PCB and scope channel mapping;
- MAX1932 code, DAC code, commanded bias, and measured bias;
- temperature and the sensor location;
- probe ratio, scope range, sample interval, time window, coupling, and trigger;
- event target and accepted-event count;
- integration/search gates and analysis version;
- whether the run is accepted, diagnostic, or rejected, with the reason.

Without this information a beautiful spectrum may still be unusable.
