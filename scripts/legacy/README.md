# Historical peak-finding scripts

These scripts are the analysis record from the manual waveform stage. They are kept because the early plots and provisional breakdown fits should be reproducible from the method actually used at the time.

The common historical path was:

1. read PicoScope CSV waveforms;
2. estimate and subtract the event baseline;
3. calculate prompt peak height and/or waveform area;
4. histogram the event values;
5. smooth the histogram with `scipy.ndimage.gaussian_filter1d`;
6. locate candidates with `scipy.signal.find_peaks`;
7. refine local centers with Gaussian fits from `scipy.optimize.curve_fit`;
8. select a nearly even sequence and calculate its mean gap;
9. fit gap versus bias and extrapolate to zero.

`README_PEAK_FINDING_CODE.md` describes the individual files. Several scripts retain absolute paths from `/home/muon/brDownVstudy` and are therefore snapshots rather than polished command-line tools.

The current primary method is in `../analysis/`. It uses peak-aligned pulse area, a pedestal-constrained multi-population fit, model-selection checks, repeated scans, and explicit systematic tests. Do not replace the current result with a legacy output merely because the historical plot looks cleaner.
