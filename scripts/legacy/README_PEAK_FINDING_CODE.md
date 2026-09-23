# p.e. peak finding and breakdown-voltage analysis code

This folder is a snapshot of the Python scripts used for the SiPM waveform
p.e. peak finding, p.e. gap extraction, histogram plotting, and breakdown
voltage fits.

The main method is:

1. Read PicoScope CSV waveforms.
2. Convert scope units to mV.
3. Estimate the baseline with robust median/MAD logic.
4. Subtract the baseline from each waveform.
5. Measure prompt peak height and prompt pulse area.
6. Build a histogram of peak heights or pulse areas.
7. Smooth the histogram with `scipy.ndimage.gaussian_filter1d`.
8. Locate candidate peaks with `scipy.signal.find_peaks`.
9. Refine each peak position with local Gaussian-plus-linear-background fits
   using `scipy.optimize.curve_fit`.
10. Select the best nearly-even sequence of peak centers.
11. Calculate the 1 p.e. spacing from adjacent peak-center differences.
12. Plot the histogram with selected p.e. peaks marked.
13. Plot 1 p.e. spacing versus bias voltage and extrapolate the x intercept
    to estimate breakdown voltage.

Important scripts:

- `robust_sipm1_pe_analysis.py`: robust waveform processing and peak-height
  p.e. spacing for SIPM1.
- `replot_sipm2_final_peak_plots.py`: SIPM2 peak-height spectra and final
  p.e. spacing plots.
- `build_clean_height_gap_linear_plot.py`: combined final peak-height
  histogram panels and breakdown-voltage line fits.
- `build_area_gap_linear_plot.py`: pulse-area version of the same p.e. gap
  analysis.
- `analyze_waveform_run.py`: generic single-run waveform analysis script.
- `build_three_panel_waveform_peak_plot.py`: one-run panel with histogram,
  waveform overlay, and height-versus-area behavior.
- `plot_checked_vbr_line_fit.py`: final checked Vbr fit plotting.
- `plot_all_waveform_histogram_runs.py`: bookkeeping panels for all waveform
  histogram runs.
- `plot_raw_pedestal_hist.py`: histogram before baseline subtraction.
- `rediscover_height_peaks_audit.py`: diagnostic script for rechecking peak
  discovery and selected p.e. gaps.

The scripts were written to run on the Ubuntu waveform PC with the data under:

`/home/muon/brDownVstudy`

Some plotting scripts still contain absolute output paths from the original
analysis folders. If a run folder name changes, update the `RUNS`, `ROOT`, or
manifest path at the top of the relevant script before running.

Minimum Python packages:

- numpy
- pandas
- matplotlib
- scipy
