# References

## Device information

1. Hamamatsu Photonics, [S13360-2050VE product page](https://www.hamamatsu.com/eu/en/product/optical-sensors/mppc/mppc_mppc-array/S13360-2050VE.html).
2. Hamamatsu Photonics, [MPPC S13360 series datasheet](https://www.hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/s13360_series_kapd1052e.pdf).
3. Hamamatsu Photonics, [MPPC technical information](https://www.hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/ssd/mppc_kapd9005e.pdf).

## SiPM characterization and breakdown voltage

4. R. Klanner, [Characterisation of SiPMs](https://doi.org/10.1016/j.nima.2018.11.083), *Nuclear Instruments and Methods in Physics Research A* 926 (2019) 36-56.
5. V. Chmill et al., [On the characterisation of SiPMs from pulse-height spectra](https://doi.org/10.1016/j.nima.2017.02.049), *Nuclear Instruments and Methods in Physics Research A* 854 (2017) 70-81.
6. V. Chmill et al., [Study of the breakdown voltage of SiPMs](https://arxiv.org/abs/1605.01692), arXiv:1605.01692.

## How these references are used

The Hamamatsu documents define the device-family operating limits and typical values. They do not replace measurement of an individual SiPM. The Klanner review gives the wider characterization framework. The Chmill papers support extracting gain from multi-p.e. spectra and determining breakdown voltage by extrapolating the gain-related peak spacing.

The analysis choices in this repository are also constrained by the measured readout chain. In particular, the acquisition-trigger study and the lower-bias exclusion rules come from this apparatus rather than being copied from a general SiPM prescription.
