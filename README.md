# Bio-signal-Processing
ECG & Bio-signal Processing — repository overview

This repository contains MATLAB code and live scripts used for ECG and bio-signal processing exercises (filter design, ensemble averaging, Wiener/adaptive filtering, and wavelet time–frequency analysis).

Top-level folders

Quick start
1. Install MATLAB R2018b or later and add the Signal Processing Toolbox. For wavelet examples also add the Wavelet Toolbox. The `awgn` function requires Communications Toolbox or you can replace it with `randn`-based noise.
2. Put any required `.mat` data files (listed in each folder README) into the relevant folder or add the folder to the MATLAB path.
3. Open the relevant `.mlx` live script in MATLAB Live Editor to run step-by-step, or run the plain `.m` scripts in the order suggested in each folder README.

Data files referenced across the repo

Repository maintenance notes

Contact


ECG & Bio-signal Processing — repository overview

This repository contains MATLAB code and live scripts for ECG and bio-signal processing exercises, including filter design, ensemble averaging, Wiener/adaptive filtering, and wavelet time–frequency analysis.

Top-level folders
- ecg-filtering-analysis: ECG template plotting, moving-average and Savitzky–Golay filtering, FIR window/Kaiser design, comb/IIR filtering, ensemble averaging. See ecg-filtering-analysis/README.md
- wiener-adaptive-filtering-ecg: Wiener filter (time and frequency domain), adaptive filtering (LMS, RLS) and ECG experiments. See wiener-adaptive-filtering-ecg/README.md
- wavelet-time-frequency-analysis: Mexican-hat CWT examples, DWT decomposition/reconstruction/denoising (Haar, DB9). See wavelet-time-frequency-analysis/README.md

Quick start
1. Install MATLAB R2018b or later and add the Signal Processing Toolbox. For wavelet examples also add the Wavelet Toolbox. The `awgn` function requires the Communications Toolbox or can be replaced with `randn`-based noise.
2. Place any required `.mat` data files (listed in each folder README) into the relevant folder or add the folder to the MATLAB path.
3. Open the relevant `.mlx` live script in MATLAB Live Editor to run interactively, or run the plain `.m` scripts in the order suggested in each folder README.

Data files referenced across the repo
- ECG_template.mat, ECG_rec.mat, ECG_with_noise.mat, ABR_rec.mat (used by ecg-filtering-analysis)
- idealECG.mat (used by wiener-adaptive-filtering-ecg)

Repository notes
- Filenames and folder names have been updated to be descriptive and consistent.

If you need commits, a release zip, or `.mlx`→`.m` conversions, I can prepare those on request.
