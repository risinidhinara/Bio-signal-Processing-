Wavelet Time–Frequency Analysis — Summary and run instructions

Contents
- wavelet_time_frequency_analysis.mlx — MATLAB Live Script combining wavelet examples and exercises
- wavelet_mexican_hat_cwt.m — Mexican hat mother/daughter wavelets, CWT via convolution, spectra
- wavelet_dwt_reconstruction_denoising.m — DWT decomposition (Haar, DB9), reconstruction, denoising, and energy checks

Overview
This folder contains scripts and a live script focused on wavelet-based time–frequency analysis and wavelet denoising. Key topics:
- Mexican hat mother wavelet and its scaled (daughter) forms
- Continuous Wavelet Transform (CWT) via convolution and scale-time spectrograms
- Discrete Wavelet Transform (DWT) using Haar and Daubechies-9: decomposition, reconstruction, and denoising
- PSD and energy verification for reconstructed signals

Data/Dependencies
- No external .mat data required for these scripts (they generate synthetic signals).
- Required MATLAB toolboxes: Wavelet Toolbox, Signal Processing Toolbox

How to run
- Open `wavelet_time_frequency_analysis.mlx` in MATLAB Live Editor and run interactively.
- Or run the scripts from MATLAB command window in this order:
  - `wavelet_mexican_hat_cwt.m`
  - `wavelet_dwt_reconstruction_denoising.m`

Notes
- Scripts generate multiple figures; run in MATLAB with sufficient display capability.
- If you want I can also extract the live script sections into plain `.m` files.

Location
- This README and the files above are in the folder: /workspaces/Bio-signal-Processing-/wavelet-time-frequency-analysis
