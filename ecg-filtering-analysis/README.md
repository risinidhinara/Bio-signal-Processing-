ECG Filtering & Analysis — Summary and run instructions

Contents
- ecg_template_filtering.m — ECG template plotting and filter comparisons (MA, MA(N), SG)
- ensemble_averaging_and_segmentation.m — ABR ensemble averaging; ECG segmentation and SNR improvement
- fir_window_analysis_and_kaiser.m — FIR/window comparisons, Kaiser design, comb filter examples
- iir_cascade_forward_backward.m — IIR Butterworth + comb cascade; forward vs forward-backward filtering
- assignment1.mlx (MATLAB Live Script)

Overview
This folder contains the code and live script for "Assignment 1" (signal-processing exercises). The live script `assignment1.mlx` includes sections on:
- Preliminaries and plotting (ECG template)
- Moving-average and Savitzky–Golay filtering
- Ensemble averaging (ABR and ECG pulse segmentation)
- FIR/IIR filter design and window comparisons
- PSD and filter performance comparisons

Required data files (referenced by the live script)
- ABR_rec.mat
- ECG_rec.mat
- ECG_template.mat
- ECG_with_noise.mat

Requirements
- MATLAB R2018b or later (to open `.mlx` live scripts)
- Signal Processing Toolbox
- (Optional) Communications Toolbox for `awgn`

How to run

- Open `assignment1.mlx` in MATLAB Live Editor and run sections interactively.
- Or run the individual scripts in order from MATLAB command window:
  - `ecg_template_filtering.m`, `ensemble_averaging_and_segmentation.m`, `fir_window_analysis_and_kaiser.m`, `iir_cascade_forward_backward.m`

Notes
- The `.mlx` file contains embedded outputs and large `output.xml` data (this is expected for live scripts).
- Running the scripts will generate several figures. Ensure required `.mat` data files are on the MATLAB path or in the same folder before running.

Location
- This README and the files above are in the folder: /workspaces/Bio-signal-Processing-/ecg-filtering-analysis
