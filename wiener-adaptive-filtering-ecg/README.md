Wiener & Adaptive Filtering (ECG) — Summary and run instructions

Contents
- wiener_filter_ecg_analysis.m — Wiener filter (time- and frequency-domain) analysis for ECG; linear model and non-stationary noise experiments
- adaptive_filtering_lms_rls.m — Adaptive filtering: LMS and RLS implementations, parameter sweeps, and ECG tests
- wiener_adaptive_analysis.mlx (MATLAB Live Script)

Overview
This folder contains scripts and a live script for Wiener and adaptive filtering experiments applied to ECG and synthetic signals. Key topics:
- Time-domain Wiener filter design and optimal FIR Wiener filter computation
- Frequency-domain Wiener filter implementation
- Non-stationary noise experiments (50 Hz/100 Hz transitions)
- Adaptive filtering using LMS and RLS; parameter sweeps and comparisons

Required data files
- idealECG.mat
- (scripts may create) r_n.mat and x_in.mat — intermediate saved variables used by adaptive scripts

Requirements
- MATLAB R2018b or later
- Signal Processing Toolbox

How to run
- Open `wiener_adaptive_analysis.mlx` in MATLAB Live Editor and run sections interactively.
- Or run scripts from MATLAB command window in this order:
  - `wiener_filter_ecg_analysis.m`
  - `adaptive_filtering_lms_rls.m`

Notes
- `adaptive_filtering_lms_rls.m` saves `r_n.mat` and `x_in.mat` which are then used by the live script; ensure you run it before running any sections of the live script that expect those files.
- Figures will be produced by the scripts; ensure required `.mat` files are in the same folder or on MATLAB path.

Location
- This README and the files above are in the folder: /workspaces/Bio-signal-Processing-/wiener-adaptive-filtering-ecg
