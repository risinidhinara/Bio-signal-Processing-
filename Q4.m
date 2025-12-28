fs = 500;   % Sampling frequency

%% (i) Low-pass Butterworth (order = 46, cutoff = 80 Hz)
M_LPF = 46; Wn_LPF = 80/(fs/2);  
[b_LPF, a_LPF] = butter(M_LPF, Wn_LPF, 'low');

% Visualize LPF characteristics
fvtool(b_LPF, a_LPF);
set(fvtool(b_LPF,a_LPF),'MagnitudeDisplay','Magnitude');
% Show phase response
set(fvtool(b_LPF,a_LPF), 'Analysis', 'phase');    

set(fvtool(b_LPF,a_LPF), 'Analysis', 'grpdelay');

% Butterworth highpass filter
M_HPF = 6; Wn_HPF = 0.5/(fs/2);
[b_HPF, a_HPF] = butter(M_HPF, Wn_HPF, 'high');

% Visualize LPF characteristics
fvtool(b_HPF, a_HPF);
set(fvtool(b_HPF,a_HPF),'MagnitudeDisplay','Magnitude');
set(fvtool(b_HPF,a_HPF), 'Analysis', 'phase'); 
set(fvtool(b_HPF,a_HPF), 'Analysis', 'grpdelay');

% Comb Filter 
fc = 50;  
q = 35;
bw = (fc / (fs/2)) / q;
[b_comb, a_comb] = iircomb(fs/fc, bw, 'notch');
fvtool(b_comb, a_comb);
set(fvtool(b_comb, a_comb),'MagnitudeDisplay','Magnitude');
set(fvtool(b_comb, a_comb), 'Analysis', 'phase'); 
set(fvtool(b_comb, a_comb), 'Analysis', 'grpdelay');

% Combined filter
LPF = dfilt.df1(b_LPF, a_LPF);
HPF = dfilt.df1(b_HPF, a_HPF);
Comb = dfilt.df1(b_comb, a_comb);
combinedFilter = dfilt.cascade(LPF, HPF, Comb);

% Plot the frequency response of the cascaded IIR filter
freqz(combinedFilter, 4096);

% 4.2 Filtering methods using IIR filters

% (i) Apply forward filtering

fwdFilter = filter(combinedFilter, nECG);

figure();
plot(t, nECG, t, fwdFilter);
title('Filtered nECG signal with forward filtering');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
legend('Noisy ECG', 'Forward Filtered ECG');
grid on;
figure('Position',[100 100 1200 600]);
plot(t, nECG, t, fwdFilter);
title('Filtered nECG signal with forward filtering');
xlabel('Time (s)');
xlim([0 2])
ylabel('Amplitude (mV)');
legend('Noisy ECG', 'Forward Filtered ECG');
grid on;
% (ii) Apply forward backward filtering

fwdbackFilter_1 = filtfilt(b_LPF, a_LPF, nECG);
fwdbackFilter_2 = filtfilt(b_HPF, a_HPF, fwdbackFilter_1);
fwdbackFilter = filtfilt(b_comb, a_comb, fwdbackFilter_2);

figure();
plot(t, nECG, t, fwdbackFilter);
title('Filtered nECG signal with forward-backward filtering');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
legend('Noisy ECG', 'Forward-Backward Filtered ECG');
grid on;
figure('Position',[100 100 1200 600]);
plot(t, nECG, t, fwdbackFilter);
title('Filtered nECG signal with forward-backward filtering');
xlim([0 2])
xlabel('Time (s)');
ylabel('Amplitude (mV)');
legend('Noisy ECG', 'Forward-Backward Filtered ECG');
grid on;

% (iii) Generate overlapping time domain plots

figure('Position',[100 100 1200 600]);
plot(t,nECG, t,ECG_CombF_dComp, t, fwdFilter, t, fwdbackFilter);
title('Overlapping time domain plots of filtered ECG');
xlim([0 2])
xlabel('Time (s)');
ylabel('Amplitude (mV)');
legend('nECG','FIR Filtered ECG', 'IIR Forward Filtered ECG', 'IIR Forward-Backward Filtered ECG');
grid on

% (iv) Generate overlapping plots of the PSDs

[p0,f0] = periodogram(nECG, window, N, fs);  % PSD estimate of nECG
[P1,f1] = periodogram(ECG_CombF_dComp, window, N, fs);  % PSD estimate of FIR filtered ECG
[P2,f2] = periodogram(fwdFilter, window, N, fs);   % PSD estimate of IIR forward filtered ECG
[P3,f3] = periodogram(fwdbackFilter, window, N, fs);   % PSD estimate of IIR forward-backward filtered ECG

figure('Position',[100 100 1200 600]);
plot(f0,10*log10(p0), f1, 10*log10(P1), f2, 10*log10(P2), f3, 10*log10(P3));   % Plotting the PSDs
title('Overlapping plots of PSDs of filtered ECG');
xlabel('Frequency (Hz)');
ylabel('Power/frequency (dB/Hz)');
legend('nECG' ,'FIR Filtered ECG', 'IIR Forward Filtered ECG', 'IIR Forward-Backward Filtered ECG');
grid on;
