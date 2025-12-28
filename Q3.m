%% Overlay Impulse Responses of FIR Filters for different M
wc = 0.4*pi;    % Cutoff frequency in radians
fc = wc/pi;     % Normalized cutoff

M_values = [5, 50, 100];   % Different filter lengths
colors = {'r','b','k'};    % Plot colors

figure('Position',[100 100 1200 600]); hold on;
for i = 1:length(M_values)
    M = M_values(i);
    b = fir1(M-1, fc, rectwin(M));   % FIR filter coefficients (impulse response)

    stem(0:length(b)-1, b, colors{i}, 'filled'); 
end

xlabel('Sample number n');
ylabel('h_n');
title('Truncated Impulse Response of FIR Lowpass Filter');
legend('M=5','M=50','M=100');
grid on;



%% (i) Effect of Rectangular Window Length M
%% FIR Lowpass filter design with Rectangular Window

wc = 0.4*pi;        % Cutoff frequency in radians
fc = wc/pi;         % Normalized cutoff (0 to 1)

M_values = [5, 50, 100];   % Different window lengths
colors = {'b','r',[1, 0.65, 0]};     % Colors for plotting, black 'k' changed to orange RGB

figure;
hold on;
grid on;
title('FIR Lowpass Filter Frequency Response with Rectangular Window');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude (absolute)');
legend_labels = {};

for i = 1:length(M_values)
    M = M_values(i);
    b = fir1(M-1, fc, rectwin(M));  % FIR filter design
    [H, w] = freqz(b, 1, 1024);     % Frequency response
    
    plot(w/pi, abs(H), 'Color', colors{i}); % Magnitude response
    legend_labels{i} = ['M = ', num2str(M)];
end

legend(legend_labels);

%% FIR Lowpass filter with Rectangular Window - Phase response only
wc = 0.4*pi;        % Cutoff frequency
fc = wc/pi;         % Normalized cutoff (0–1)

M_values = [5, 50, 100];   % Filter lengths
Hd = cell(length(M_values),1);

for i = 1:length(M_values)
    M = M_values(i);
    b = fir1(M-1, fc, rectwin(M));   % FIR filter coefficients
    Hd{i} = dfilt.dffir(b);          % Convert to filter object
end

% Plot only the phase response
fvtool(Hd{:}, 'Analysis', 'phase');
legend('M=5','M=50','M=100');

%% (iii) Window comparison: Rectangular, Hann, Hamming, Blackman (M=50)
% Length of the window
M = 50;

% Create windows
rect_win = rectwin(M);
hann_win = hann(M);
hamming_win = hamming(M);
blackman_win = blackman(M);

% Time-domain plot
figure;
plot(rect_win); hold on;
plot(hann_win);
plot(hamming_win);
plot(blackman_win);
title('Time domain');
xlabel('Samples');
ylabel('Amplitude');
legend('Rectangular', 'Hanning', 'Hamming', 'Blackman');
grid on;



wc = 0.4*pi;

h_rect  = fir1(M-1, wc/pi, rectwin(M));
h_hann  = fir1(M-1, wc/pi, hann(M));
h_hamm  = fir1(M-1, wc/pi, hamming(M));
h_black = fir1(M-1, wc/pi, blackman(M));

Hd_rect  = dfilt.dffir(h_rect);
Hd_hann  = dfilt.dffir(h_hann);
Hd_hamm  = dfilt.dffir(h_hamm);
Hd_black = dfilt.dffir(h_black);

% Linear magnitude
fvtool(Hd_rect,Hd_hann,Hd_hamm,Hd_black, ...
       'Analysis','magnitude', 'MagnitudeDisplay','Magnitude');
legend('Rectangular','Hann','Hamming','Blackman');
set(gcf,'Name','Linear Magnitude, M=50');
% Logarithmic magnitude (in dB)
fvtool(Hd_rect,Hd_hann,Hd_hamm,Hd_black, ...
       'Analysis','magnitude', 'MagnitudeDisplay','Magnitude (dB)');
legend('Rectangular','Hann','Hamming','Blackman');
set(gcf,'Name','Log Magnitude (dB), M=50');

% Impulse response (window morphology)
fvtool(Hd_rect,Hd_hann,Hd_hamm,Hd_black, 'Analysis','impulse');
legend('Rectangular','Hann','Hamming','Blackman');
set(gcf,'Name','Impulse Response (Morphology), M=50');

% Impulse response (window morphology)
fvtool(Hd_rect,Hd_hann,Hd_hamm,Hd_black, 'Analysis','phase');
legend('Rectangular','Hann','Hamming','Blackman');
set(gcf,'Name','Impulse Response (Morphology), M=50');



data = load('ECG_with_noise.mat');

N = length(nECG);   % Number of samples in the file
fs = 500;
t = (0 : N-1) / fs; % Calculating time points

% Plotting the loaded ECG signal
figure();
plot(t, nECG)
title('nECG Signal')
xlabel('Time (s)')
ylabel('Amplitude (mV)')
grid on

% Plot time-domain signal
figure;
plot(t, nECG);
xlabel('Time (s)'); ylabel('Amplitude');
title('Noisy ECG Signal');
xlim([0 2])

window = rectwin(N);    % Using Rectangular window

% Plotting the power spectral density (PSD) estimate of the nECG
figure();
periodogram(nECG, window, N, fs)

% (iii) Calculate the relevant β and M

[M_HPF, Wn_HPF, beta_HPF, ftype_HPF] = kaiserord([0.2 0.7], [0 1], [0.001 0.001], fs);
[M_LPF, Wn_LPF, beta_LPF, ftype_LPF] = kaiserord([60 100], [1 0], [0.001 0.001], fs);

fprintf('The optimum filter parameters for HPF are M = %d and β = %f', M_HPF, beta_HPF);
fprintf('The optimum filter parameters for LPF are M = %d and β = %f', M_LPF, beta_LPF);

% (iv) Visualise the windows (highpass and lowpass)

% Implement the filters with Kaiser window
a = 1;  % Denominator coefficient
% Numerator coefficients
b_HPF = fir1(M_HPF, Wn_HPF, ftype_HPF, kaiser(M_HPF+1, beta_HPF));
b_LPF = fir1(M_LPF, Wn_LPF, ftype_LPF, kaiser(M_LPF+1, beta_LPF));

% High Pass Filter Response
figure();
freqz(b_HPF, a);
title('High Pass Filter Magnitude Response');
% Low Pass Filter Response
figure();
freqz(b_LPF, a);
title('Low Pass Filter Magnitude Response');

% Apply HPF
ECG_HPF = filter(b_HPF, a, nECG);

groupDelay_H = ceil(M_HPF/2); % Grp delay = order / 2
% Compensating for group delay
ECG_HPF_dComp = [ECG_HPF((groupDelay_H+1) : length(ECG_HPF)), ECG_HPF(1: groupDelay_H)];

% Effect of HPF in time domain
figure();
plot(t, nECG, t, ECG_HPF_dComp);
title('High pass filtered nECG signal with delay compensation');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
legend('Noisy ECG', 'Delay Compensated ECG');
grid on;

% Define duration for plotting (0 to 2 seconds)
t_end = 2;                          % seconds
idx = t <= t_end;               

% Effect of HPF in time domain 
figure('Position',[100 100 1200 600]);
plot(t(idx), nECG(idx)); hold on;
plot(t(idx), ECG_HPF_dComp(idx));
title('High pass filtered nECG signal with delay compensation (0–3 s)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
legend('Noisy ECG', 'Delay Compensated ECG');
grid on;

% Apply LPF
ECG_LPF = filter(b_LPF, a, ECG_HPF_dComp);

groupDelay_L = ceil(M_LPF/2);
% Compensating for group delay
ECG_LPF_dComp = [ECG_LPF((groupDelay_L+1) : length(ECG_LPF)), ECG_LPF(1: groupDelay_L)];
% Effect of LPF in time domain
figure;
plot(t, nECG, t, ECG_LPF_dComp);
title('Low pass filtered nECG signal with delay compensation');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
legend('Noisy ECG', 'Delay Compensated ECG');
grid on;


% Effect of LPF in time domain
figure('Position',[100 100 1200 600]);
plot(t, nECG, t, ECG_LPF_dComp);
title('Low pass filtered nECG signal with delay compensation');
xlim([0 2])
xlabel('Time (s)');
ylabel('Amplitude (mV)');
legend('Noisy ECG', 'Delay Compensated ECG');
grid on;

% FIR comb filter design with 3 cut off frequencies
f_stop = [50, 100, 150];

% Calculate zeros on the unit circle for each frequency and their conjugates
zeros = [];
for f = f_stop
    zero = exp(-1i * 2 * pi * f / fs);
    zeros = [zeros, zero, conj(zero)];
end

% Calculate the polynomial coefficients from the zeros
coefficients = poly(zeros);
% Normalize the coefficients so that the gain at DC is 1
b_comb = coefficients / sum(coefficients);

% Comb Filter Response
freqz(b_comb, a, 4096, fs);
title('Comb Filter Magnitude Response');
% Apply Comb Filter
ECG_CombF = filter(b_comb, a, ECG_LPF_dComp);

groupDelay_C = ceil((length(b_comb)-1)/2);  % Grp delay = order / 2
% Compensating for group delay
ECG_CombF_dComp = [ECG_CombF((groupDelay_C+1) : length(ECG_CombF)), ECG_CombF(1: groupDelay_C)];

% Effect of Comb filtering in time domain
figure();
plot(t, nECG, t, ECG_CombF_dComp);
title('Comb filtered nECG signal with delay compensation');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
legend('Noisy ECG', 'Delay Compensated ECG');
grid on;


figure('Position',[100 100 1200 600]);
plot(t, nECG, t, ECG_CombF_dComp);
title('Comb filtered nECG signal with delay compensation');
xlabel('Time (s)');
xlim([0 2])
ylabel('Amplitude (mV)');
legend('Noisy ECG', 'Delay Compensated ECG');
grid on;
% (vi) Plot the magnitude response of the combined three filters and the PSD of the final filtered ECG

% Assuming b_HPF, b_LPF, and b_comb are already defined
% Assuming the 'a' parameter for FIR filters is 1.

% Frequency response of the filters using fvtool
fvtool(b_HPF, 1, b_LPF, 1, b_comb, 1);
title('Magnitude Response of Cascaded Highpass, Lowpass, and Comb Filters');
legend('Highpass Filter', 'Lowpass Filter', 'Comb Filter');

% To plot the combined response of the cascaded filter:
H_cascaded = conv(b_HPF, conv(b_LPF, b_comb));
figure;
fvtool(H_cascaded, 1);
title('Magnitude Response of the Combined Cascaded Filter');
% Effect of cascaded filter in time domain
figure();
plot(t, nECG, t, ECG_CombF_dComp);
title('Filtered nECG signal with delay compensation');
xlim([0 2])
xlabel('Time (s)');
ylabel('Amplitude (mV)');
legend('Noisy ECG', 'Delay Compensated ECG');
grid on;
[P1,f1] = periodogram(nECG, window, N, fs);  % PSD estimate of nECG
[P2,f2] = periodogram(ECG_CombF_dComp, window, N, fs);   % PSD estimate of delay compensated ECG

figure();
plot(f1, 10*log10(P1), f2, 10*log10(P2));   % Plotting the PSD of noisy ECG and delay compensated ECG
title('PSD estimate of nECG vs delay compensated ECG');
xlabel('Frequency (Hz)');
ylabel('Power/frequency (dB/Hz)');
legend('Noisy ECG', 'Filtered Delay Compensated ECG');
grid on;