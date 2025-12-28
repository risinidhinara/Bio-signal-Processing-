%% 1. Clear workspace and load data
clear all; clc; close all;

load('ABR_rec.mat');   % ABR_rec(:,1)=stimulus, ABR_rec(:,2)=ABR+EEG
fs = 40000;            % 40 kHz

%% 2. Plot stimulus train and ABR train
figure('Position',[100 100 1200 600]);
plot(ABR_rec);
legend('Stimuli','ABR train');
xlabel('Samples'); ylabel('Amplitude (mV)');
title('ABR Recording (Stimulus + ABR+EEG)');

%% 3. Detect stimuli using threshold
thresh = find(ABR_rec(:,1) > 50);   % index where stimulus is >50 mV

% Extract only the first sample of each stimulus
j = 1;
for i = 1:length(thresh)-1
    if (thresh(i+1) - thresh(i)) > 1
        stim_point(j,1) = thresh(i+1); 
        j = j + 1;
    end
end

%% 4. Epoch extraction
% Window: -2ms → +10ms = total 12 ms
% Convert to samples: -2ms*fs = -80, +10ms*fs = +400, so window = 480 samples
preSamp  = 80;      % samples before stimulus
postSamp = 399;     % samples after stimulus

epochs = zeros(preSamp+postSamp+1, length(stim_point)); 

for i = 1:length(stim_point)
    epochs(:,i) = ABR_rec((stim_point(i)-preSamp):(stim_point(i)+postSamp),2);
end

%% 5. Ensemble averaging
ensmbl_avg = mean(epochs,2);

%% 6. Plot ensemble averaged ABR
time_ms = (-preSamp:postSamp)/fs*1000;   % convert to ms

figure;
plot(time_ms, ensmbl_avg, 'k','LineWidth',1.5);
xlabel('Time (ms)'); ylabel('Voltage (uV)');
title(['Ensemble averaged ABR from ', num2str(size(epochs,2)), ' epochs']);
grid on;

MSEk = progressiveMSE(ensmbl_avg, epochs, length(stim_point));


% Plotting the MSE values
figure('Position',[100 100 1200 600]);
plot(1 : length(stim_point),MSEk, '-');
xlim([-100 1100]);
ylim([-0.5 5]);
xlabel('Number of Epochs (k)');
ylabel('Mean Squared Error (MSE)');
title('MSE vs. Number of Epochs');
grid on

clear; clc; close all;

load('ECG_rec.mat');   % loads ECG_rec
fs = 128;              % Hz
ECG_rec = ECG_rec(:);  % ensure column vector
t = (0:length(ECG_rec)-1)/fs;   % time axis in seconds


figure();
plot(t, ECG_rec);   % Plotting ECG signal
title('Segment of the ECG Recording');
xlim([0 1.8])
ylim([-1 2])
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on

%% (iii) Extract one PQRST cycle as ECG_template
% (Here, pick an interval manually – adjust indices as needed)
% Example: suppose one PQRST lies between samples 200 and 360
idx_start = 230;
idx_end   = 340;
ECG_template = ECG_rec(idx_start:idx_end);

t_template = (0:length(ECG_template)-1)/fs;

figure;
plot(t_template, ECG_template);
xlabel('Time (s)'); ylabel('Amplitude (mV)');
title('Extracted PQRST waveform (ECG\_template)');
grid on;

%% (iv) Add Gaussian white noise (5 dB) to ECG_rec
nECG = awgn(ECG_rec, 5, 'measured');   % noisy ECG

figure;
plot(t,ECG_rec); hold on;
plot(t,nECG,'r');
xlim([0 1.8])
ylim([-1 2])
xlabel('Time (s)');
ylabel('Amplitude (mV)');
legend('Original ECG','Noisy ECG (5 dB)');
title('Comparison of ECG and noisy ECG');
grid on;

Segmenting ECG into separate epochs and ensemble averaging

% Normalize data
norm_ECG_template = (ECG_template - mean(ECG_template)) / std(ECG_template);
norm_nECG = (nECG - mean(nECG)) / std(nECG);

% Calculate normalized cross-correlation
[correlation, lags] = xcorr(norm_nECG, norm_ECG_template); 
lag_time = lags/fs;

%% (ii) Plot cross-correlation
figure('Position',[100 100 1200 600]);
plot(lag_time, correlation);
xlabel('Lag (s)'); ylabel('Normalised correlation');
title('Cross-correlation between template and noisy ECG');
xlim([0 60])
grid on;

%% (iii) Segment ECG pulses using correlation peaks

% Define threshold for correlation peaks
threshold = 0.6;   % adjust if needed (0.6–0.8 works well)

% Find peaks in correlation above threshold
[~, locs] = findpeaks(correlation, ...
                      'MinPeakHeight', threshold, ...
                      'MinPeakDistance', fs*0.5);

% Convert correlation peak locations (indices in correlation array) to lag values
startPoints = lags(locs);

% Initialize epoch matrix
epochLen = length(ECG_template);
ECG_epochs = [];

% Extract ECG pulses aligned with template
for i = 1:length(startPoints)
    idx = startPoints(i);
    if idx > 0 && (idx + epochLen - 1) <= length(nECG)
        ECG_epochs(:, end+1) = nECG(idx : idx+epochLen-1);
    end
end

%% Plot first 5 segmented ECG epochs
figure;
plot(ECG_epochs(:,1:5));
title('First 5 segmented ECG epochs');
xlabel('Samples (n)');
ylabel('Amplitude (mV)');
legend('Epoch 1','Epoch 2','Epoch 3','Epoch 4','Epoch 5');
grid on;

%% (iv) Calculate and plot SNR improvement with ensemble averaging

% Number of available ECG epochs
M = size(ECG_epochs,2);

% Preallocate array for progressive SNR values
progressiveSNRs = zeros(1, M);

% Ensure ECG_template is a column vector
ECG_template = ECG_template(:);

for k = 1:M
    % Ensemble average of first k epochs
    ECG_ensmbl_avg = mean(ECG_epochs(:,1:k), 2);

    % Noise = difference from template
    noise = ECG_ensmbl_avg - ECG_template;

    % Power calculations
    signalPower = mean(ECG_ensmbl_avg .^ 2);
    noisePower  = mean(noise .^ 2);

    % SNR in dB
    progressiveSNRs(k) = 10 * log10(signalPower / noisePower);
end

% Plot progressive SNR
figure('Position',[100 100 1200 600]);
plot(1:M, progressiveSNRs, 'b-o','LineWidth',1.2);
xlabel('Number of Epochs (k)');
ylabel('Signal-to-Noise Ratio (dB)');
title('SNR vs. Number of Ensemble-Averaged Epochs');
grid on;

%% (v) Compare noisy ECG vs ensemble-averaged ECG (20 and 80 pulses)

% Pick one noisy pulse (e.g., 3rd)
noisy_pulse = ECG_epochs(:,3);

% Compute ensemble averages
avg20 = mean(ECG_epochs(:,1:20), 2);
avg80 = mean(ECG_epochs(:,1:75), 2);

% Time axis (in seconds)
N = size(ECG_epochs,1);
time_axis = (0:N-1)/fs;

% Plot comparison
figure;
plot(time_axis, noisy_pulse, 'r'); hold on;
plot(time_axis, avg20, 'b');
plot(time_axis, avg80, 'k');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
legend('Single noisy ECG pulse','Average of 20 pulses','Average of 75 pulses');
title('Comparison of Noisy ECG vs Ensemble Averaged ECG');
grid on;