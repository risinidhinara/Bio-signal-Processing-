%% BM4152 Assignment 3 - Section 1.2 (iv)
% wavelet_construction.m



% Section 1.2 (iv & v) - Mexican hat daughter wavelets


clearvars; close all; clc;


%% Parameters
fs = 250;                   % Sampling frequency
N = 3000;                   % Data length
t = (-N:N)/fs;              % Time vector
s = 0.01:0.1:2;             % Scaling factors

%% Mexican hat wavelet 
m = @(tt) (1/sqrt(2*pi))*(1 - tt.^2).*exp(-tt.^2/2);
E = trapz(t,m(t).^2);
m_norm = @(tt) m(tt)/sqrt(E);

%% Plot all daughter wavelets in one figure
figure; hold on;
set(gcf,'Position',[100 100 1200 600]);
for k = 1:length(s)
    wavelt = (1/sqrt(s(k))) * m_norm(t/s(k));
    plot(t,wavelt,'DisplayName',['s=' num2str(s(k))]);
end
xlabel('Time (s)'), ylabel('\psi_s(t)')
title('Mexican Hat Daughter Wavelets (Time Domain)')
xlim([-8 8]);
ylim([-3 4]);
legend show; grid on; hold off;


%% Verification table: mean, energy, effective support
means = zeros(size(s));
energies = zeros(size(s));
support_width = zeros(size(s));

for k = 1:length(s)
    wavelt = (1/sqrt(s(k))) * m_norm(t/s(k));
    
    % Zero mean & unity energy
    means(k) = trapz(t,wavelt);        
    energies(k) = trapz(t,wavelt.^2);  
    
    % Compact support: width where |ψ| > threshold
    thresh = 1e-3 * max(abs(wavelt));
    idx = find(abs(wavelt) > thresh);
    support_width(k) = t(idx(end)) - t(idx(1));
end

% Round means to 2 decimal places for display
means = round(means, 2);

T = table(s', means', energies', support_width', ...
          'VariableNames', {'Scale','Mean','Energy','SupportWidth'});
disp(T)

%% Part (vi) - Frequency spectra of Mexican Hat daughter wavelets

figure('Name','Full Spectrum: Mexican Hat Daughter Wavelets','NumberTitle','off'); 
set(gcf,'Position',[100 100 1200 600]);
hold on;

for k = 1:length(s)
    wavelt = (1/sqrt(s(k))) * m_norm(t/s(k));   % Daughter wavelet
    
    % FFT-based spectrum
    Fwavelt = fft(wavelt)/length(wavelt);
    hz = linspace(0, fs/2, floor(length(wavelt)/2)+1);
    plot(hz, 2*abs(Fwavelt(1:length(hz))), 'DisplayName',['s=' num2str(s(k))]);
end

xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('Full Frequency Spectra of Mexican Hat Daughter Wavelets')
legend show, grid on; hold off;
xlim([0 80]);


% Zoomed-in spectrum plot (0-10 Hz)
figure('Name','Zoomed Spectrum: Mexican Hat Daughter Wavelets','NumberTitle','off'); 
set(gcf,'Position',[150 150 1200 600]);
hold on;

for k = 1:length(s)
    wavelt = (1/sqrt(s(k))) * m_norm(t/s(k));   % Daughter wavelet
    
    % FFT-based spectrum
    Fwavelt = fft(wavelt)/length(wavelt);
    hz = linspace(0, fs/2, floor(length(wavelt)/2)+1);
    plot(hz, 2*abs(Fwavelt(1:length(hz))), 'DisplayName',['s=' num2str(s(k))]);
end

xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('Zoomed Frequency Spectra (0–10 Hz) of Mexican Hat Daughter Wavelets')
legend show, grid on; hold off;
xlim([0 6]);


%% Parameters  
n = 1:3*N;         
t = n/fs;          

%% Create waveform using actual time
x = zeros(size(n));

% First segment: sin(0.5*pi*n/fs)
x(1:3*N/2) = sin(0.5*pi*(n(1:3*N/2)/fs));

% Second segment: sin(1.5*pi*n/fs)
x(3*N/2+1:end) = sin(1.5*pi*(n(3*N/2+1:end)/fs));

%% Plot the waveform
figure;
set(gcf,'Position',[150 150 1200 600]);
plot(t, x);
xlabel('Time (s)');
ylabel('Amplitude');
title('Piecewise Signal x[n]');
grid on;
%% Mexican Hat mother wavelet (normalized)
wavelet_mother = @(tt) (1/sqrt(2*pi))*(1 - tt.^2).*exp(-tt.^2/2);
E = trapz((-N:N)/fs, wavelet_mother((-N:N)/fs).^2);   % Energy
wavelet_norm = @(tt) wavelet_mother(tt)/sqrt(E);       % Normalized
s_values = 0.01:0.01:2; 
%% Continuous Wavelet Transform (CWT) using convolution
coeffs = zeros(length(s_values), length(x));

for k = 1:length(s_values)
    s = s_values(k);
    wavelt = (1/sqrt(s)) * wavelet_norm((-N:N)/(fs*s));  % Scaled wavelet
    coeffs(k,:) = conv(x, wavelt, 'same');               % Convolution for translation
end

%% Plot Scale-Time Spectrogram
figure('Name','CWT: Mexican Hat','NumberTitle','off');
h = pcolor(t, s_values, coeffs);
set(h, 'EdgeColor', 'none');
colormap jet;
xlabel('Time (s)');
ylabel('Scale');
title('Continuous Wavelet Transform (Mexican Hat)');
colorbar;
set(gca,'YDir','normal');  % Ensure scale increases upwards