%2. Adaptive filtering



% 2.1 LMS method

clear; clc;

phi = pi/6; % Phase shift
snr = 15;   % SNR in dB units

N = 5000;   % Number of samples
t = linspace(0, 5, N)'; % Time vector with fs = 500 Hz
s = sawtooth(2 * pi * 2 * t(1:N,1), 0.5);  % Sawtooth signal of 0.5 width
n1 = 0.2 * sin(2 * pi * 50 * t(1:N/2,1) - phi);   % Sinusoid at 50Hz
n2 = 0.3 * sin(2 * pi * 100 * t(N/2+1:N,1) - phi);    % Sinusoid at 100Hz
nwg = s - awgn(s, snr, 'measured');   % Gaussian white noise

non_stationary_noise = [n1; n2];    % Concatenate n1 and n2 noise 
noise = non_stationary_noise + nwg;

x_in = noise + s;   % Input signal
r_n = 0.6* (nwg + sin(2*pi*50*t + pi/4) + sin(2*pi*100*t + pi/6));

% Save the x_in and r_n vectors to a .mat file
save('r_n', 'r_n'); save('x_in', 'x_in');

% (a) Implement the LMS method

order_lms = 21; % LMS filter order
step_lms = 0.01;    % Step size 
[~, ~, e_lms] = LMS_Filter(r_n, x_in, step_lms, order_lms);

% (b) Plot the signals 𝑦𝑖(𝑛), 𝑥(𝑛), 𝑒(𝑛) = 𝑦̂(𝑛) and the absolute error |𝑦𝑖(𝑛) − 𝑦̂(𝑛)|

% Plot the signal yi(n)
figure('Position', [100, 100, 1000, 800]);   % Set figure size
subplot(4, 1, 1);
plot(t, s);
title('Sawtooth Waveform y_{i}(n) with a Width of 0.5');
ylim([-2, 2]);
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;

% Plot the input signal x(n)
subplot(4, 1, 2);
plot(t, x_in);
title('Input Signal x(n)');
ylim([-2, 2]);
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;

% Plot the error signal e(n)
subplot(4, 1, 3);
plot(t, e_lms);
title('Error Signal e(n) of LMS Filter');
ylim([-2, 2]);
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;

% Plot the absolute error signal |yi(n) - e(n)|
subplot(4, 1, 4);
plot(t, abs(s - e_lms));
title('Absolute Error |y_{i}(n) - e(n)| of LMS Filter');
ylim([0, 0.4]);
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;

% (c) Explore the rate of adaptation by varying the rate of convergence 𝜇 and the order of the adaptive filter (𝑀 − 1)

% Define the range of order and step to test
order_range = 1:1:20; % Different filter orders
step_range = 0.001:0.001:0.1;  % Step sizes from 0.001 to 0.1

% Initialize the matrix to store MSE values
MSE_matrix = zeros(length(step_range), length(order_range));

% Loop over all combinations of N and L
for i = 1 : length(order_range)
    for j = 1 : length(step_range)
        [~, ~, e_lms] = LMS_Filter(r_n, x_in, step_range(j), order_range(i));
        MSE_matrix(j, i) = mean((s - e_lms).^2);
    end
end

% Plot the MSE matrix
figure();
surf(order_range, step_range, MSE_matrix);
xlabel('Filter Order M');
ylabel('Step Size \mu');
zlabel('MSE');
title('MSE for Various LMS Filter Parameters');
colorbar;
grid on;

% Find the minimum MSE and corresponding indices
[minMSE, idx] = min(MSE_matrix(:));
[j, i] = ind2sub(size(MSE_matrix), idx);
optOrder_lms = order_range(i);
optStep = step_range(j);
fprintf('The optimum filter parameters which gives the minimum MSE are for M = %d and \x03BC = %f. It gives an MSE value of %f', optOrder_lms, optStep, minMSE);

[~, ~, e_lms] = LMS_Filter(r_n, x_in, optStep, optOrder_lms);

% Plot the error signal e(n)
figure('Position', [100, 100, 1000, 400]);
subplot(2, 1, 1);
plot(t, e_lms);
title('Error Signal e(n) of Optimal LMS Filter');
ylim([-2, 2]);
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;

% Plot the absolute error signal |yi(n) - e(n)|
subplot(2, 1, 2);
plot(t, abs(s - e_lms));
title('Absolute Error |y_{i}(n) - e(n)| of Optimal LMS Filter');
ylim([0, 0.4]);
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;



% 2.1 RLS method

% (a) Implement the RLS method
order_rls = 7;   % RLS filter order
lambda    = 0.95; % Forgetting factor

load('r_n'); load('x_in');

[~, ~, e_rls] = RLS_Filter(r_n, x_in, lambda, order_rls);

% (b) Compare the performance of LMS and RLS algorithms using a plot
figure('Position', [100, 100, 1000, 800]);
subplot(4,1,1);
plot(t, s); title('Sawtooth Waveform y_{i}(n)'); ylim([-2,2]);
xlabel('Time (s)'); ylabel('Amplitude (mV)'); grid on;

subplot(4,1,2);
plot(t, x_in); title('Input Signal x(n)'); ylim([-2,2]);
xlabel('Time (s)'); ylabel('Amplitude (mV)'); grid on;

subplot(4,1,3);
plot(t, e_rls); title('Error Signal e(n) of RLS Filter'); ylim([-2,2]);
xlabel('Time (s)'); ylabel('Amplitude (mV)'); grid on;

subplot(4,1,4);
plot(t, abs(s - e_rls));
title('|y_{i}(n) - e(n)| of RLS Filter'); ylim([0,0.4]);
xlabel('Time (s)'); ylabel('Amplitude (mV)'); grid on;

% (c) Explore λ and M for adaptation
order_range  = 1:1:30;
lambda_range = 0.97:0.001:1;
MSE_matrix   = zeros(length(lambda_range), length(order_range));

for i = 1:length(order_range)
    for j = 1:length(lambda_range)
        [~, ~, e_rls] = RLS_Filter(r_n, x_in, lambda_range(j), order_range(i));
        MSE_matrix(j,i) = mean((s - e_rls).^2);
    end
end

figure();
surf(order_range, lambda_range, MSE_matrix);
xlabel('Filter Order M'); ylabel('Forgetting Factor \lambda'); zlabel('MSE');
title('MSE for Various RLS Filter Parameters'); colorbar; grid on;

[minMSE, idx] = min(MSE_matrix(:));
[j,i] = ind2sub(size(MSE_matrix), idx);
optOrder_rls = order_range(i);
optLambda    = lambda_range(j);
fprintf('Optimum RLS: M = %d, λ = %.3f, MSE = %.6f\n', ...
    optOrder_rls, optLambda, minMSE);

[~, ~, e_rls] = RLS_Filter(r_n, x_in, optLambda, optOrder_rls);

figure('Position', [100,100,1000,400]);
subplot(2,1,1);
plot(t, e_rls); title('Error Signal e(n) of Optimal RLS Filter');
ylim([-2,2]); xlabel('Time (s)'); ylabel('Amplitude (mV)'); grid on;

subplot(2,1,2);
plot(t, abs(s - e_rls));
title('|y_{i}(n) - e(n)| of Optimal RLS Filter'); ylim([0,0.4]);
xlabel('Time (s)'); ylabel('Amplitude (mV)'); grid on;

% (d) Test LMS and RLS on idealECG
load('idealECG.mat');
idealECG = idealECG(1:5000);
fs = 500;
t_ecg = (0:N-1)'/fs;
x_in = idealECG' + noise;

[~, ~, ECG_LMS] = LMS_Filter(r_n, x_in, optStep, optOrder_lms);
[~, ~, ECG_RLS] = RLS_Filter(r_n, x_in, optLambda, optOrder_rls);

figure('Position',[100,100,1000,600]);
subplot(3,1,1); plot(t_ecg(1:500), idealECG(1:500)); title('Ideal ECG');
xlabel('Time (s)'); ylabel('Amplitude (mV)'); grid on;
subplot(3,1,2); plot(t_ecg(1:500), x_in(1:500)); title('Noisy ECG');
xlabel('Time (s)'); ylabel('Amplitude (mV)'); grid on;
subplot(3,1,3); plot(t_ecg(1:500), ECG_LMS(1:500)); title('LMS Filtered ECG');
xlabel('Time (s)'); ylabel('Amplitude (mV)'); grid on;

figure('Position',[100,100,1000,600]);
subplot(3,1,1); plot(t_ecg(1:500), idealECG(1:500)); title('Ideal ECG');
xlabel('Time (s)'); ylabel('Amplitude (mV)'); grid on;
subplot(3,1,2); plot(t_ecg(1:500), x_in(1:500)); title('Noisy ECG');
xlabel('Time (s)'); ylabel('Amplitude (mV)'); grid on;
subplot(3,1,3); plot(t_ecg(1:500), ECG_RLS(1:500)); title('RLS Filtered ECG');
xlabel('Time (s)'); ylabel('Amplitude (mV)'); grid on;

%% --- Local function inside script ---

function [w,y,e] = RLS_Filter(r_n, x_n, lambda, M)
    K = M + 1; N = length(x_n);
    P = 1e3*eye(K);       % large value for inverse correlation matrix
    w = zeros(K,1);       
    y = zeros(N,1);       
    e = zeros(N,1);       

    for n = K:N
        r = r_n(n:-1:n-K+1);             % input vector
        k = (P*r) / (lambda + r'*P*r);   % Kalman gain vector
        e(n) = x_n(n) - w'*r;            % a priori error
        w = w + k*e(n);                  % weight update
        P = (P - k*r'*P)/lambda;         % matrix inverse update
        y(n) = w'*r;                     % filter output (noise est.)
    end
end



%% --- Local LMS function inside script ---
function [w, y, e] = LMS_Filter(r_n, x_n, mu, M)
    K = M + 1;              % number of filter coefficients
    N = length(x_n);
    w = zeros(K,1);         % weight vector
    y = zeros(N,1);         % output
    e = zeros(N,1);         % error
    for n = K:N
        r = r_n(n:-1:n-K+1);    % input vector
        y(n) = w' * r;          % filter output
        e(n) = x_n(n) - y(n);   % error
        w = w + 2*mu*e(n)*r;    % weight update
    end
end
