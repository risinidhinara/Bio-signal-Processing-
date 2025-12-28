clear; clc;
load('ECG_template.mat');
x = ECG_template(:);        
fs = 500;                   % Hz, given in assignment
t = (0:length(x)-1)/fs;

figure('Position',[100 100 1200 600]);
plot(t,x,'b');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
title('ECG Template Signal (fs = 500 Hz)');
grid on;

nECG = awgn(x,5,'measured');

figure('Position',[100 100 1200 600]);
plot(t,nECG,'r');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
title('Noisy ECG (5 dB AWGN)');
grid on;

figure('Position',[100 100 1200 600]);
periodogram(nECG,[],[],fs);
title('PSD of Noisy ECG (nECG)');


N = 3;  % Filter length

% Zero-padding the signal
xPadded = [zeros(N-1,1); nECG];

% Initialize the filtered output
ma3ECG_1 = zeros(size(nECG));

% Apply the moving average filter
for n = N:length(xPadded)
    ma3ECG_1(n-(N-1)) = sum(xPadded(n-N+1 : n)) / N;
end

% Plot the filtered ECG signal
figure('Position',[100 100 1200 600]);
plot(t, nECG, 'r'); hold on;
plot(t,x,'k--','LineWidth',1);
plot(t, ma3ECG_1, 'b');
xlabel('Time [s]');
ylabel('Amplitude [mV]');
legend('Noisy ECG','Original ECG','MA(3) Filtered');
title('Comparison of Original, Noisy, and MA(3) Filtered ECG');
grid on;

groupDelay = ceil((N-1)/2);
fprintf('The group delay for MA(3) is %d', groupDelay);


ma3ECG_1_shifted = [ma3ECG_1(groupDelay+1:end); zeros(groupDelay,1)];

% Plot comparison
figure('Position',[100 100 1200 600]);
plot(t, nECG, 'r'); hold on;
plot(t,x,'k--','LineWidth',1);
plot(t, ma3ECG_1_shifted, 'b');
xlabel('Time [s]');
ylabel('Amplitude [mV]');
legend('Noisy ECG', 'Original ECG','MA(3) Filtered (Delay Compensated)');
title('Comparison of nECG and MA(3) Filtered ECG (Delay Compensated)');

[pxx1,f1] = periodogram(nECG, [], [], fs);
[pxx2,f2] = periodogram(ma3ECG_1, [], [], fs);

figure('Position',[100 100 1200 600]);
plot(f1, 10*log10(pxx1), 'r'); hold on;
plot(f2, 10*log10(pxx2), 'b');
xlabel('Frequency [Hz]');
ylabel('Power/Frequency [dB/Hz]');
legend('Noisy ECG', 'MA(3) Filtered ECG');
title('PSD Comparison: nECG vs MA(3) ECG');
grid on;

N3 = 3;
b3 = ones(1,N3)/N3;  % numerator coefficients
a = 1;               % denominator
ma3ECG_2 = filter(b3,a,nECG);

delay3 = (N3-1)/2;
ma3ECG_2 = [ma3ECG_2(delay3+1:end); zeros(delay3,1)];


% (ii) Plot nECG, ECG_template, ma3ECG_2
figure('Position',[100 100 1200 600]);
plot(t,nECG,'r'); hold on;
plot(t,x,'k--','LineWidth',1);
plot(t,ma3ECG_2,'b');
xlabel('Time (s)'); ylabel('Amplitude (mV)');
legend('Noisy ECG','Original ECG Template','MA(3) Filtered ECG (ma3ECG\_2)');
title('Comparison of Noisy, Original and MA(3) filtered ECG');

% (iii) Inspect MA(3) filter with fvtool
% Magnitude response
h1 = fvtool(b3,a,'Analysis','magnitude');
set(h1,'WindowStyle','normal');
set(h1,'Position',[100 400 250 250]); % [x y width height]


% Phase response
h2 = fvtool(b3,a,'Analysis','phase');
set(h2,'WindowStyle','normal');
set(h2,'Position',[550 400 250 250]);


% Pole-zero plot
h3 = fvtool(b3,a,'Analysis','polezero');
set(h3,'WindowStyle','normal');
set(h3,'Position',[1000 400 250 250]);


N10 = 10;
b10 = ones(1,N10)/N10;
a = 1;  
ma10ECG = filter(b10,a,nECG);

% Compensate delay
delay10 = floor((N10-1)/2); 
ma10ECG = [ma10ECG(delay10+1:end); zeros(delay10,1)];
fvtool(b3,a,b10,a);


fvtool(b3,a,b10,a,'Analysis','polezero');

fvtool(b3,a,b10,a,'Analysis','phase');

figure('Position',[100 100 1200 600]);
plot(t,nECG,'r'); hold on;
plot(t,x,'k--');
plot(t,ma3ECG_2,'b');
plot(t,ma10ECG,'c','LineWidth',1);
xlabel('Time (s)'); ylabel('Amplitude (mV)');
legend('Noisy ECG','Original ECG','MA(3) Filtered','MA(10) Filtered');
title('Comparison of Noisy, Template, MA(3) and MA(10)');

N_range = 1:50;               % test from N=2 to N=50
MSE_values = zeros(size(N_range));

for N = 1 : 50
    
    b= ones(1,N)/N;
    a = 1;  
    maECG = filter(b,a,nECG);
    
    % Compensate delay
    delay = floor((N-1)/2); 
    maECG = [maECG(delay+1:end); zeros(delay,1)];
    MSE_values(N) = mean((x - maECG).^2);
end

%% Plot MSE vs. N
figure('Position',[100 100 1200 600]);
plot(N_range, MSE_values, 'b-o','LineWidth',1.2);
xlabel('Window length N');
ylabel('MSE');
title('MSE vs MA(N) Filter Order');
grid on;

% Find optimum N
[~, idx] = min(MSE_values);
optN = N_range(idx);
disp(['Optimum MA filter order: N = ', num2str(optN)]);


%% Apply SG(3,11) filter
N = 3;       % polynomial order
L = 5;       % half-window
Lprime = 2*L + 1;   % window length = 11

sg310ECG = sgolayfilt(nECG, N, Lprime);



%% Plot nECG, ECG_template, and sg310ECG
figure('Position',[100 100 1200 600]);
plot(t, x, 'k--', 'LineWidth',1); hold on;   % original ECG
plot(t, nECG, 'r');                          % noisy ECG
plot(t, sg310ECG, 'b'); % SG(3,11) filtered ECG
xlabel('Time [s]');
ylabel('Amplitude [mV]');
legend('Original ECG','Noisy ECG','SG(3,11) Filtered ECG');
title('Comparison of ECG signals with SG(3,11) filter');
grid on;

%% Search for optimum SG(N,L)
N_values = 1:20;     % polynomial orders to test
L_values = 1:30;    % half-window sizes (window length = 2L+1)

MSE_sg = zeros(length(N_values), length(L_values));

for i = 1:length(N_values)
    for j = 1:length(L_values)
        N_poly = N_values(i);
        L_half = L_values(j);
        Lprime = 2*L_half + 1;

        % Ensure valid condition N <= L'-1
        if N_poly <= Lprime - 1
            y_sg = sgolayfilt(nECG, N_poly, Lprime);
            % Compute MSE against clean ECG
            MSE_sg(i,j) = mean((x - y_sg).^2);
        else
            MSE_sg(i,j) = NaN; % invalid region
        end
    end
end

%% Plot MSE surface
figure;
surf(L_values, N_values, MSE_sg);
xlabel('Half-window length L');
ylabel('Polynomial order N');
zlabel('MSE');
title('MSE Surface for SG(N,L) filter');

[minMSE, idx] = min(MSE_sg(:));
[i_opt, j_opt] = ind2sub(size(MSE_sg), idx);

optN = N_values(i_opt);
optL = L_values(j_opt);
optLprime = 2*optL + 1;

disp(['Optimum SG filter: N = ', num2str(optN), ', L = ', num2str(optL), ...
      ' (Window length L'' = ', num2str(optLprime), ')']);

%% Apply SG(2,27) filter
N = 2;       % polynomial order
L = 13;       % half-window
Lprime = 2*L + 1;   % window length = 27

sgoptimumECG = sgolayfilt(nECG, N, Lprime);



%% Plot nECG, ECG_template, and sg310ECG
figure('Position',[100 100 1200 600]);
plot(t, x, 'k--', 'LineWidth',1); hold on;   % original ECG
plot(t, nECG, 'r');                          % noisy ECG
plot(t, ma12ECG, 'b','LineWidth',1); % SG(3,11) filtered ECG
plot(t, sgoptimumECG, 'g','LineWidth',1);
xlabel('Time [s]');
ylabel('Amplitude (mV)');
legend('Original ECG','Noisy ECG','MA(12) Filtered ECG','SG optimal Filtered ECG');
title('Comparison of ECG signals with MA and SC optimal filters');
grid on;

N12 = 12;
b12 = ones(1,N12)/N12;  % numerator coefficients
a = 1;               % denominator
ma12ECG = filter(b12,a,nECG);

delay12 = floor((N12-1)/2);
ma12ECG = [ma12ECG(delay12+1:end); zeros(delay12,1)];




%% Plot nECG, ECG_template, and sg310ECG
figure('Position',[100 100 1200 600]);
plot(t, x, 'k--', 'LineWidth',1); hold on;   % original ECG
plot(t, nECG, 'r');                          % noisy ECG
plot(t, ma12ECG, 'b','LineWidth',1); % SG(3,11) filtered ECG
plot(t, sgoptimumECG, 'g','LineWidth',1);
xlabel('Time (s)');
ylabel('Amplitude (mV)');
legend('Original ECG','Noisy ECG','MA(12) Filtered ECG','SG optimal Filtered ECG');
title('Comparison of ECG signals with MA and SC optimal filters');
grid on;