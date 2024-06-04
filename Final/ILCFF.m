% Simulation validation for the manuscript entitled "Iterative Learning of Dynamic Inverse Filters for Feedforward Tracking Control"
% This is a simplified version (2019-7-29 by Cheng-Wei Chen)

clc; clear all; close all;
%% Simulation settings
SimIter       = 10000;              % Number of iterations performed to learn the target impulse
SimLength     = 10000;              % The length of the target impulse
SimDelay      = SimLength / 2;      % Delays to cover the anti-causal parts 

%% Prepare the controlled plant G(z)
Fs = 10000;                         % Control sampling rate [Hz]
Ts = 1/Fs;                         

P = zpk([1.01],[1 1],-0.01,Ts);     % Open-loop plant
C = 100*tf([0.2 -0.199],[1 0],Ts);  % Closed-loop controller

G = feedback(P*C,1);                % Closed-loop system
[numG, denG] = tfdata(G, 'v');     
pole(G)

%% Prepare the filter M(z) that shapes the bandwidth of inversion
Order = 200;                       % The order of the FIR filter N(z)
Fp  = 300;                         % The bandwidth of N(z)
Rp  = 0.00057565;                  % Corresponds to 0.01 dB peak-to-peak ripple
Rst = 1e-4;                        % Corresponds to 80 dB stopband attenuation
Ncoef = firceqrip(Order,Fp/(Fs/2),[Rp Rst],'passedge'); 
N = dsp.FIRFilter('Numerator',Ncoef); 

[numN, denN] = tf(N);
M = tf(numN,denN,Ts)*tf(numN,denN,Ts)'*tf([zeros(1, SimDelay),1],[1,zeros(1,SimDelay)],Ts);
M = M / dcgain(M);
[numM, denM] = tfdata(M, 'v');

%% Generate the target impulse r_M(k)
impM = impulse(M,SimLength*Ts) * Ts;  

H = zpk([],[],1/max(impM(2:end)),Ts);               % Scaling such that max(r_M)=1;
[numH, denH] = tfdata(H, 'v');

r_M = filter(numH, denH, fliplr(impM(2:end)))';
          

%% Learning the input u(k) which generates the target impulse r_M(k) at the output
u = zeros(1, SimLength);
y = zeros(1, SimLength);
e = zeros(1, SimLength);
Y = zeros(SimIter, SimLength);
U = zeros(SimIter, SimLength);
E = zeros(SimIter, SimLength);
E_max = zeros(1,SimIter);
E_rms = zeros(1,SimIter);
U_max = zeros(1,SimIter);
U_rms = zeros(1,SimIter);

disp('ILC Learning in Progress...');

% ILC algorithm
for ii = 1:SimIter   
    disp(['ILC iteration #' num2str(ii)]);
    
    U(ii, :) = u;
    
    y = filter(numG, denG, u);   % Inject the input to the controlled plant
    Y(ii, :) = y;
    
    e = r_M - y;   % Calculate the tracking error
    E(ii,:) = e;   
      
    u = u + 0.08*fliplr( filter(numG, denG, fliplr(e)) );   % Reversed-time based ILC update law
           
    E_max(ii) = max(abs(E(ii,:)));
    E_rms(ii) = sqrt(mean(E(ii,:).^2));
    U_max(ii) = max(abs(U(ii,:)));
    U_rms(ii) = sqrt(mean(U(ii,:).^2));
  
end


%% Construct the inverse filter F(z) from the learned input u(k)
numF = filter(denH, numH, u);
denF = [1, zeros(1, size(numF,2))];
F = tf(numF, denF, Ts);


%% Evaluation in frequency domain
fff = linspace(0,1/(2*Ts),10001);
www = fff*2*pi;

figure;

Freq_M = squeeze(freqresp(M,www));
Freq_G = squeeze(freqresp(G,www));
Freq_F = squeeze(freqresp(F,www));
Freq_Err = Freq_M - Freq_G.*Freq_F;

semilogx(fff,20*log10(abs(Freq_G)));hold on;
semilogx(fff,20*log10(abs(Freq_M)));
semilogx(fff,20*log10(abs(Freq_Err)));
ylabel('Magnitude (dB)');
xlim([1 5000]);xlabel('Frequency (Hz)');
legend('|G(z)|','|M(z)|','|M(z)-G(z)F(z)|')
grid on;


%% Evaluation by feedforward tracking
% 50 Hz triangular wave
RefFreq = 50;                         
T = 100*(1/RefFreq);
t = 0:1/Fs:T-1/Fs;
r = sawtooth(2*pi*RefFreq*t,1/2);

y_noF = filter(numG, denG, r);
u_withF = filter(numF, denF, r);
y_withF = filter(numG, denG, [u_withF(SimDelay:end) zeros(1,SimDelay-1)]);

figure;
plot(t,r); hold on
plot(t,y_noF);
plot(t,y_withF);
xlabel('Time [s]');xlim([t(end/2) t(end/2)+10*(1/RefFreq)]);
ylabel('Position output');
legend('r','y (without feedforward)','y (with feedforward)');
grid on;

%% Generating Fig. 4c
Color = [251  180  185;
         247  104  161;
         174  1    126;
         0    0    0    ]/255;
     
numF1 = filter(denH, numH, U(10, :));
numF2 = filter(denH, numH, U(100, :));
numF3 = filter(denH, numH, U(1000, :));
numF4 = filter(denH, numH, U(10000, :));

denF1 = [1, zeros(1, size(numF1,2))];
denF2 = [1, zeros(1, size(numF2,2))];
denF3 = [1, zeros(1, size(numF3,2))];
denF4 = [1, zeros(1, size(numF4,2))];

F1 = tf(numF1, denF1, Ts);  
F2 = tf(numF2, denF2, Ts);  
F3 = tf(numF3, denF3, Ts);  
F4 = tf(numF4, denF4, Ts);  
     
Freq_F1 = squeeze(freqresp(F1,www));
Freq_F2 = squeeze(freqresp(F2,www));
Freq_F3 = squeeze(freqresp(F3,www));
Freq_F4 = squeeze(freqresp(F4,www));

Freq_Err1 = Freq_M - Freq_G.*Freq_F1;
Freq_Err2 = Freq_M - Freq_G.*Freq_F2;
Freq_Err3 = Freq_M - Freq_G.*Freq_F3;
Freq_Err4 = Freq_M - Freq_G.*Freq_F4;

figure;
semilogx(fff,20*log10(abs(Freq_Err1)),'Color',Color(1,:),'Linewidth',2);hold on;
semilogx(fff,20*log10(abs(Freq_Err2)),'Color',Color(2,:),'Linewidth',2);
semilogx(fff,20*log10(abs(Freq_Err3)),'Color',Color(3,:),'Linewidth',2);
semilogx(fff,20*log10(abs(Freq_Err4)),'Color',Color(4,:),'Linewidth',2);
title('|M(z)-G(z)F(z)|');
ylabel('Magnitude (dB)');
xlim([1 5000]);xlabel('Frequency (Hz)');
legend('Case 1: 10 iterations','Case 2: 100 iterations','Case 3: 1,000 iterations','Case 4: 10,000 iterations');
grid on;