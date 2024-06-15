clc; clear; close all;

%% Reference Model Design(The impulse signal convolution with lowpass filter)
% zero phase low pass filter design
Fs = 2e4; Ts = 1/Fs;
amp = 1;

SimBandwidth = 400;		% passband-edge frequency
SimPadding = 100;

N = 46;						% FIR filter order
Rp  = 10^(0.01/2/20)-1;		% Corresponds to 0.01 dB peak-to-peak ripple
Rst = 1e-4;					% Corresponds to 80 dB stopband attenuation

NUM = firceqrip(N, SimBandwidth/(Fs/2), [Rp, Rst], 'passedge')'; % SimBandwidth/(Fs/2) is normalized frequency; NUM = vector of coeffs
% fvtool(NUM)

%% Analyze NUM
FIR = tf(NUM', [1, zeros(1, 47)], Ts) * tf(NUM', [1, zeros(1, 47)], Ts)';
figure; bode(FIR);
% disp(bandwidth(FIR)/2/pi);
Q = tf([1 2 1], [4 0], Ts);
figure; bode((1-FIR)*Q);

%% reference impulse
imp = conv(NUM, NUM);
r = amp * imp / max(imp);

% padding
r = [zeros(SimPadding, 1); r; zeros(SimPadding, 1)];

t = (1:length(r)) * Ts;
figure;
subplot(3, 1, 1);
plot(t, r, '.-');
xlabel('Time [s]'); ylabel('Pos [mm]');
grid on;
subplot(3, 1, 2);
plot((t(1:end-1)+t(2:end))/2, diff(r)/Ts, '.-');
xlabel('Time [s]'); ylabel('Vel [mm/s]');
grid on;
subplot(3, 1, 3);
plot(t(2:end-1), diff(r,2)/Ts^2, '.-');
xlabel('Time [s]'); ylabel('Acc [mm/s^2]');
grid on;

%% Write the reference impulse to "ref_imp.csv"
% filename = sprintf('imp-order%d-padding%d-%dHz.csv', N, SimPadding, SimBandwidth);
% csvwrite(filename, r);
