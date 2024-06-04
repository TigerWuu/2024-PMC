%% Lab2: ILC
clc; clear; close all;

%% ---------------------------------------- Step 1 ------------------------------------
% Load your plant model (The model obtained from Lab1)
% Hint: you can use 'load' command to load your plant model P(z)
fs = 2e4;   %Sampling frequency
Ts = 1/fs;  %Sampling time
load("sys_est.mat");

%% ---------------------------------------- Step 2 ------------------------------------
% Generate the triangular wave reference signal with amplitude of 1 and duration of 2 seconds
% You can decide the true amplitude on LabView code (default 1 [deg])
Duration = 2;	% second
f = 400;		% Hz
padding = 0.1;	% second, front and back

t = (0:Ts:Duration)';
tri = sawtooth(2*pi*f*t + pi/2 + pi/50, 0.5);
t = (0:Ts:(Duration + padding*2))';
r = [zeros(fs*padding, 1); tri; zeros(fs*padding, 1)];

% Write reference to .csv file
% filename = sprintf('tri-%dHz-%ds-padding%.1fs.csv', f, Duration, padding);
% writematrix(filename, r);

%% ---------------------------------------- Step 3 ------------------------------------
% Creat the Q Filter here
% Hint: 1.Use zplpf.m to design Q-filter, or design it directly by yourselves (e.g. with zpetc)
%		2.Use command [numQ, denQ] = tfdata(Q, 'v') to generate numerator vector of Q
% Record numQ and delayQ by yourselves
z = tf('z',Ts);
Q1 = (z^(-1)+2+z)/4;
Q2 = (z^(-1)+4+z)/6;
Q = Q1;

%% ---------------------------------------- Step 4 ------------------------------------
% Load the L Filter here
load("Cr.mat");
load("F.mat");
L = F;

[numQ, denQ] = tfdata(Q, 'v');
[numL, denL] = tfdata(L, 'v');
% Save vector L
% writematrix(numL, "ILC-numL.csv");
% Record delayQ by yourselves

%% ---------------------------------------- Step 5 ------------------------------------
% Simulate tracking results, analyze, and plot some stuff
iterations = 20;
u0 = zeros(length(r),1);
y = zeros(length(r),1);

for i = 1:iterations
    e = r-y;
    % a = lsim(L, e, t);
    % a = filter(numL, denL, e);
    a = conv(e, numL/sum(numL), 'same');
    u1 = u0 + a;

    u1 = conv(u1, numQ/4, 'same');

    y = lsim(sys_est, u1, t);
    u0 = u1;
    plot(t,r,t,y);
    pause(0.1);
end


