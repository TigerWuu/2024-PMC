%% Lab2: ILC
clc; clear; close all;

%% ---------------------------------------- Step 1 ------------------------------------
% Load your plant model (The model obtained from Lab1)
% Hint: you can use 'load' command to load your plant model P(z)
fs = 2e4;   %Sampling frequency
Ts = 1/fs;  %Sampling time
load("");

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
filename = sprintf('tri-%dHz-%ds-padding%.1fs.csv', f, Duration, padding);
writematrix(filename, r);

%% ---------------------------------------- Step 3 ------------------------------------
% Creat the Q Filter here
% Hint: 1.Use zplpf.m to design Q-filter, or design it directly by yourselves (e.g. with zpetc)
%		2.Use command [numQ, denQ] = tfdata(Q, 'v') to generate numerator vector of Q
% Record numQ and delayQ by yourselves

%% ---------------------------------------- Step 4 ------------------------------------
% Creat the L Filter here
% Hint: 1.Use zpetc.m to generate the L filter
% 		2.Use command [numL, denL] = tfdata(L, 'v') to generate numerator vector of L

% Save vector L
writematrix(numL, "ILC-numL.csv");
% Record delayQ by yourselves

%% ---------------------------------------- Step 5 ------------------------------------
% Simulate tracking results, analyze, and plot some stuff
