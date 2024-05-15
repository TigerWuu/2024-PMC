%% Lab2: RC
clc; clear; close all;

%% ---------------------------------------- Step 1 ------------------------------------
% Load your plant model (The model obtained from Lab1)
% Hint: you can use 'load' command to load your plant model P(z)
fs = 2e4;	%Sampling frequency
Ts = 1/fs;	%Sampling time
load("");

%% ---------------------------------------- Step 2 ------------------------------------
% Generate the triangular wave reference signal with amplitude of 1 and duration of 2 seconds
% You can decide the true amplitude on LabView code (default 1 [deg])
Duration = 2;	% second
f = 400;		% Hz
padding = 0;	% second, front and back

t = (0:Ts:Duration)';
tri = sawtooth(2*pi*f*t + pi/2 + pi/50, 0.5);
t = (0:Ts:(Duration + padding*2))';
r = [zeros(fs*padding, 1); tri; zeros(fs*padding, 1)];

% Write reference to .csv file
filename = sprintf('tri-%dHz-%ds-padding%.1fs.csv', f, Duration, padding);
writematrix(filename, r);

%% ---------------------------------------- Step 3 ------------------------------------
% Creat the RC controller
% Hint: 1.Use zpetc.m file to creat feedforward filter F(z).
%		2.Use zplpf.m to design Q-filter, or design it directly by yourselves
%		2.Please generate Cr(z) as in ref[1].

%The controller must be saved as the variable Cr
%i.e.: Cr = 

%% ---------------------------------------- Step 4 ------------------------------------
% Save the controller
filename2 = sprintf('RC-ZPETC-%dHz-Q.dat', f);
tf2host(Cr, filename2);

%% ---------------------------------------- Step 5 ------------------------------------
% Simulate tracking results, analyze, and plot some stuff
