clc; clear; 
%% Feedforward Controller Design
%% Load the ILCFF data here and generate your FF
% where FFdelay should be 47.
fs = 2e4;
Ts = 1/fs;
FFdelay = 47;
Cr= ...

%% save
filename = sprintf('FF-ILCFF-TR-delta-0.02s-0.20deg-beta1.0-alpha1.00-delay%d.dat',FFdelay);
tf2host(Cr, filename);


