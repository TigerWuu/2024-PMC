clc;
clear; 
%% Feedforward Controller Design
%% Load the ILCFF data here and generate your FF
% where FFdelay should be 47.
fs = 2e4;
Ts = 1/fs;
FFdelay = 47;

load("sys_est.mat");
filename = 'imp-order46-padding100-4000Hz.csv';
r = readmatrix(filename);
t = (1:length(r)) * Ts;

z = tf('z',Ts);
% Q1 = (z^(-1)+2+z)/4;
% Q2 = (z^(-1)+4+z)/6;
% Q = Q1;
% [numQ, denQ] = tfdata(Q, 'v');
%% Construct learning filter F by ILC
iterations = 100;
u0 = zeros(length(r),1);
y = zeros(length(r),1);
alpha = 1;
for i = 1:iterations
    e = r-y;
    e = flipud(e);
    e = lsim(sys_est, e, t);
    a = flipud(e);
    
    u1 = u0 + alpha*a;
    % u1 = conv(u1, numQ/4, 'same');

    y = lsim(sys_est, u1, t);
    u0 = u1;
    plot(t,r,t,y);
    pause(0.1);
end

%% Construct feedforward controller(F) from u1
Cr = 0;
N = length(u1);
for k = -(N-1)/2:(N-1)/2
    Cr = u1(k+(N-1)/2+1)*z^(-k)+Cr;
end

numF = u1';
denF = [1, zeros(1, size(numF,2))];
F = tf(numF, denF, Ts);

%% Save
filename = sprintf('FF-ILCFF-TR-delta-0.02s-0.20deg-beta1.0-alpha1.00-delay%d.dat',FFdelay);
tf2host(Cr, filename);
save("Cr.mat", "Cr");
save("F.mat", "F");
