clc;
clear; 
%% initial
fs = 2e4;
Ts = 1/fs;
FFdelay = 47;

z = tf('z',Ts);
load("sys_est.mat");
%% Feedforward Controller Design (from experiments)
hz = 4000;
% filename = './experiments/0604/4000Hz/Group5-2024-06-04-14-01-04-TR-(imp-order46-padding100-4000Hz)-1.0deg-beta1.0-alpha0.10-Q121-iter100.csv';
filename = './experiments/0608/2000Hz/Group5-2024-06-08-12-49-29-TR-(imp-order46-padding100-2000Hz)-1.0deg-beta1.0-alpha0.50-Q000-iter75.csv';
data = readmatrix(filename);
% u = data(:,2);
u = data(100:194,2); % cause the FFdelay is desiged as 47
% filename = 'imp-order46-padding100-4000Hz.csv';
filename = 'imp-order46-padding100-2000Hz.csv';
r = readmatrix(filename);



Cr = 0;
N = length(u);
for k = -(N-1)/2:(N-1)/2
    Cr = u(k+(N-1)/2+1)*z^(-k)+Cr;
end
Cr = Cr/sum(r(100:194));

numF = u';
denF = [1, zeros(1, size(numF,2))];
F = tf(numF, denF, Ts);
F = F/sum(r(100:194));
% save data
% filename = sprintf('FF-ILCFF-TR-delta-0.02s-0.20deg-beta1.0-alpha0.50-delay%d-exp-%d-normalized-F.dat',FFdelay, hz);
% tf2host(F, filename);

save("Cr.mat", "Cr");
% save("F.mat", "F");
%% Load the ILCFF data here and generate your FF
% where FFdelay should be 47.


t = (1:length(r)) * Ts;

Q1 = (z^(-1)+2+z)/4;
Q2 = (z^(-1)+4+z)/6;
Q = Q1;
[numQ, denQ] = tfdata(Q, 'v');
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
    % pause(0.01);
end



%% Construct feedforward controller(F) from u1
u1 = u1(100:194); % cause the FFdelay is desiged as 47

Cr = 0;
N = length(u1);
for k = -(N-1)/2:(N-1)/2
    Cr = u1(k+(N-1)/2+1)*z^(-k)+Cr;
end
Cr = Cr/sum(r(100:194));

numF = u1'; %% normalized ?
% numF = numF/sum(numF);
denF = [1, zeros(1, size(numF,2))];
F = tf(numF, denF, Ts);
F = F/sum(r(100:194));
% F = F/dcgain(F);
% Save
% filename = sprintf('FF-ILCFF-TR-delta-0.02s-0.20deg-beta1.0-alpha0.50-delay%d-sim-%d-normalized-F.dat',FFdelay, hz);
% tf2host(F, filename);

save("Cr.mat", "Cr");
% save("F.mat", "F");
