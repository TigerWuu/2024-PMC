close all;
clear;
clc;
%% load data
FFdelay = 47;
% F_sim_2000 = load("./filters/F_sim_2000.mat").F;
% F_exp_2000 = load("./filters/F_exp_2000.mat").F;
% F_sim_4000 = load("./filters/F_sim_4000.mat").F;
% F_exp_4000 = load("./filters/F_exp_4000.mat").F;
% F_sim_2000_un = load("./filters/F_sim_2000_un.mat").F;
% F_exp_2000_un = load("./filters/F_exp_2000_un.mat").F;
% F_sim_4000_un = load("./filters/F_sim_4000_un.mat").F;
% F_exp_4000_un = load("./filters/F_exp_4000_un.mat").F;
F_sim_2000 = load("./filters/Cr_sim_2000.mat").Cr;
F_exp_2000 = load("./filters/Cr_exp_2000.mat").Cr;
F_sim_4000 = load("./filters/Cr_sim_4000.mat").Cr;
F_exp_4000 = load("./filters/Cr_exp_4000.mat").Cr;
F_sim_2000_un = load("./filters/Cr_sim_2000_un.mat").Cr;
F_exp_2000_un = load("./filters/Cr_exp_2000_un.mat").Cr;
F_sim_4000_un = load("./filters/Cr_sim_4000_un.mat").Cr;
F_exp_4000_un = load("./filters/Cr_exp_4000_un.mat").Cr;

load("sys_est.mat");
fs = 2e4;   %Sampling frequency
Ts = 1/fs;  %Sampling time
% reference setup
Duration = 2;
padding = 1;	% second, front and back
f = 400;		% Hz

% low pass filter
z = tf('z',Ts);
Q = (z^(-1)+2+z)/4;
[numQ, denQ] = tfdata(Q, 'v');
%% plot style
titlefont = 20;
tickfont = 18;
legendfont = 16;

plotoptions = bodeoptions;
plotoptions.Grid = 'on';
plotoptions.FreqUnits = 'Hz';   
% plotoptions.Title.String = 'Bode Plot (Complex Curve Fitting)';
plotoptions.Title.FontSize = titlefont;
plotoptions.Title.FontWeight = "bold";
plotoptions.Xlabel.FontSize= tickfont;
plotoptions.Ylabel.FontSize= tickfont;
plotoptions.Ticklabel.FontSize= 10;
%% bode plot of inverse filter : 4 cases
figure();
hold on;
bodemag(F_sim_2000_un*sys_est, plotoptions);
bodemag(F_exp_2000_un*sys_est, plotoptions);
bodemag(F_sim_4000_un*sys_est, plotoptions);
bodemag(F_exp_4000_un*sys_est, plotoptions);
legend("$F_{sim}:2k$","$F_{exp}:2k$","$F_{sim}:4k$","$F_{exp}:4k$",'interpreter','latex','fontsize',legendfont,'location','best');

figure();
hold on;
bodemag(F_sim_2000*sys_est, plotoptions);
bodemag(F_exp_2000*sys_est, plotoptions);
bodemag(F_sim_4000*sys_est, plotoptions);
bodemag(F_exp_4000*sys_est, plotoptions);
legend("$F_{sim}:2k$","$F_{exp}:2k$","$F_{sim}:4k$","$F_{exp}:4k$",'interpreter','latex','fontsize',legendfont,'location','best');
% figure();
% plot(t, sine_r, t, sine_y, t, sine_y_hat,'-','linewidth',1.5); hold on;
% xlabel('Times[s]','interpreter','latex','fontsize',tickfont);
% ylabel('Magnitude[m]','interpreter','latex','fontsize',tickfont);
% title('400Hz sinewave tracking', 'fontsize',titlefont);
% legend("$r$","$y$","$\hat{y}$",'interpreter','latex','fontsize',legendfont,'location','best')
% grid on; grid minor;
%% stability analysis
filters = [F_sim_2000, F_exp_2000, F_sim_4000, F_exp_4000];
for fil=filters
    bodemag(Q*(1-fil*sys_est),plotoptions); hold on;
end
title('$Q(z)|1-F(z)G(z)|$','interpreter','latex','fontsize',titlefont);
legend("$F(z)=F_{sim}:2k$","$F(z)=F_{exp}:2k$","$F(z)=F_{sim}:4k$","$F(z)=F_{exp}:4k$",'interpreter','latex','fontsize',legendfont,'location','best')
grid on; grid minor;
ylim([-128 5])


%% data-based feedforward 
% experiment
filename = './experiments/0610/Group5-2024-06-10-17-43-43-FF-(FF-ILCFF-TR-delta-0.02s-0.20deg-beta1.0-alpha0.50-delay47-sim-2000-normalized-F).csv';
data_sim_2000 = readmatrix(filename);
filename = './experiments/0610/Group5-2024-06-10-11-14-04-FF-(FF-ILCFF-TR-delta-0.02s-0.20deg-beta1.0-alpha0.50-delay47-exp-2000-normalized-F).csv';
data_exp_2000 = readmatrix(filename);
filename = './experiments/0610/Group5-2024-06-10-17-44-57-FF-(FF-ILCFF-TR-delta-0.02s-0.20deg-beta1.0-alpha0.50-delay47-sim-4000-normalized-F).csv';
data_sim_4000 = readmatrix(filename);
filename = './experiments/0610/Group5-2024-06-10-11-15-01-FF-(FF-ILCFF-TR-delta-0.02s-0.20deg-beta1.0-alpha0.50-delay47-exp-4000-normalized-F).csv';
data_exp_4000 = readmatrix(filename);
filename = './tri-400Hz-2s-padding1.0s.csv';
r = readmatrix(filename);
t = (0:Ts:(Duration + padding*2))'; 
y_sim_2000 = data_sim_2000(:,3);
y_exp_2000 = data_exp_2000(:,3);
y_sim_4000 = data_sim_4000(:,3);
y_exp_4000 = data_exp_4000(:,3);
% shift FFdelay?
t = t(1:end-FFdelay);
r = r(1:end-FFdelay);
y_sim_2000 = y_sim_2000(FFdelay+1:end);
y_exp_2000 = y_exp_2000(FFdelay+1:end);
y_sim_4000 = y_sim_4000(FFdelay+1:end);
y_exp_4000 = y_exp_4000(FFdelay+1:end);

figure();
plot(t, r, '-','linewidth',1.5); hold on;
plot(t, y_exp_2000, t, y_exp_4000, '-','linewidth',1.5);
xlabel('Times[s]','interpreter','latex','fontsize',tickfont);
ylabel('Magnitude[m]','interpreter','latex','fontsize',tickfont);
title('Data-based feedforward tracking', 'fontsize',titlefont);
legend("r","$y_{exp}:2k$","$y_{exp}:4k$",'interpreter','latex','fontsize',legendfont,'location','best')
grid on; grid minor;
xlim([0.9995 1.0025])
ylim([-1.1 1.1])

% simulation
F_sim_2000_d = F_sim_2000*z^(-FFdelay);
[numF, denF] = tfdata(F_sim_2000_d, 'v');
u = filter(numF, denF, r);
y_sim_2000 = lsim(sys_est, [u(FFdelay+1:end)' zeros(1,FFdelay)], t);
% u = conv(r, numF/sum(numF), 'same');
% y_sim_2000 = lsim(sys_est, u, t);
F_exp_2000_d = F_exp_2000*z^(-FFdelay);
[numF, denF] = tfdata(F_exp_2000_d, 'v');
u = filter(numF, denF, r);
y_exp_2000 = lsim(sys_est, [u(FFdelay+1:end)' zeros(1,FFdelay)], t);

F_sim_4000_d = F_sim_4000*z^(-FFdelay);
[numF, denF] = tfdata(F_sim_4000_d, 'v');
u = filter(numF, denF, r);
y_sim_4000 = lsim(sys_est, [u(FFdelay+1:end)' zeros(1,FFdelay)], t);

F_exp_4000_d = F_exp_4000*z^(-FFdelay);
[numF, denF] = tfdata(F_exp_4000_d, 'v');
u = filter(numF, denF, r);
y_exp_4000 = lsim(sys_est, [u(FFdelay+1:end)' zeros(1,FFdelay)], t);

figure();
plot(t, r, '-','linewidth',1.5); hold on;
plot(t, y_exp_2000,t, y_exp_4000, '-','linewidth',1.5);
xlabel('Times[s]','interpreter','latex','fontsize',tickfont);
ylabel('Magnitude[m]','interpreter','latex','fontsize',tickfont);
title('Data-based feedforward tracking', 'fontsize',titlefont);
legend("r","$y_{exp}:2k$","$y_{exp}:4k$",'interpreter','latex','fontsize',legendfont,'location','best')
grid on; grid minor;
xlim([0.9995 1.0025])
ylim([-1.1 1.1])

%% Model-free ILC
% experiment
files = ["exp-2000","exp-4000"];
iters = ["0", "1", "7"];
filename = './tri-400Hz-2s-padding1.0s.csv';
r = readmatrix(filename);
t = (0:Ts:(Duration + padding*2))'; 
for file = files
    figure()
    plot(t, r, '-','linewidth',1.5); hold on;
    for i=iters
        hold on;
        filename = './experiments/0610/Group5-(FF-ILCFF-TR-delta-0.02s-0.20deg-beta1.0-alpha0.50-delay47-'+file+'-normalized-F)/iter'+i+'.csv';
        data = readmatrix(filename);
        % t = data(:,1);
        y = data(:,3);
        plot(t, y, '-.','linewidth',1.5); hold on;
    end
    xlabel('Times[s]','interpreter','latex','fontsize',tickfont);
    ylabel('Magnitude[m]','interpreter','latex','fontsize',tickfont);
    title('Model-free ILC tracking', 'fontsize',titlefont);
    legend("r","iter: "+iters(1),"iter: "+iters(2),"iter: "+iters(3),'interpreter','latex','fontsize',legendfont,'location','best')
    % legend("r","$y_{exp}:2k$","$y_{exp}:4k$",'interpreter','latex','fontsize',legendfont,'location','best')

    grid on; grid minor;
    xlim([0.9995 1.0025])
    ylim([-1.1 1.1])
end

% simulation
F_exp_2000_d = F_exp_2000*z^(-FFdelay);
[numF_2, denF_2] = tfdata(F_exp_2000_d, 'v');
F_exp_4000_d = F_exp_4000*z^(-FFdelay);
[numF_4, denF_4] = tfdata(F_exp_4000_d, 'v');

iterations = 10;
u0 = zeros(length(r),1);
y = zeros(length(r),1);

figure();
plot(t, r, '-','linewidth',1.5); hold on;
y0 = lsim(sys_est, r, t);
plot(t, y0, '-.','linewidth',1.5);
for i = 1:iterations
    e = r-y;
    % a = lsim(L, e, t);
    a = filter(numF_4, denF_4, e);
    % a = conv(e, numL/sum(numL), 'same');
    u1 = u0 + [a(FFdelay+1:end); zeros(FFdelay,1)];
    % u1 = u0 + a;
    u1 = conv(u1, numQ/4, 'same');

    y = lsim(sys_est, u1, t);
    u0 = u1;
    if i == 1 || i == 7
        plot(t, y, '-.','linewidth',1.5);
    end
end
xlabel('Times[s]','interpreter','latex','fontsize',tickfont);
ylabel('Magnitude[m]','interpreter','latex','fontsize',tickfont);
title('Model-free ILC tracking', 'fontsize',titlefont);
legend("r","iter: 0","iter: 1","iter: 7",'interpreter','latex','fontsize',legendfont,'location','best')
% legend("r","$y_{exp}:2k$","$y_{exp}:4k$",'interpreter','latex','fontsize',legendfont,'location','best')

grid on; grid minor;
xlim([0.9995 1.0025])
ylim([-1.1 1.1])

%% Model-free RC
% experiment
filename = './experiments/0610/Group5-2024-06-10-17-47-29-RCFF-(FF-ILCFF-TR-delta-0.02s-0.20deg-beta1.0-alpha0.50-delay47-sim-2000-normalized-F).csv';
data_sim_2000 = readmatrix(filename);
filename = './experiments/0610/Group5-2024-06-10-11-25-17-RCFF-(FF-ILCFF-TR-delta-0.02s-0.20deg-beta1.0-alpha0.50-delay47-exp-2000-normalized-F).csv';
data_exp_2000 = readmatrix(filename);
filename = './experiments/0610/Group5-2024-06-10-17-48-14-RCFF-(FF-ILCFF-TR-delta-0.02s-0.20deg-beta1.0-alpha0.50-delay47-sim-4000-normalized-F).csv';
data_sim_4000 = readmatrix(filename);
filename = './experiments/0610/Group5-2024-06-10-11-27-07-RCFF-(FF-ILCFF-TR-delta-0.02s-0.20deg-beta1.0-alpha0.50-delay47-exp-4000-normalized-F).csv';
data_exp_4000 = readmatrix(filename);
filename = './tri-400Hz-2s-padding1.0s.csv';
r = readmatrix(filename);
t = (0:Ts:(Duration + padding*2))'; 

y_sim_2000 = data_sim_2000(:,3);
y_exp_2000 = data_exp_2000(:,3);
y_sim_4000 = data_sim_4000(:,3);
y_exp_4000 = data_exp_4000(:,3);
% shift FFdelay?
t = t(1:end-FFdelay);
r = r(1:end-FFdelay);
y_sim_2000 = y_sim_2000(FFdelay+1:end);
y_exp_2000 = y_exp_2000(FFdelay+1:end);
y_sim_4000 = y_sim_4000(FFdelay+1:end);
y_exp_4000 = y_exp_4000(FFdelay+1:end);

figure();
plot(t, r, '-','linewidth',1.5); hold on;
plot(t, y_exp_2000,t, y_exp_4000, '-','linewidth',1.5);
xlabel('Times[s]','interpreter','latex','fontsize',tickfont);
ylabel('Magnitude[m]','interpreter','latex','fontsize',tickfont);
title('Model-free RC tracking', 'fontsize',titlefont);
legend("r","$y_{exp}:2k$","$y_{exp}:4k$",'interpreter','latex','fontsize',legendfont,'location','best')
grid on; grid minor;
xlim([1.0025 1.01]);
ylim([-1.1 1.1])

% simulation
Np = fs/f;
N1 = FFdelay;
N2 = 1;
filters = [F_exp_2000, F_exp_4000];
figure();
plot(t, r, '-','linewidth',1.5); hold on;
for F = filters
    Cr = feedback(z^(-Np+N1+N2)*Q*z^(-N2),-z^(-N1))*F*z^(-N1);
    sys_closed_RC = feedback(Cr*sys_est,1);
    y = lsim(sys_closed_RC, r, t);
    plot(t, y , '-','linewidth',1.5);
end
xlabel('Times[s]','interpreter','latex','fontsize',tickfont);
ylabel('Magnitude[m]','interpreter','latex','fontsize',tickfont);
title('Model-free RC tracking', 'fontsize',titlefont);
legend("r","$y_{exp}:2k$","$y_{exp}:4k$",'interpreter','latex','fontsize',legendfont,'location','best')
grid on; grid minor;
xlim([1.0025 1.01]);
ylim([-1.1 1.1]);