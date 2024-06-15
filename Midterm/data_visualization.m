close all;

%%
t = out.tout;

stepRespons = get(out.logsout, 'step').Values.Data;
step_r = stepRespons(:,1);
step_y = stepRespons(:,2);
step_y_hat = stepRespons(:,3);

sineRespons = get(out.logsout, 'sine').Values.Data;
sine_r = sineRespons(:,1);
sine_y = sineRespons(:,2);
sine_y_hat = sineRespons(:,3);

%% plot style
titlefont = 20;
tickfont = 18;
legendfont = 16;

%% compare y_real vs. y_sim
filename = './experiments/0325/Group5-2024-03-25-11-30-34-StfbObsv_only-0.20deg-step-10.0s.csv';
data = readmatrix(filename)';
save("experiment.mat", "data");

t_r = data(1,:);
step_y_r = data(4,:);

figure();
plot(t_r, step_y_r, '-','linewidth',1.5); hold on;
plot(t, step_y_hat, t, step_r,'-','linewidth',1.5); 
xlabel('Times[s]','interpreter','latex','fontsize',tickfont);
ylabel('Magnitude[m]','interpreter','latex','fontsize',tickfont);
title('Step response', 'fontsize',titlefont);
legend("$y_{experiment}$","$\hat{y}$","$r$",'interpreter','latex','fontsize',legendfont,'location','best')
grid on; grid minor;

%% step response (simulation)
figure();
plot(t, step_r, t, step_y, t, step_y_hat,'-','linewidth',1.5); hold on;
xlabel('Times[s]','interpreter','latex','fontsize',tickfont);
ylabel('Magnitude[m]','interpreter','latex','fontsize',tickfont);
title('Step response', 'fontsize',titlefont);
legend("$r$","$y$","$\hat{y}$",'interpreter','latex','fontsize',legendfont,'location','best')
grid on; grid minor;


%% sinewave tracking (simulation)
figure();
plot(t, sine_r, t, sine_y, t, sine_y_hat,'-','linewidth',1.5); hold on;
xlabel('Times[s]','interpreter','latex','fontsize',tickfont);
ylabel('Magnitude[m]','interpreter','latex','fontsize',tickfont);
title('400 Hz Sinewave tracking', 'fontsize',titlefont);
legend("$r$","$y$","$\hat{y}$",'interpreter','latex','fontsize',legendfont,'location','best')
grid on; grid minor;

%% step response (experiment)
filename = './experiments/0520/Group5-2024-05-20-11-11-32-StfbObsv_only-0.20deg-step-10.0s_Kpp.csv';
data = readmatrix(filename);
t = data(:,1);
step_r = data(:,3);
step_y = data(:,4);
step_y_hat = data(:,5);

figure();
plot( t, step_y, t, step_y_hat,t, step_r,'-','linewidth',1.5); hold on;
xlabel('Times[s]','interpreter','latex','fontsize',tickfont);
ylabel('Magnitude[m]','interpreter','latex','fontsize',tickfont);
title('Step response tracking', 'fontsize',titlefont);
legend("$y$","$\hat{y}$","$r$",'interpreter','latex','fontsize',legendfont,'location','best')
grid on; grid minor;

%% sinewave tracking - nonaggressive poles (experiment)
filename = './experiments/0520/Group5-2024-05-20-11-13-57-StfbObsv_withIMP-0.20deg-400HzSine-10.0s_noaggressive.csv';
data = readmatrix(filename);
t = data(:,1);
sine_r = data(:,3);
sine_y = data(:,4);
sine_y_hat = data(:,5);

figure();
plot(t, sine_r, t, sine_y, t, sine_y_hat,'-','linewidth',1.5); hold on;
xlabel('Times[s]','interpreter','latex','fontsize',tickfont);
ylabel('Magnitude[m]','interpreter','latex','fontsize',tickfont);
title('400Hz sinewave tracking', 'fontsize',titlefont);
legend("$r$","$y$","$\hat{y}$",'interpreter','latex','fontsize',legendfont,'location','best')
grid on; grid minor;

%% sinewave tracking - aggressive poles (experiment)
filename = './experiments/0521/Group5-2024-05-21-10-09-05-StfbObsv_withIMP-0.20deg-400HzSine-10.0s-aggressive.csv';
data = readmatrix(filename);
t = data(:,1);
sine_r = data(:,3);
sine_y_f = data(:,4);
sine_y_hat_f = data(:,5);

figure();
plot(t, sine_r, t, sine_y_f, t, sine_y_hat_f,'-','linewidth',1.5); hold on;
xlabel('Times[s]','interpreter','latex','fontsize',tickfont);
ylabel('Magnitude[m]','interpreter','latex','fontsize',tickfont);
title('400Hz sinewave tracking', 'fontsize',titlefont);
legend("$r$","$y$","$\hat{y}$",'interpreter','latex','fontsize',legendfont,'location','best')
grid on; grid minor;