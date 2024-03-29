clc;
clear;
close all;
%% initialization
fs = 20000;
Ts = 1/fs;
w = {10*2*pi, 10000*2*pi}; % hz to rad/s because the unit of w is rad/s
load("Pa.mat");
plotoptions = bodeoptions;
plotoptions.Grid = 'on';
plotoptions.FreqUnits = 'Hz';   
plotoptions.Title.String = 'Bode Plot (Complex Curve Fitting)';
plotoptions.Title.FontSize = 15;
plotoptions.Title.FontWeight = "bold";
plotoptions.Xlabel.FontSize= 12;
plotoptions.Ylabel.FontSize= 12;
plotoptions.Ticklabel.FontSize= 10;

%% concatenate all the u and y
ufreq = [10, 20, 30, 50, 100, 200, 300, 500, 1000, 2000, 3000, 5000];
u = [];
y = [];
Gks= [];
for i=ufreq
    filename = strcat(strcat('./Group11-2024-03-04-11-10-51-sine-sweep-0.08deg/',num2str(i)),'Hz.csv');
    data = readmatrix(filename);
    datalen = length(data(:,1));
    t = data(:,1);
    u = fft(data(:,2));
    y = fft(data(:,3));
    ub = max(u(1:datalen/2));
    yb = max(y(1:datalen/2));
    Gk = yb./ub;
    % f = fs/2*(0:length(y)-1)/length(y);
    Gks = cat(2, Gks, Gk);
    % u = cat(1,u,data(:,2));
    % y = cat(1,y,data(:,3));
    % scatter(i,angle(ub))
    % plot(t,data(:,2),t,data(:,3))
    % plot(angle(y))
    % pause
end

%% bodeplot
figure();
% magnitude
subplot(211)
% plot(ufreq, mag2db(abs(Gks)),'-','linewidth',1.5); hold on;
semilogx(ufreq,mag2db(abs(Gks)),'-','linewidth',1.5); hold on;
grid on; grid minor;
ylabel('Magnitude [dB]','interpreter','latex','fontsize',12)
% phase
subplot(212)
% plot(ufreq,rad2deg(angle(Gks)),'-','linewidth',1.5); hold on;
semilogx(ufreq,rad2deg(angle(Gks)),'-','linewidth',1.5); hold on;
grid on; grid minor;
ylabel('Phase [deg]','interpreter','latex','fontsize',12)
xlabel('Frequency [Hz]','interpreter','latex','fontsize',12)
legend("frd",'interpreter','latex','fontsize',12,'location','best')

%% complex curve fitting by System Indentification HW3
figure();
n = 5;
[a,b] = complexfit(Gks(1:12)', n, ufreq, fs);
a = cat(1,[1],a);
sys_hw = tf(b',a',Ts);

plotoptions.Title.String = 'Bode Plot (Complex Curve Fitting)';
% plotoptions.PhaseWrapping = "on";
% plotoptions.PhaseWrappingBranch = -360;
bodeplot(sys_hw, w, plotoptions);
save("sys_hw.mat", "sys_hw");
% [mag,phase,wout] = bode(sys_hw, w);
% % magnitude
% subplot(211)
% semilogx(wout,mag2db(abs(squeeze(mag))),'-','linewidth',1.5); hold on;
% grid on; grid minor;
% ylabel('Magnitude [dB]','interpreter','latex','fontsize',12)
% title('Bode Plot (Complex Curve Fitting)','interpreter','latex','fontsize',20);
% % phase
% subplot(212)
% semilogx(wout,squeeze(phase),'-','linewidth',1.5); hold on;
% grid on; grid minor;
% ylabel('Phase [deg]','interpreter','latex','fontsize',12)
% xlabel('Frequency [Hz]','interpreter','latex','fontsize',12)

% legend("frd",'interpreter','latex','fontsize',12,'location','best')
%% complex curve fitting by matlab toolbox
figure();
sys_frd=frd(Gks,ufreq,Ts,'FrequencyUnit','Hz');
model=idtf([NaN NaN 0 0 0],[1 NaN NaN NaN NaN NaN]);
sys_est = tfest(sys_frd,model,'Ts', Ts);
% sys_est = c2d(sys_est,Ts); % why??

plotoptions.Title.String = 'Bode Plot (tfest)';
% plotoptions.PhaseWrapping = "on";
% plotoptions.PhaseWrappingBranch = -360;


[mag,phase,wout] = bode(sys_est, w); % the unit of wout is rad/s
wout = wout/2/pi; % rad/s to Hz
% magnitude
subplot(211)
semilogx(wout,mag2db(abs(squeeze(mag))),'-','linewidth',1.5); hold on;
semilogx(ufreq,mag2db(abs(Gks)),'.','MarkerSize',10);
grid on; grid minor;
ylabel('Magnitude [dB]','interpreter','latex','fontsize',12)
title('Bode Plot (Complex Curve Fitting)','interpreter','latex','fontsize',20);
% phase
subplot(212)
semilogx(wout,squeeze(phase),'-','linewidth',1.5); hold on;
semilogx(ufreq,rad2deg(angle(Gks)),'.','MarkerSize',10);
grid on; grid minor;
ylabel('Phase [deg]','interpreter','latex','fontsize',12)
xlabel('Frequency [Hz]','interpreter','latex','fontsize',12)
% bodeplot(sys_est, w, plotoptions);
save("sys_est.mat", "sys_est");
%% Paper's tf
figure();
numerator = [-0.001647,0.05631,0,0,0];
denominator = [1,-2.804,3.599,-2.594,1.042,-0.1882];

sys_true_d = tf(numerator,denominator,Ts);
sys_true_c = d2c(sys_true_d);

plotoptions.Title.String = 'Bode Plot (Paper)';
% plotoptions.PhaseWrapping = "on";
% plotoptions.PhaseWrappingBranch = -360;

bodeplot(sys_true_d, w, plotoptions);

%% TA's tf
figure();
plotoptions.Title.String = 'Bode Plot (TA)';
bodeplot(Pa, w, plotoptions);

%% compare : frequency response
% true vs. toolbox
figure();
hold on;
bodeplot(sys_true_d, w, plotoptions);
bodeplot(sys_est, w, plotoptions);
legend;

% hw3 vs. toolbox 
figure();
hold on;
bodeplot(sys_hw, w, plotoptions);
bodeplot(sys_est, w, plotoptions);
legend;

% hw3 vs. true
figure();
hold on;
bodeplot(sys_hw, w, plotoptions);
bodeplot(sys_true_d, w, plotoptions);
legend;

% Ta vs. true
figure();
hold on;
bodeplot(Pa, w, plotoptions);
bodeplot(sys_true_d, w, plotoptions);
legend;

% Ta vs.hw3
figure();
plotoptions.Title.String = 'Bode Plot : TA vs. CCF';
hold on;
bodeplot(Pa, w, plotoptions);
bodeplot(sys_hw, w, plotoptions);
legend;


% Ta vs.hw3
figure();
plotoptions.Title.String = 'Bode Plot : TA vs. tfest';
hold on;
bodeplot(Pa, w, plotoptions);
bodeplot(sys_est, w, plotoptions);
legend;
% [Gk , f] = frf_spa(y,u,Ts,"H1",true);
% plot(t,y,t,u);
% plot(u)


%% compare : step response
% Ta vs. tfest
figure();
hold on;
tf = 0.001;
[y,t] = step(Pa,tf);
plot(t, squeeze(y), 'LineWidth',2);
[y,t] = step(sys_est,tf);
plot(t, squeeze(y), 'LineWidth',2)

grid on; grid minor;
ylabel('Amplitude[deg]','interpreter','latex','fontsize',12);
xlabel('Time[s]','interpreter','latex','fontsize',12);
title('Step: Ta vs. tfest','interpreter','latex','FontWeight','bold','fontsize',15);
legend("Ta's","tfest",'interpreter','latex','fontsize',12,'location','best');

% Ta vs. ccf
figure();
hold on;
tf = 6*0.001;
[y,t] = step(Pa,tf);
plot(t, squeeze(y), 'LineWidth',2);
[y,t] = step(sys_hw,tf);
plot(t, squeeze(y), 'LineWidth',2)

grid on; grid minor;
ylabel('Amplitude[deg]','interpreter','latex','fontsize',12);
xlabel('Time[s]','interpreter','latex','fontsize',12);
title('Step: Ta vs. CCF','interpreter','latex','FontWeight','bold','fontsize',15);
legend("Ta's","CCF",'interpreter','latex','fontsize',12,'location','best');