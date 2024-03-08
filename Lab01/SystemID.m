clc;
clear;

%% initialization
fs = 20000;
Ts = 1/fs;
w = {10, 10000};
load("Pa.mat");

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
plotoptions = bodeoptions;
plotoptions.Grid = 'on';
plotoptions.FreqUnits = 'Hz';
plotoptions.Title.String = 'Bode Plot of Transfer Function';
% plotoptions.PhaseWrapping = "on";
% plotoptions.PhaseWrappingBranch = -360;

bodeplot(sys_hw, w, plotoptions);
%% complex curve fitting by matlab toolbox
figure();
sys_frd=frd(Gks,ufreq);
model=idtf([NaN NaN NaN NaN NaN],[1 NaN NaN NaN NaN NaN]);
model.Structure.den.Minimum=0;
sys_est = tfest(sys_frd,model);
sys_est = c2d(sys_est,Ts); % why??

plotoptions = bodeoptions;
plotoptions.Grid = 'on';
plotoptions.FreqUnits = 'Hz';
plotoptions.Title.String = 'Bode Plot of Transfer Function';
% plotoptions.PhaseWrapping = "on";
% plotoptions.PhaseWrappingBranch = -360;

bodeplot(sys_est, w, plotoptions);
save("sys_est.mat", "sys_est");
%% true tf
figure();
numerator = [-0.001647,0.05631,0,0,0];
denominator = [1,-2.804,3.599,-2.594,1.042,-0.1882];

sys_true_d = tf(numerator,denominator,Ts);
sys_true_c = d2c(sys_true_d);

plotoptions = bodeoptions;
plotoptions.Grid = 'on';
% plotoptions.FreqUnits = 'Hz';
plotoptions.Title.String = 'Bode Plot of Transfer Function';
% plotoptions.PhaseWrapping = "on";
% plotoptions.PhaseWrappingBranch = -360;

bodeplot(sys_true_d, w, plotoptions);

%% compare
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
bodeplot(sys_est, w, plotoptions);
legend;

% [Gk , f] = frf_spa(y,u,Ts,"H1",true);
% plot(t,y,t,u);
% plot(u)