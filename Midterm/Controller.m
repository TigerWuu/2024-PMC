clc;
clear;
close all;
%%
plotoptions = bodeoptions;
plotoptions.Grid = 'on';
plotoptions.FreqUnits = 'Hz';   
plotoptions.Title.String = 'Sensitivity';
plotoptions.Title.FontSize = 15;
plotoptions.Title.FontWeight = "bold";
plotoptions.Xlabel.FontSize= 12;
plotoptions.Ylabel.FontSize= 12;
plotoptions.Ticklabel.FontSize= 10;

%% initialization
load("sys_est.mat");
% load("Pa.mat");
% sys_est = Pa;

fs = 20000;
Ts = 1/fs;

%% tf2ss
b = sys_est.Numerator;
a = sys_est.Denominator;
[A,B,C,D] = tf2ss(b,a);

sys_est_ss = ss(A,B,C,D,Ts);
% step(sys_est_ss);
% pzmap(sys_est_ss)
% zgrid;
%% design L
% poleplacement
% poles = [0.9, 0.91, 0.92, 0.93, 0.94]; % bad
% poles = [0.5, 0.6, 0.7, 0.8, 0.9];
% poles = [0.5, 0.51, 0.52, 0.53, 0.54]; % how to choose
poles = [0.5, 0.55, 0.62, 0.63, 0.6]; % how to choose
% poles = [0.1, 0.2, 0.3, 0.4, 0.5];
% poles = [0.3, 0.31, 0.32, 0.33, 0.34];
% poles = [0.1, 0.11, 0.12, 0.13, 0.14];
% L = [0;0;0;0;0];
L = place(A',C',poles)';


eig(A-L*C)
%% design K (LQR)
Q = 1*eye(5);
R = 1;
N = zeros(5,1);
[K,S,e] = dlqr(A,B,Q,R,N);

eig(A-B*K)
%% design K (poleplacement)
% poles = [0.5, 0.6, 0.7, 0.8, 0.9];
% poles = [0.9, 0.91, 0.92, 0.93, 0.94];
poles = [0.85, 0.8, 0.1, 0.11, 0.15];
% poles = [0.85, 0.8, 0.7, 0.65, 0.6];
K = place(A,B,poles);
% K = [0 ,0, 0, 0, 0];
% writematrix(K,'./design/K_pp.csv');

eig(A-B*K)
%% design M 
w0 = 400/fs*2*pi; % to normalized frequency(but why?
Ma = [0 , 1 ; -1, 2*cos(w0)];
Mb = [0;1];
Mc = [1,0];
Md = 0;
M = ss(Ma, Mb, Mc, Md, Ts);
[a,b] = ss2tf(Ma, Mb, Mc, Md);
M_tf = tf(a,b,Ts);


A_aug = [A , zeros(length(A),length(Ma));Mb*C, Ma];
B_aug = [B;zeros(length(Ma),1)];
C_aug = [C , zeros(1, length(Ma))];

% poles = [0.85, 0.8, 0.1, 0.11, 0.15, -0.45, -0.4]; % aggressive
poles = [0.85, 0.8, 0.7, 0.65, 0.6, -0.6, -0.55]; % non-aggressive
Kaug = place(A_aug , B_aug , poles);
K = Kaug(1:length(A));
Km = Kaug(length(A)+1:end);

%% Loop gain analysis (with internal model)
Lg = ss(A-L*C-B*K,L,K,0,Ts)*sys_est_ss+ ... 
    ss(A-L*C-B*K,L,K,0,Ts)*sys_est_ss*ss(Ma,Mb,-Km,0,Ts)- ...
    sys_est_ss*ss(Ma,Mb,-Km,0,Ts);

% bode(Lg);
% sensitivity analysis
S = feedback(1,Lg);
hold on;
bodeplot(S, plotoptions);
%% Loop gain analysis (without internal model)
Lg = ss(A-L*C-B*K,L,K,0,Ts)*sys_est_ss;


% bode(Lg);
% sensitivity analysis
S = feedback(1,Lg);
bodeplot(S, plotoptions);
%% save designed parameters
% writematrix(A,'./design/A.csv');
% writematrix(B,'./design/B.csv');
% writematrix(C,'./design/C.csv');
% 
% % writematrix(L,'./design/L.csv');
% % writematrix(K,'./design/K.csv');
% % 
% % writematrix(Ma,'./design/Ma.csv');
% % writematrix(Mb,'./design/Mb.csv');
% % writematrix(-Km,'./design/Mc.csv');
% 
% writematrix(L,'./design/La.csv');
% writematrix(K,'./design/Ka.csv');
% 
% writematrix(Ma,'./design/Maa.csv');
% writematrix(Mb,'./design/Mba.csv');
% writematrix(-Km,'./design/Mca.csv');



