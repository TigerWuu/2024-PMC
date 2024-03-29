clc;
clear;
close all;

%% initialization
load("sys_est.mat");
load("sys_hw.mat");
load("Pa.mat");
% sys_est = Pa;
% sys_est = sys_hw;

fs = 20000;
Ts = 1/fs;
hz = 400;

%% open-loop system
t= 0:Ts:1;
r = 0.08*sin(2*pi*hz*t);
lsim(sys_est, r, t);

%% controller design : P controller
% spec ?
rlocus(sys_est);
zgrid;

%% closed-loop : P controller
kp = 0.9;
sys_est_cls_P = feedback(sys_est*kp, 1);
lsim(sys_est_cls_P, r, t);

%% controller design : PD controller
% spec ?
z = tf('z',Ts);
sys_est_L_pd = ((z-1)*sys_est)/(z*(1+kp*sys_est));
rlocus(sys_est_L_pd);
zgrid;

%% closed-loop : PD controller
kd = 1;
sys_est_cls_PD = feedback(sys_est_L_pd*kd, 1);
lsim(sys_est_cls_PD, r, t);

%% controller design : PID controller
% spec ?
z = tf('z',Ts);
sys_est_L_pid = ((1-z^(-1))^(-1)*sys_est)/(z*(1+kp*sys_est+sys_est*kd*(1-z^(-1))));
rlocus(sys_est_L_pid);
zgrid;

%% closed-loop : PID controller
ki = 0.4;
sys_est_cls_PID = feedback(sys_est_L_pid*ki, 1);
lsim(sys_est_cls_PID, r, t);
%% pid close bode
plotoptions = bodeoptions;
plotoptions.FreqUnits = 'Hz';   
hold on;
bodeplot(sys_est_cls_PD, plotoptions);
bodeplot(sys_est, plotoptions);


