clc;
clear;
close all;

%% initialization
load("sys_est.mat");
load("Pa.mat");
% sys_est = Pa;

fs = 20000;
Ts = 1/fs;
hz = 400;

%% open-loop system
t= 0:Ts:1;
r = 0.08*sin(2*pi*hz*t);
lsim(sys_est, r, t);

%% controller design : k controller
% spec ?
rlocus(sys_est);
zgrid;

%% closed-loop : P controller
kp = 1.5;
sys_est_cls_P = feedback(sys_est*kp, 1);
lsim(sys_est_cls_P, r, t);

%% controller design : k controller
% spec ?
z = tf('z',Ts);
sys_est_L = ((z-1)*sys_est)/(z*(1+kp*sys_est));
rlocus(sys_est_L);
zgrid;

%% closed-loop : PD controller
kd = 10;
sys_est_cls_PD = feedback(sys_est_L*kd, 1);
lsim(sys_est_cls_PD, r, t);