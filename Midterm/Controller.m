clc;
clear;
close all;

%% initialization
load("sys_est.mat");
fs = 20000;
Ts = 1/fs;

%% tf2ss
b = sys_est.Numerator;
a = sys_est.Denominator;
[A,B,C,D] = tf2ss(b,a);

sys_est_ss = ss(A,B,C,D,Ts);

%% design L
% poleplacement
L = place(A',C',[0.5, 0.6, 0.7, 0.8, 0.9])';

%% design K
Q = eye(5);
R = 1;
N = zeros(5,1);
[K,S,e] = dlqr(A,B,Q,R,N);