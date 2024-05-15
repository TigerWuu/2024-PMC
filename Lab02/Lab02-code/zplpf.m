function [HH] = zplpf(N,Wp,FIR_tap,Max_step,Ts)
% Zero-phase low-pass filter
% using N order butterworth
% cut off frequency(Hz): Wp
% truncate at FIR_tap steps


[b,a]=butter(N,Wp*Ts,'low');
IIR_filter=tf(b,a,Ts);

h=impulse(IIR_filter,Max_step);
H=h(1:FIR_tap);
H=H';   %convert to row vector
H=H/sum(H);  % make d.c. gain =1
FIR_H=tf(H,[1,zeros(1,FIR_tap-1)],Ts);
FIR_Hconj=tf(fliplr(H),[1,zeros(1,FIR_tap-1)],Ts);
   
HH=FIR_H*FIR_Hconj;  % MM=H(z^-1)*H(z)