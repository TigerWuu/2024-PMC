clc;close all;clear;

ref_amp = 30;
ref_freq = 50;
ref_delay = pi/4;


fs = 1000;
t = 0:1/fs:1-1/fs;
x = ref_amp*cos(2*pi*ref_freq*t - ref_delay)


y = fft(x);
z = fftshift(y);

ly = length(y);
f = (-ly/2:ly/2-1)/ly*fs;



figure;
subplot(2,1,1);
stem(f,abs(z)/(fs/2))
xlabel 'Frequency (Hz)'
ylabel 'Magnitude'
xlim([0 fs/2]);
grid


subplot(2,1,2);
tol = 1e-6;
z(abs(z) < tol) = 0;
theta = angle(z);
stem(f,theta/pi)
xlabel 'Frequency (Hz)'
ylabel 'Phase / \pi'
xlim([0 fs/2]);
grid














