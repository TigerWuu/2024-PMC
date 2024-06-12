function [F,delay, gamma] = zpetc(sysG)
%--------------------------------------------------------------------------
% sysG is stable plant G, F is the zerophase error compensator scaled for unity d.c. gain, delay is the
% non-causal steps of F, gamma is the scaling gain required for making 0 <=GF <=1
% nh: order of zero-phase low pass filter
%--------------------------------------------------------------------------
% sysG = Gyrz;
[B,A,Ts] = tfdata(sysG,'v');
[z,p,k] = tf2zp(B,A);
m = length(z);
n = length(p);
d = n-m;

% seperate stable and unstable zeros
Bplus = [1];    %LHP zero
Bminus = [1];   %RHP zero
nu = 0;         %number of unstable zeros
Zcritic = 1.;   %unit circle for stability boundary
for i = 1:length(z)
    if abs(z(i)) > Zcritic
        Bminus = conv(Bminus, [1 -z(i)]);
        nu = nu+1;
    else 
        Bplus = conv(Bplus, [1 -z(i)]);
    end
end 

sigma = sum(Bminus);        %dc gain
ZnuBmi = fliplr(Bminus);    %conjugate; z^-1=z

b = sigma^2*k;              %dc gain of ZhuBmi=sigma, of A/Bplus = sigma/k = sigma / dcgain(sysG)

[mag,pha]=bode(tf(Bminus,[1,zeros(size(Bminus))],1));
gamma = sigma^2/max(mag)^2;


% construct ZPETC controllers
numZP = conv(A,ZnuBmi)/b;
denZP = conv(real(Bplus),[1 zeros(1,d+2*nu)]);   %d+2*nu: make F to be non-causal, relative order=0;
F = tf(numZP,denZP,Ts);

delay = nu+d;