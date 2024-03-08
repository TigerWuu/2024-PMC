


function [Gk, f] = frf_spa(y,u,Ts,algorithm,plotOption)
Fs = 1/Ts; % Fs is sampling frequency
M = length(u);
% why is using it better??
N = M/10;
% apply window funciton ? 
% w = blackmanharris(N);
w=1;

% if we don't use??
nfft = 2^nextpow2(N);

% calculate u, y, uy, yu correlation
% why biased ???
Ru = xcorr(u,"biased");
Ry = xcorr(y,"biased");
Ruy = xcorr(u,y,"biased");
Ryu = xcorr(y,u,"biased");

Ru = Ru(N:2*N-1);
Ry = Ry(N:2*N-1);
Ruy = Ruy(N:2*N-1);
Ryu = Ryu(N:2*N-1);

% calculate fourier transform of the correlation of u, y, uy, yu

fu = fft(Ru.*w, nfft);
fy = fft(Ry.*w, nfft);
fuy = fft(Ruy.*w, nfft);
fyu = fft(Ryu.*w, nfft);


% calculate the estimated FRF 
if algorithm == "H1"
   Gk =  fyu./fu;
   Gk = Gk(1:(nfft/2+1)); % w of dft represents normalized frequency
   % Gk = Gk(N:2*N-1);
else
   Gk =  fy./fuy;
   Gk = Gk(1:(nfft/2+1));
   % Gk = Gk(N:2*N-1);
end

% calculate the f ??? relation???
% Gk(w) : w represents normalized frequency, real frequency = sampling frequency * w
% f = Fs*(0:(N/2))/N;
f = Fs*(0:(nfft/2)) / nfft;

% plot the bode plot
if plotOption == true
    subplot(211)
    plot(f, mag2db(abs(Gk)),'-','linewidth',1.5); hold on;
    % semilogx(f,mag2db(abs(Gk)),'-','linewidth',1.5); hold on;
    % semilogx(f,mag2db(abs(H_true)),'-','linewidth',2); hold on;
    grid on; grid minor;
    ylabel('Magnitude [dB]','interpreter','latex','fontsize',12)
    
    subplot(212)
    plot(f,rad2deg(angle(Gk)),'-','linewidth',1.5); hold on;
    % semilogx(f,rad2deg(angle(Gk)),'-','linewidth',1.5); hold on;
    % semilogx(f,rad2deg(angle(H_true)),'-','linewidth',2); hold on;
    grid on; grid minor;
    ylabel('Phase [deg]','interpreter','latex','fontsize',12)
    xlabel('Frequency [Hz]','interpreter','latex','fontsize',12)
    legend(algorithm,'interpreter','latex','fontsize',12,'location','best')
end

end