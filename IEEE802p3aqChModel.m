Tbaud = 1/10e9;
OSR = 20; % ensure OSR is even
Ts = Tbaud/OSR;

t = 0:Ts:40*Tbaud;

% defining the gaussian filter
trf2080 = 47e-12;
Tgauss = trf2080/(2*erfinv(0.8-0.2));
hgauss = exp(-(t-mean(t)).^2/Tgauss^2) / (Tgauss*sqrt(pi));

% 4th-order 7.5-GHz Bessell filter
[Bc,Ac] = besself(4,2*pi*7.5e9/.655);
[Bd,Ad] = bilinear(Bc,Ac,1/Ts);
hBT = impulse(tf(Bd,Ad,Ts),t);

% ISI generator
% AISI = [0.254 0.453 0.155 0.138]; % mostly post-cursor 
% AISI = [0.158 0.176 0.499 0.167]; % mostly pre-cursor 
AISI = [0.000 0.513 0.000 0.487]; % split-pulse
hISI = upfirdn(AISI,1,round(0.75*OSR),1);

himp = Ts*conv(hgauss,conv(hISI,hBT));
hpulse = conv(ones(1,OSR),himp)';
