Tbaud = 1/10e9;
OSR = round(Tbaud/Ts); % ensure OSR is even
Tbaud = OSR*Ts;
fbaud = 1/Tbaud;
jittrms = 0.01*Tbaud;
rxBW = 0.8; txBW = 0.8;
[RXnum,RXden] = butter(1,2*rxBW/OSR);
[TXnum,TXden] = butter(1,2*txBW/OSR);
delaynum = [zeros(1,floor(OSR/2)) 1];
delayden = 1;
ic = zeros(1,10); ic(2) = 1;
%ic = [0.8310 0  -1.8837 0   3.8833  0 -1.6524   0 0.2004 0];
mu1  = .3e-3*[1 1 1 1 1 1 1 1 1 1];
mu2  = .3e-3*[1 0 1 1 1 1 1 1 1 1];
a = 0e-5; % lms integrator loss (0 = ideal LMS)
signdata = 2; % 1 = full precision, 2 = sign-data
signerror = 2; % 1 = full precision, 2 = sign-data
hstep1 = filter(TXnum,TXden,hstep);
hstep2 = filter(RXnum,RXden,hstep1);
hpulse = hstep2(OSR+1:end) - hstep2(1:end-OSR);
hpulse = hpulse(1:end-mod(length(hpulse),OSR));

% determine the noise power for a given SNR
SNR = 30; % dB
T = round(OSR/2); % tap spacing
phase = mod(find(hpulse==max(hpulse),1)-1,T) + 1;
fpulsesampled = hpulse(phase:T:end); % sample the pulse response
nrms = sqrt(sum(fpulsesampled.^2)/2*10^(-SNR/10));
noise_num = gausswin(round(OSR/3))';
noise_num = noise_num./sqrt(sum(noise_num.^2));

% loop filter design
Kvco = 200e6; % (Hz/V)
% bang-bang loop
alpha = 1e-8; % integral path
stab_factor = 1e3;
beta = 2*alpha/(stab_factor*Tbaud); % proportional (bang-bang) path 
LFnumd = [(alpha+beta) -beta];
LFdend = [1 -1];
init_ph = -1.7; % rx vco initial phase (in radians)