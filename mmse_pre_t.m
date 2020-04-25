
% tap spacing
T = OSR;

% partial respose target
target = [1 0 0];

SNR = 30; % dB

% matched filter
% wmf = hpulse(end:-1:1)./10^-(SNR/10); fpulse = conv(hpulse,wmf); % include matched filter  (assuming white noise)
fpulse = hpulse; % skip matched filter

% equivalent discrete-time channel, F = H H* / N
phase = mod(find(fpulse==max(fpulse),1)-1,T) + 1;
fpulsesampled = fpulse(phase:T:end); % sample the pulse response

% equalizer length
n = length(fpulsesampled); % length of channel response
N = n + 100;
N = 5;

% make channel convolution matrix
clear F
for k = 1:N
    F(k,:) = [zeros(1,k-1) fpulsesampled zeros(1,N-k)];
end

% desired response; centre equalizer tap is the "main" one
Eh = zeros(1,n+N-1);
delay = round(find(fpulsesampled==max(fpulsesampled),1) + N/2) - 1;
Eh(delay:delay+length(target)-1) = target;

% optimal equalizer
A = F*F' + 10^-(SNR/10)*eye(N);
c = Eh * F' * inv(A)

c2 = upfirdn(c,1,T,1);
ch = filter(TXnum,TXden,filter(RXnum,RXden,conv(c2,h)));
chpulse = conv(ch,ones(1,OSR));

nrms = sqrt(sum(fpulsesampled.^2)/2*10^(-SNR/10));
n = conv(nrms*randn(1,10000),sqrt([.2 .2 .2 .2 .2]));
n2 = conv(n,upfirdn(c,1,T,1)); n2 = n2(T*N/2:T*N/2+1e4-1)';
