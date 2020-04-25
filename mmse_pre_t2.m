
% tap spacing
T = OSR/2;

% partial respose target
target = [1];

SNR = 30; % dB

% matched filter
% wmf = hpulse(end:-1:1)./10^-(SNR/10); fpulse = conv(hpulse,wmf); % include matched filter  (assuming white noise)
fpulse = hpulse; % skip matched filter

% equivalent discrete-time channel, F = H H* / N
phase = mod(find(fpulse==max(fpulse),1)-1,T) + 1;
fpulsesampled = fpulse(phase:T:end); % sample the pulse response
n = length(fpulsesampled); % length of channel response
if mod(n,2)
    n = n - 1;
    fpulsesampled = fpulsesampled(1:n);
end

% equalizer length (must be even)
% N = n + 100;
N = 12;

% make channel convolution matrix
clear F
for k = 1:N/2
    F(2*k-1,:) = [zeros(1,k-1) fpulsesampled(2:2:end) zeros(1,N/2-k)];
    F(2*k,:) = [zeros(1,k-1) fpulsesampled(1:2:end) zeros(1,N/2-k)];
end

% desired response; approx centre equalizer tap is the "main" one
Eh = zeros(1,(n+N)/2-1);
delay = round(find(fpulsesampled==max(fpulsesampled),1)/2 + N/4);
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
