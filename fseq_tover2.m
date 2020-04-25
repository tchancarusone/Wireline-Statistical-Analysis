
function [c,cheq] = fseq_tover2(N,fpulsesampled,SNR,target,eqdelay)
% function [c,cheq] = fseq_tover2(N,fpulsesampled,SNR,target,eqdelay)
%
% c = fractionally spaced equalizer taps
% cheq = equalized baud-rate channel response
%
% N = number of fractionally-spaced equalizer taps (must be even)
% fpulsesampled = sampled channel pulse response (length must be even)
% SNR in dB
% target = partial response equalizer target (1 for no partial response)
% eqdelay = approximate desired delay through equalizer (in # of taps)

n = length(fpulsesampled); % length of channel response must be even
if mod(n,2)
    n = n - 1;
    fpulsesampled = fpulsesampled(1:n);
end

% equalizer length must be even
if mod(N,2)
    N = N + 1;
end

% make channel convolution matrix
clear F
for k = 1:N/2
    F(2*k-1,:) = [zeros(1,k-1) fpulsesampled(2:2:end) zeros(1,N/2-k)];
    F(2*k,:) = [zeros(1,k-1) fpulsesampled(1:2:end) zeros(1,N/2-k)];
end

% desired response
Eh = zeros(1,(n+N)/2-1);
delay = round(find(fpulsesampled==max(fpulsesampled),1)/2 + eqdelay/2);
Eh(delay:delay+length(target)-1) = target;

% optimal equalizer
A = F*F' + 10^-(SNR/10)*eye(N);
c = Eh * F' * inv(A);

cheq = c*F;
