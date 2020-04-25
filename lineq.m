
function [c,cheq] = lineq(N,fpulsesampled,SNR,target,eqdelay)
% function [c,cheq] = lineq(N,fpulsesampled,SNR,target,eqdelay)
%
% c = baud-spaced equalizer taps
% cheq = equalized baud-rate channel response
%
% N = number of fractionally-spaced equalizer taps 
% fpulsesampled = sampled channel pulse response
% SNR in dB
% target = partial response equalizer target (1 for no partial response)
% eqdelay = approximate desired delay through equalizer (in # of taps)

n = length(fpulsesampled); 

% make channel convolution matrix
clear F
for k = 1:N
    F(k,:) = [zeros(1,k-1) fpulsesampled zeros(1,N-k)];
end

% desired response
Eh = zeros(1,n+N-1);
delay = round(find(fpulsesampled==max(fpulsesampled),1) + eqdelay);
Eh(delay:delay+length(target)-1) = target;

% optimal equalizer
A = F*F' + 10^-(SNR/10)*eye(N);
c = Eh * F' * inv(A);

cheq = c*F;
