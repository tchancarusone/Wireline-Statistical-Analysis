
function [c,dfe,cheq] = dfe_lineq(N,Ndfe,fpulsesampled,SNR,eqdelay)
% function [c,dfe,cheq] = dfe_lineq(N,Ndfe,fpulsesampled,SNR,eqdelay)
% 
% c = baud-spaced equalizer taps
% dfe = DFE taps
% cheq = equalized baud-rate channel response
%
% N = number of fractionally-spaced equalizer taps (must be even)
% Ndfe = number of (baud-spaced) DFE taps
% fpulsesampled = sampled channel pulse response (length must be even)
% SNR in dB
% eqdelay = approximate desired delay through equalizer (in # of taps)

n = length(fpulsesampled); % length of channel response 

% desired response
Eh = zeros(1,n+N-1);
delay = round(find(fpulsesampled==max(fpulsesampled),1) + eqdelay);
Eh(delay) = 1;

% make channel convolution matrix
clear F
for k = 1:N
    F(k,:) = [zeros(1,k-1) fpulsesampled zeros(1,N-k)];
end
% ignore columns that will be taken care of by the DFE
F(:,delay+1:delay+Ndfe) = 0;

% optimal equalizer
A = F*F' + 10^-(SNR/10)*eye(N);
c = Eh * F' * inv(A);

cheq = conv(fpulsesampled,c);
dfe = cheq(delay+1:delay+Ndfe);
cheq(delay+1:delay+Ndfe) = 0;
