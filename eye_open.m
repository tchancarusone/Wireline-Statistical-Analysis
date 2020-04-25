function [eo,ISIpen,kp,prffe,k] = eye_open(ch,ffe,delay,OSR)
% function [eo,ISIpen,kp,prffe,k] = eye_open(ch,ffe,delay,OSR)
%
% ch = baud-rate pulse response
% ffe = linear equalizer tap weights
% delay = equalizer tap spacing (in terms of samples)
% OSR = ratio of baud rate to simulation sampling rate
%
% eo = worst case eye opening vector (length OSR)

Nlevels = 2;
prffe = conv(ch,upfirdn(ffe,1,delay,1));
prffe = prffe(1:length(ch));

for k = 1:OSR
   
   % obtain baud-rate channel pulse response
   h = prffe(k:OSR:end);
   
   peak = find(h==max(h)); peak = peak(1);
   sig = h(peak)/(Nlevels-1);
   isi = sum(abs([h(1:peak-1) h(peak+1:end)]));
   eo(k) = sig-isi;
end
k = find(eo==max(eo)); k = k(1);
h = prffe(k:OSR:end);
peak = find(h==max(h)); peak = peak(1);
A0 = sum(h);
ISIpen = eo(k)/A0;
kp = -sign(h); kp(peak) = 1; kp = kp(end:-1:1);
