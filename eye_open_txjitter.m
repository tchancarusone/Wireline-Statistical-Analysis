function [eo,ISIpen,kp,prffe,k] = eye_open_txjitter(ch,ffe,delay,OSR,txjitter)

prffe = conv(ch,upfirdn(ffe,1,delay,1));
prffe = prffe(1:length(ch));
sffe = zeros(1,length(ch));
for k = 1:floor(length(ch)/OSR)
    sffe = sffe + [zeros(1,(k-1)*OSR) prffe(1:length(ch)-(k-1)*OSR)];
end
hffe = sffe - [0 sffe(1:end-1)];

for k = 1:OSR
   
   % obtain baud-rate channel pulse response
   s = sffe(k:OSR:end);
   pr = prffe(k:OSR:end);
   h = hffe(k:OSR:end);
   
   peak = find(pr==max(pr)); peak = peak(1);
   kp = -sign(pr); kp(peak) = 1; kp = kp(end:-1:1);
   ahat = kp - [0 kp(1:end-1)];
   eo(k) = ahat(end:-1:1)*s' - abs(ahat(end:-1:1))*txjitter*h';
end
k = find(eo==max(eo)); k = k(1);
h = prffe(k:OSR:end);
peak = find(h==max(h)); peak = peak(1);
A0 = sum(h);
ISIpen = eo(k)/A0;
kp = -sign(h); kp(peak) = 1; kp = kp(end:-1:1);
