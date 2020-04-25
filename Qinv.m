function y = Qinv(x)
% y = Qinv(x)
%
% y chosen to satisfy
%     x = Qfun(y)

y = sqrt(2)*erfcinv(2*x);