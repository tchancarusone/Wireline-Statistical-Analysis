function Q = Qfun(x)
% function Q = Qfun(x)
%
% Q = [ 1 / sqrt(2*pi) ] * integral from x to inf exp(y^2 / 2)

Q = 0.5*erfc(x/sqrt(2));