function [prb]= lsrprb(N,xstart,varargin)
% Generates PRB-N using the LINEAR SHIFT REGISTER
% function [prb]= lsrprb(N,xstart,length)
%
%   N = length of linear feedback shift register (LSR)
%   xstart = vector of initial binary values in length-N LSR
%   length = length of sequence to be generated (default is 2^N-1)

if nargin == 3
    len = varargin{1};
else
    len = 2^N-1;
end

if N==7
   indx=[7 6];
elseif N==9
   indx=[9 5];
elseif N==10
   indx=[10 7];
elseif N==11
   indx=[11 9];   
elseif N==13
   indx=[13 4 3 1];
elseif N==15
   indx=[15 14];
elseif N==20
   indx=[20 17];
elseif N==23
   indx=[23 18];
elseif N==31
   indx=[31 28];   
end

x= xstart;
for i=1:len
   prb(i)= x(N);
   newbit= mod(sum(x(indx)),2);
   x= [newbit x(1:N-1)];
end 
