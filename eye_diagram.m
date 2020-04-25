function eye_diagram(y,OSR,varargin)
% function eye_diagram(y,OSR,Tbaud,gaussian_noisepower)
%
% Plots an eye diagram for the waveform stored in variable y.
% OSR = the number of points in y per baud interval (e.g. y-vector
% is 10 Gbaud/sec data sampled at 5 ps => OSR = 20)
% Tbaud is an optional argument for the baud rate (default is 1)
% If gaussian_noisepower is not specified, it defaults to zero

if nargin == 3
    Tbaud = varargin{1};
elseif nargin == 4
    Tbaud = varargin{1};
    noisepower = varargin{2};
    y = y + sqrt(noisepower)*randn(size(y));
else
    Tbaud = 1;
end
t = linspace(0,2*Tbaud,2*OSR+1);

yplot = reshape(y(1:2*OSR*floor(length(y)/(2*OSR))),2*OSR,floor(length(y)/(2*OSR)));
plot(t(1:end-1),yplot,'k-');
grid on