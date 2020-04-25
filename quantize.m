function q = quantize(x,bot,top,Nlevels)
% function q = quantize(x,min,max,Nlevels)
%
% x = quantizer input
% min = lower saturation limit of quantizer
% max = upper saturation limit of quantizer
% Nlevels = number of equally-spaced quantization levels

outputs = linspace(bot,top,Nlevels);

e = outputs - x;
[y,z] = min(abs(e));
q = outputs(z);
