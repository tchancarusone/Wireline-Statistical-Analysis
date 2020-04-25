***
Some New functions for MMSE equalization and statistical eye analysis
Plus the ABCD wireline channel modeling functions
***

mmse_pre

frequency-domain approach [ref, Lee & Messerschmidt]:

input: hpulse is the oversampled channel response, including Tx and Rx bandlimiting

output: "infinite" length equalizer and wmf

1. compute folded signal to noise, FSN

2. G = frequency-domain target

3. WMF, wmf = whitened matched filter

4. BREQ, br = baud-rate equalizer freq-response = G / (FSN + 1/A)

5. pre = conv(br,wmf)




mmse_pre_2

input: hpulse is the oversampled channel response, including Tx and Rx bandlimiting

output: "infinite" length equalizer and wmf

1. wmf = whitened matched filter = hpulse(-t)

2. F = sample conv(hpulse,wmf) & take fft

3. g, G = target

4. C, c = equalizer = [G conj(F)] / [conj(F) + 1/A]



mmse_pre_t2

input: hpulse is the oversampled channel response, including Tx and Rx bandlimiting

output: finite length fractionally-spaced equalizer

1. sample pulse response at its peak

2. N = the equalizer length (must be even)

3. F = channel convolution matrix for fractional tap spacing

4. Eh = target

5. c = equalizer = Eh * F' * inv(F * F' + SNR * eye(N))



mmse_pre_t

input: hpulse is the oversampled channel response, including Tx and Rx bandlimiting

output: finite length baud-rate equalizer (can also be used for jitter-aware equalizer if T/2 spacing is used with Nyquist-2 target)

1. sample pulse response at its peak

2. N = the equalizer length (must be even)

3. F = channel convolution matrix for baud-rate tap spacing

4. Eh = target

5. c = equalizer = Eh * F' * inv(F * F' + SNR * eye(N))

************
Matlab ABCD two-port analysis functions
Version: 2.0
Date: April 1, 2007
Author: Tony Chan Carusone, ECE Dept, U. of Toronto
************

The files in this package are intended for simulating 2-port electrical systems and channel models in Matlab.  They use frequency-dependent transmission (ABCD) matrices to describe the behavior of two-port networks.  Help for each particular function is accessible by typing "help <function-name>" at the matlab command prompt.  Here are some general instructions.

All vectors passed to functions in this package should be row vectors.

A frequency vector, f, should be defined storing the frequency points at which you would like the two-port network defined.  The frequency points need not be linearly spaced (or even in increasing order), although this is useful if you later wish to perform an FFT to get a time-domain step response.

Two-port network descriptions are stored in structures with 4 elements:  x.A, x.B, x.C, x.D

Each element of the structure is a row vector, of length equal to length(f), containing the A, B, C, and D parameters (respectively) at each frequency in f.

The structures may be passed to functions without reference to the underlying ABCD elements.  For example:

%% --------
f = linspace(0,1e9,1000);
x = impedance(10,f);
y = admittance(0.5,f);
z = series(x,y);
%% --------

This structures x & y are two-port networks consisting of simple impedances and admittances.  The structure z is a series combination of x & y.

To get the frequency response of a two-port network, there are 2 options:

1. Use the "terminate" or "src_term" functions, which allow you to terminate the 2-port in a fixed impedances.

2. If the two-port has an open-circuit output, the frequency response for V2/V1 is simply "1./z.A"

Functions are also provided to create 2-port descriptions from (frequency-dependent) transmission line RLGC parameters or from S-parameters.

The file "chmodel.m" provides an example of creating a channel model, in both the frequency-domain and time-domain, using the package.  It starts from a lossy tansmission line description including skin effect and dielectric loss, calculates frequency-dependent RLGC parameters, creates transmission matrices for the transmission line and a simple package model, then combines them and plots the resulting channel response in the time and frequency domains.

LIST OF FILES:
--------------

admittance.m
chmodel.m
freq2impulse.m
impedance.m
rlgc.m
series.m
src_term.m
terminate.m
sparam.m
