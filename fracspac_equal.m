%  fracspac_equal.m               Jack Kurzweil

% Computing Fractionally-Spaced Ts/2 Bandpass Equalizer Performance

%		This script takes as its input the complex Ts/2 (fractionally
%spaced samples of the baseband equivalent of a bandpass channel impulse
%response and computes the optimum transversal equalizer tap 
%coefficients.

%		The user is able to specify the number of equalizer tap coefficients
%and the signal-to-noise ratio.


%		The input to the script is a file such as 

%    fschan = [.05 .10 1.00 .06 -.03; .20 -.48 -.05 .12 .05]'

%where the two columns are the in-phase and quadrature channel samples.

%		The user may specify their own files and would then run the script
%according to the procedure:

%    	>>fschan = [.05 .10 1.00 .06 -.03; .20 -.48 -.05 .12 .05]'
%		>>fracspac_equal		

%    	Accompanying this script is a collection of channel models based 
%on the B3 telephone line model.  These files are generically labeled
%b3chxf.m where x = 1,2,3,4,5 representing five different sampling
%phases and 'f' indicates a fractionally-spaced Ts/2 system.  There is
%also a channel model labeled notch_f.m that is the time domain equivalent
%of a channel having an in-band 20db notch.  These channel models can used with the
%then the following procedure:
%                       >>load b3ch3f.m
%                       >>fschan = b3ch3f
%                       >>fracspac_equal





fschan;

clear Hp Hq H0p H0q H1p H1q H0 H1 Ap Aq C0 C1 E0 E1

x = size( fschan );
g = x(1,1);


%separate into real and imaginary vectors
for k = 1:g
   Hp(1,k) = fschan(k,1);
   Hq(1,k) = fschan(k,2);
end
Hp;
Hq;


%separate into the on-beat 0 and the off-beat 1 vectors
for j = 1:(g/2)
   H0p(1,j) = Hp(1,(2*j - 1) );
   H1p(1,j) = Hp(1,(2*j) );   
   
   H0q(1,j) = Hq(1,(2*j - 1) );
   H1q(1,j) = Hq(1,(2*j) );   
end

%these quantities are used to construct the convolution matrix
Hp = [H0p;H1p];
Hq = [H0q;H1q];

%these quantities are used to compute the on-beat and off-beat outputs
H0 = H0p + i*H0q;
H1 = H1p + i*H1q;

disp('The channel length is')
g
disp('The equalizer length should be an even number')
%m = input('equalizer length m = ')   
m = 14;
Ap = Hp;
Aq = Hq;
for j = 1:((m/2)-1) 
    Ap = [Ap zeros(2*j,1); zeros(2,j) Hp];
    Aq = [Aq zeros(2*j,1); zeros(2,j) Hq]; 
end
Ap;
Aq;


A = Ap+i*Aq;
Auto = A * A';

B = eye(size(Auto));
%SNR = input('SNR in db = ')
SNR = 30
sigsq = exp(-(SNR/10.0)*log(10));

B = sigsq*B;

Auto = Auto + B;


Eh = zeros(1, (m+g-2)/2 );
x=length(Eh);
Eh(1,ceil(x/2)) = 1;
Eh;
B = Eh*A';

Cconj = B/Auto;

MSE = 1 - (Cconj * B'); 

%separate into the on-beat 0 and the off-beat 1 vectors
for j = 1:(m/2)
   C0(1,j) = Cconj(1,(2*j - 1) );
   C1(1,j) = Cconj(1,(2*j) );   
end
C0;
C1;
%pause

E0 = conv(H0,C0);
E1 = conv(H1,C1);
error = E0 + E1 - Eh;
error = error';


disp('The mean square ISI is')
ISI = norm(error)^2 
disp('The ISI plus noise is')
MSE

A1 = (norm(C0)^2 + norm(C1)^2);
A2 = (norm(H0)^2 ); 
A3 = (norm(E0)^2 + norm(E1)^2); 
disp('The noise amplification is:')
namp = (MSE-ISI)/sigsq

