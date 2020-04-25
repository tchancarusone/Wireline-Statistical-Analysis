function [A,B] = eye_stats(himp,ffe,dfe,Nlevels,txjitter,rxjitter,Ts,OSR,TapSpacing,N,nrms)
% function [A,B] = eye_stats(himp,ffe,dfe,Nlevels,txjitter,rxjitter,Ts,OSR,TapSpacing,Nbins,ndist)
%
% himp = channel impulse response not including equalizer
% ffe = linear equalizer coefficients
% dfe = feedback equalizer coefficients (0 for no DFE)
% Nlevels = number of modulations levels (2=binary, 4=4-PAM, etc.)
% txjitter = transmitter rms jitter (in sec.)
% rxjitter = receiver rms jitter (in sec.)
% Ts = sampling time-step of the step response
% OSR = Baud Time / Ts (oversampling ratio)
% TapSpacing = Delay per tap / Ts
% Nbins = number of vertical bins in eye diagrams
% nrms = rms of AWGN at equalizer input
%
% A = PDF eye of equalizer output
% B = CDF eye of equalizer output

% bounded noise assumption on transmitter and receiver jitter and noise
numstdevs = 5000;

% equalized pulse response
hpulse = conv(himp,ones(1,OSR));
pr = conv(hpulse,upfirdn(ffe,1,TapSpacing,1));
himp = [himp' zeros(1,length(pr)-length(himp))]';

% define vertical bins
minmax = 1.5*sum(abs(pr))/OSR;
bins = linspace(-minmax,minmax,N);

% length of pulse response
M = floor(length(pr)/OSR);
delay = ceil(find(abs(pr-1) == min(abs(pr-1)),1)/OSR);

A = zeros(N,OSR);
for k = 1:OSR
    k;
    % compute sampled pulse response
    p = pr(k:OSR:M*OSR);
    
    % include effect of the DFE assuming no error propagation
    %delay = find(abs(p-1) == min(abs(p-1)),1);
    p(delay+1:delay+length(dfe)) = p(delay+1:delay+length(dfe)) - dfe;
    
    % compute sampled impulse response
    h = himp(k+OSR:OSR:M*OSR);
    
    for l = 1:Nlevels^M
        u = ((dec2bin(l-1,M) == '1') - .5)*2;
        du = u(2:end) - u(1:end-1);
        
        % compute output of channel with ISI only (no jitter or noise)
        x = u*p(end:-1:1);
        b = min(find(abs(bins-x) == min(abs(bins-x)),1));
        
        % compute the voltage noise due to transmitter jitter
        if txjitter ~= 0
            txjitter_noise_pwr = 0;
            for m = 1:length(du)
                txjitter_noise_pwr = txjitter_noise_pwr + ((txjitter/Ts)*abs((du(m)*h(end+1-m))))^2;
            end
            if (txjitter_noise_pwr ~= 0)
                txjitter_noise_dist = exp(-(bins.^2./txjitter_noise_pwr)./2);
                txjitter_noise_dist(find(abs(bins)>numstdevs*sqrt(txjitter_noise_pwr))) = 0; % bound tx jitter
                if (sum(txjitter_noise_dist) ~= 0)
                    txjitter_noise_dist = txjitter_noise_dist./sum(txjitter_noise_dist);
                else
                    txjitter_noise_dist = zeros(1,N); txjitter_noise_dist(ceil(N/2)) = 1;
                end
            else
                txjitter_noise_dist = zeros(1,N); txjitter_noise_dist(ceil(N/2)) = 1;
            end
        else
            txjitter_noise_dist = zeros(1,N); txjitter_noise_dist(ceil(N/2)) = 1;
        end
        
        if b > ceil(N/2)
            A(:,k) = A(:,k) + [zeros(b-ceil(N/2),1); txjitter_noise_dist(1:end+ceil(N/2)-b)']./(Nlevels^M);
        else
            A(:,k) = A(:,k) + [txjitter_noise_dist(ceil(N/2)+1-b:end)'; zeros(ceil(N/2)-b,1)]./(Nlevels^M);
        end
        
    end
    
end

% include noise
if nrms == 0
    ndist = zeros(1,N);
    ndist(ceil(N/2)) = 1;
else
    nrmseq = nrms*sqrt(sum(ffe.^2));
    ndist = exp(-(bins.^2)./(nrmseq^2));
    ndist(find(abs(bins)>numstdevs*nrmseq)) = 0;
    ndist = ndist./sum(ndist);
end

for k = 1:OSR
    temp = conv(A(:,k),ndist);
    A(:,k) = temp(ceil(N/2):end-floor(N/2));
end   

% include receiver jitter
if (rxjitter == 0)
    rxjitter_dist = zeros(1,OSR); rxjitter_dist(round(OSR/2)) = 1;
else
    tnorm = (1-floor(OSR/2):ceil(OSR/2))./(rxjitter/Ts);
    rxjitter_dist = exp(-(tnorm.^2)./2);
    rxjitter_dist(find(abs(tnorm)>numstdevs)) = 0;
    rxjitter_dist = rxjitter_dist./sum(rxjitter_dist);
end
for k = 1:N
    temp = conv([A(k,:) A(k,:)],rxjitter_dist);
    A(k,:) = temp(OSR+1:2*OSR);
end

% create CDF from PDF
B = zeros(floor(N/2),OSR);
for k = 1:OSR
    for l = 1:floor(N/2)
        B(l,k) = sum(A(l:floor(N/2),k));
    end
end
B = [B; B(end:-1:1,:)];
