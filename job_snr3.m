IEEE802p3aqChModel

% wmf = hpulse(end:-1:1)./10^-(SNR/10); fpulse = conv(hpulse,wmf); % include matched filter  (assuming white noise)
fpulse = hpulse; % skip matched filter

% 802.3aq spec calls for BER < 1e-12 with input SNRs as follows:
% short channel: SNR = 20*log10(26.3) = 28.4 dB
% pre-cursor:    SNR = 20*log10(45.6) = 33.2 dB
% symmetrical:   SNR = 20*log10(31.4) = 31.4 dB
% post-cursor:   SNR = 20*log10(33.4) = 33.4 dB
%
% The signal to noise ratio measurement in the spec is the ratio of the
% maximum received amplitude squared (i.e. for long string of 1's or 0's)
% to the received noise power

% tight lower bound on the BER of an ideal MLSE receiver
% BERMLSE = Qfun(10^(SNR/20));

%txjitter = 0.033*Tbaud; % max permitted rms transmitter jitter from spec
%rxjitter = 5e-12; % not specified
txjitter = 0;
rxjitter = 0;
Nbins = 100;

Nlist = 1:4;
SNRspec = 25;
clear BER_dfe0 BER_dfe1 BER_dfe2 BER_dfe3

% ------
% baud-spaced equalizer with no FBE
% ------

for k = 1:length(Nlist)
    N = Nlist(k)
    
    % tap spacing
    T = OSR;

    SNRbest = -Inf;
    
    for eqdelay = 0:N % approx delay through equalizer in # of taps
    for phase = 1:OSR
        
    % sample the channel response
    fpulsesampled = fpulse(phase:T:end); % sample the pulse response

    % signal to noise ratio at the receiver input
    % is penalized by the ISI
    SNR = SNRspec + 10*log10(sum(fpulsesampled.^2)/2); % dB
    npwr = 10^(-SNRspec/10);

    [c,cheq] = lineq(N,fpulsesampled,SNR,1,eqdelay);
    
    % find resulting ISI penalty and noise amplification at 
    % the equalizer output
    sig = max(cheq);
    desired = zeros(1,length(cheq));
    desired(find(cheq==max(cheq),1)) = sig;
    ISI = cheq-desired;
    mse = sum(ISI.^2) + npwr*sum(c.^2);
    if 10*log10(sig^2/mse) > SNRbest
        SNRbest = 10*log10(sig^2/mse);
        c_br = c;
        cheq_br = cheq;
    end
    end
    end
    
    [A,B] = eye_stats(himp(380:520),c_br,0,2,txjitter,rxjitter,Ts,OSR,T,Nbins,sqrt(npwr));
    BER_dfe0(k) = min(min(B));
    SNRbest0(k) = SNRbest;
    B0 = B;
end
    
% ------
% baud-spaced equalizer with 1-tap DFE
% ------

Ndfe = 1
    
for k = 1:length(Nlist)
    N = Nlist(k)
    
    % tap spacing
    T = OSR;
    
    SNRbest = -Inf;

    for eqdelay = 0:N % approx delay through equalizer in # of taps
    for phase = 1:OSR
        
    % sample the channel response
    fpulsesampled = fpulse(phase:T:end); % sample the pulse response

    % signal to noise ratio at the receiver input
    % is penalized by the ISI
    SNR = SNRspec + 10*log10(sum(fpulsesampled.^2)); % dB
    npwr = 10^(-SNRspec/10);

    [c,dfe,cheq] = dfe_lineq(N,Ndfe,fpulsesampled,SNR,eqdelay);
    
    % find resulting SNR at equalizer output (ratio of signal to noise &
    % ISI power)
    delay = find(cheq == max(cheq));
    eyeamp = sum(cheq);
    sig = cheq(delay);
    desired = zeros(1,length(cheq));
    desired(delay) = sig;
    ISI = cheq-desired;
    mse = sum(ISI.^2) + npwr*sum(c.^2);
    if 10*log10(sig^2/mse) > SNRbest
        SNRbest = 10*log10(sig^2/mse);
        c_dfe = c;
        dfe_taps = dfe';
        eqdelay_dfe1(k) = eqdelay;
        phase_dfe1(k) = phase;
    end
    end
    end

    [A,B] = eye_stats(himp(380:520),c_dfe,dfe_taps,2,txjitter,rxjitter,Ts,OSR,T,Nbins,sqrt(npwr));
    BER_dfe1(k) = min(min(B));
    SNRbest1(k) = SNRbest;
end

    
% ------
% baud-spaced equalizer with 2-tap DFE
% ------

Ndfe = 2
    
for k = 1:length(Nlist)
    N = Nlist(k)
    
    % tap spacing
    T = OSR;
    
    SNRbest = -Inf;

    for eqdelay = 0:N % approx delay through equalizer in # of taps
    for phase = 1:OSR
        
    % sample the channel response
    fpulsesampled = fpulse(phase:T:end); % sample the pulse response

    % signal to noise ratio at the receiver input
    % is penalized by the ISI
    SNR = SNRspec + 10*log10(sum(fpulsesampled.^2)); % dB
    npwr = 10^(-SNRspec/10);

    [c,dfe,cheq] = dfe_lineq(N,Ndfe,fpulsesampled,SNR,eqdelay);
    
    % find resulting SNR at equalizer output (ratio of signal to noise &
    % ISI power)
    delay = find(cheq == max(cheq));
    eyeamp = sum(cheq);
    sig = cheq(delay);
    desired = zeros(1,length(cheq));
    desired(delay) = sig;
    ISI = cheq-desired;
    mse = sum(ISI.^2) + npwr*sum(c.^2);
    if 10*log10(sig^2/mse) > SNRbest
        SNRbest = 10*log10(sig^2/mse);
        c_dfe = c;
        dfe_taps = dfe';
        eqdelay_dfe2(k) = eqdelay;
        phase_dfe2(k) = phase;
    end
    end
    end

    [A,B] = eye_stats(himp(380:520),c_dfe,dfe_taps,2,txjitter,rxjitter,Ts,OSR,T,Nbins,sqrt(npwr));
    BER_dfe2(k) = min(min(B));
    SNRbest2(k) = SNRbest;
end

% ------
% baud-spaced equalizer with 3-tap DFE
% ------

Ndfe = 3
    
for k = 1:length(Nlist)
    N = Nlist(k)
    
    % tap spacing
    T = OSR;
    
    SNRbest = -Inf;

    for eqdelay = 0:N % approx delay through equalizer in # of taps
    for phase = 1:OSR
        
    % sample the channel response
    fpulsesampled = fpulse(phase:T:end); % sample the pulse response

    % signal to noise ratio at the receiver input
    % is penalized by the ISI
    SNR = SNRspec + 10*log10(sum(fpulsesampled.^2)); % dB
    npwr = 10^(-SNRspec/10);

    [c,dfe,cheq] = dfe_lineq(N,Ndfe,fpulsesampled,SNR,eqdelay);
    
    % find resulting SNR at equalizer output (ratio of signal to noise &
    % ISI power)
    delay = find(cheq == max(cheq));
    eyeamp = sum(cheq);
    sig = cheq(delay);
    desired = zeros(1,length(cheq));
    desired(delay) = sig;
    ISI = cheq-desired;
    mse = sum(ISI.^2) + npwr*sum(c.^2);
    if 10*log10(sig^2/mse) > SNRbest
        SNRbest = 10*log10(sig^2/mse);
        c_dfe = c;
        dfe_taps = dfe';
        eqdelay_dfe3(k) = eqdelay;
        phase_dfe3(k) = phase;
        fpulsesampled_dfe3 = fpulsesampled;
    end
    end
    end

    [A,B] = eye_stats(himp(380:520),c_dfe,dfe_taps,2,txjitter,rxjitter,Ts,OSR,T,Nbins,sqrt(npwr));
    BER_dfe3(k) = min(min(B));
    SNRbest3(k) = SNRbest;
end

semilogy(Nlist,BER_dfe0,'k-o',Nlist,BER_dfe1,'g-<',Nlist,BER_dfe2,'b-d',Nlist,BER_dfe3,'r-s','linewidth',2); set(gca,'linewidth',3,'fontsize',18,'fontname','times'); grid on
