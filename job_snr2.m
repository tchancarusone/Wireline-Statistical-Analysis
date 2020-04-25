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
Nbins = 30;

N = 6;
Ndfe = 2;
SNRspec = 20:1:35;

clear BER_fs BER_dfefs BER_dfe BER_mlse BER_br

% ------
% T/2-spaced equalizer
% ------
    
for k = 1:length(SNRspec)

    % tap spacing
    T = OSR/2;

    SNRbest = -Inf;
    
    for eqdelay = 0:2*N % approx delay through equalizer in # of taps
    for phase = 1:OSR
        
    % sample the channel response
    fpulsesampled = fpulse(phase:T:end); % sample the pulse response

    % signal to noise ratio at the receiver input
    % is penalized by the ISI
    SNR = SNRspec(k) + 10*log10(sum(fpulsesampled.^2)/2); % dB
    npwr = 10^(-SNRspec(k)/10);

    [c,cheq] = fseq_tover2(N,fpulsesampled,SNR,1,eqdelay);
    
    % find resulting ISI penalty and noise amplification at 
    % the equalizer output
    delay = find(cheq == max(cheq));
    sig = cheq(delay);
    desired = zeros(1,length(cheq));
    desired(delay) = sig;
    ISI = cheq-desired;
    mse = sum(ISI.^2) + npwr*sum(c.^2);
    if 10*log10(sig^2/mse) > SNRbest
        SNRbest = 10*log10(sig^2/mse);
        c_fs = c;
    end
    end
    end
    
    [A,B] = eye_stats(himp(380:520),c_fs,0,2,txjitter,rxjitter,Ts,OSR,T,Nbins,sqrt(npwr));
    BER_fs(k) = min(min(B));
    if BER_fs(k) < 1e-20
        break
    end
end

% ------
% T/2-spaced equalizer with DFE
% ------
    
for k = 1:length(SNRspec)

    % tap spacing
    T = OSR/2;
    
    SNRbest = -Inf;

    for eqdelay = 0:2*N % approx delay through equalizer in # of taps
    for phase = 1:OSR
        
    % sample the channel response
    fpulsesampled = fpulse(phase:T:end); % sample the pulse response

    % signal to noise ratio at the receiver input
    % is penalized by the ISI
    SNR = SNRspec(k) + 10*log10(sum(fpulsesampled.^2)/2); % dB
    npwr = 10^(-SNRspec(k)/10);

    [c,dfe,cheq] = dfe_tover2ffe(N,Ndfe,fpulsesampled,SNR,eqdelay);
    
    % find resulting SNR at equalizer output (ratio of signal to noise &
    % ISI power)
    delay = find(cheq == max(cheq));
    sig = cheq(delay);
    desired = zeros(1,length(cheq));
    desired(delay) = sig;
    ISI = cheq-desired;
    mse = sum(ISI.^2) + npwr*sum(c.^2);
    if 10*log10(sig^2/mse) > SNRbest
        SNRbest = 10*log10(sig^2/mse);
        c_dfefs = c;
        dfefs_taps = dfe';
    end
    end
    end

    [A,B] = eye_stats(himp(380:520),c_dfefs,dfefs_taps,2,txjitter,rxjitter,Ts,OSR,T,Nbins,sqrt(npwr));
    BER_dfefs(k) = min(min(B));
    if BER_dfefs(k) < 1e-20
        break
    end
    
end

% ------
% baud-spaced equalizer
% ------

for k = 1:length(SNRspec)

    % tap spacing
    T = OSR;

    SNRbest = -Inf;
    
    for eqdelay = 0:N % approx delay through equalizer in # of taps
    for phase = 1:OSR
        
    % sample the channel response
    fpulsesampled = fpulse(phase:T:end); % sample the pulse response

    % signal to noise ratio at the receiver input
    % is penalized by the ISI
    SNR = SNRspec(k) + 10*log10(sum(fpulsesampled.^2)/2); % dB
    npwr = 10^(-SNRspec(k)/10);

    [c,cheq] = lineq(round(N/2),fpulsesampled,SNR,1,eqdelay);
    
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
    end
    end
    end
    
    [A,B] = eye_stats(himp(380:520),c_br,0,2,txjitter,rxjitter,Ts,OSR,T,Nbins,sqrt(npwr));
    BER_br(k) = min(min(B));
    if BER_br(k) < 1e-20
        break
    end
end
    
% ------
% baud-spaced equalizer with DFE
% ------
    
for k = 1:length(SNRspec)

    % tap spacing
    T = OSR;
    
    SNRbest = -Inf;

    for eqdelay = 0:N % approx delay through equalizer in # of taps
    for phase = 1:OSR
        
    % sample the channel response
    fpulsesampled = fpulse(phase:T:end); % sample the pulse response

    % signal to noise ratio at the receiver input
    % is penalized by the ISI
    SNR = SNRspec(k) + 10*log10(sum(fpulsesampled.^2)); % dB
    npwr = 10^(-SNRspec(k)/10);

    [c,dfe,cheq] = dfe_lineq(round(N/2),Ndfe,fpulsesampled,SNR,eqdelay);
    
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
    end
    end
    end

    [A,B] = eye_stats(himp(380:520),c_dfe,dfe_taps,2,txjitter,rxjitter,Ts,OSR,T,Nbins,sqrt(npwr));
    BER_dfe(k) = min(min(B));
    if BER_dfe(k) < 1e-20
        break
    end
    
end

semilogy(SNRspec(1:length(BER_br)),BER_br,'k-<',SNRspec(1:length(BER_fs)),BER_fs,'b-d',SNRspec(1:length(BER_dfe)),BER_dfe,'g-o',SNRspec(1:length(BER_dfefs)),BER_dfefs,'r-s','linewidth',2); set(gca,'linewidth',3,'fontsize',18,'fontname','times'); ylim([1e-20 1]); grid on
