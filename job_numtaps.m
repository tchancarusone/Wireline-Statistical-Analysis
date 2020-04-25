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
SNRspec = 28.4; % dB
npwr = 10^(-SNRspec/10);

% tight lower bound on the BER of an ideal MLSE receiver
% BERMLSE = Qfun(10^(SNR/20));

Ntaps = 2:2:14;
ISIpen_br = Inf*ones(1,length(Ntaps));
ISIpen_fs = Inf*ones(1,length(Ntaps));
ISIpen_fsdfe1 = Inf*ones(1,length(Ntaps));
ISIpen_fsdfe2 = Inf*ones(1,length(Ntaps));
clear Namp_br Namp_fs Namp_fsdfe1 Namp_fsdfe2
clear phase_br phase_fs phase_fsdfe1 phase_fsdfe2
clear eqdelay_br eqdelay_fs eqdelay_fsdfe1 eqdelay_fsdfe2

% ------
% T/2-spaced equalizer
% ------
    
for k = 1:length(Ntaps)
    N = Ntaps(k);

    % tap spacing
    T = OSR/2;

    for eqdelay = 0:2*N % approx delay through equalizer in # of taps
    for phase = 1:OSR
        
    % sample the channel response
    fpulsesampled = fpulse(phase:T:end); % sample the pulse response

    % signal to noise ratio at the receiver input
    % is penalized by the ISI
    SNR = SNRspec + 10*log10(sum(fpulsesampled.^2)); % dB

    [c,cheq] = fseq_tover2(N,fpulsesampled,SNR,1,eqdelay);
    % [c,cheq] = lineq(N,fpulsesampled,SNR,[.5 1 .5],eqdelay); % jitter-aware
    
    % find resulting ISI penalty and noise amplification at 
    % the equalizer output
    eyeamp = sum(cheq);
    sig = max(cheq);
    desired = zeros(1,length(cheq));
    desired(find(cheq==max(cheq),1)) = sig;
    ISI = sum(abs(cheq-desired));
    if (sig>ISI) & -10*log10((sig-ISI)/eyeamp) < ISIpen_fs(k)
        ISIpen_fs(k) = -10*log10((sig-ISI)/eyeamp); 
        Namp_fs(k) = sqrt(sum(c.^2));
        eqdelay_fs(k) = eqdelay;
        phase_fs(k) = phase;
    end
    end
    end
end
    
% ------
% T/2-spaced equalizer with 1-tap DFE
% ------
    
Ndfe = 1;
for k = 1:length(Ntaps)
    N = Ntaps(k);

    % tap spacing
    T = OSR/2;

    for eqdelay = 0:2*N % approx delay through equalizer in # of taps
    for phase = 1:OSR
        
    % sample the channel response
    fpulsesampled = fpulse(phase:T:end); % sample the pulse response

    % signal to noise ratio at the receiver input
    % is penalized by the ISI
    SNR = SNRspec + 10*log10(sum(fpulsesampled.^2)); % dB

    [c,dfe,cheq] = dfe_tover2ffe(N,Ndfe,fpulsesampled,SNR,eqdelay);
    
    % find resulting SNR at equalizer output (ratio of signal to noise &
    % ISI power)
    delay = round(find(fpulsesampled==max(fpulsesampled),1)/2 + eqdelay/2);
    cheq(delay+1:delay+Ndfe) = cheq(delay+1:delay+Ndfe) - dfe;
    eyeamp = sum(cheq);
    sig = cheq(delay);
    desired = zeros(1,length(cheq));
    desired(delay) = sig;
    ISI = sum(abs(cheq-desired));
    if (sig>ISI) & -10*log10((sig-ISI)/eyeamp) < ISIpen_fsdfe1(k)
        ISIpen_fsdfe1(k) = -10*log10((sig-ISI)/eyeamp); 
        Namp_fsdfe1(k) = sqrt(sum(c.^2));
        eqdelay_fsdfe1(k) = eqdelay;
        phase_fsdfe1(k) = phase;
    end
    end
    end
end
    
% ------
% T/2-spaced equalizer with 2-tap DFE
% ------
    
Ndfe = 2;
for k = 1:length(Ntaps)
    N = Ntaps(k);

    % tap spacing
    T = OSR/2;

    for eqdelay = 0:2*N % approx delay through equalizer in # of taps
    for phase = 1:OSR
        
    % sample the channel response
    fpulsesampled = fpulse(phase:T:end); % sample the pulse response

    % signal to noise ratio at the receiver input
    % is penalized by the ISI
    SNR = SNRspec + 10*log10(sum(fpulsesampled.^2)); % dB

    [c,dfe,cheq] = dfe_tover2ffe(N,Ndfe,fpulsesampled,SNR,eqdelay);
    
    % find resulting SNR at equalizer output (ratio of signal to noise &
    % ISI power)
    delay = round(find(fpulsesampled==max(fpulsesampled),1)/2 + eqdelay/2);
    cheq(delay+1:delay+Ndfe) = cheq(delay+1:delay+Ndfe) - dfe;
    eyeamp = sum(cheq);
    sig = cheq(delay);
    desired = zeros(1,length(cheq));
    desired(delay) = sig;
    ISI = sum(abs(cheq-desired));
    if (sig>ISI) & -10*log10((sig-ISI)/eyeamp) < ISIpen_fsdfe2(k)
        ISIpen_fsdfe2(k) = -10*log10((sig-ISI)/eyeamp); 
        Namp_fsdfe2(k) = sqrt(sum(c.^2));
        eqdelay_fsdfe2(k) = eqdelay;
        phase_fsdfe2(k) = phase;
    end
    end
    end
end

% ------
% baud-rate equalizer
% ------
    
for k = 1:length(Ntaps)
    N = Ntaps(k)/2;

    % tap spacing
    T = OSR;

    % baud-rate equalizer
    for eqdelay = 0:2*N % approx delay through equalizer in # of taps
    for phase = 1:1:OSR

    % sample the channel response
    fpulsesampled = fpulse(phase:T:end); % sample the pulse response

    [c,cheq] = lineq(N,fpulsesampled,SNR,1,eqdelay);

    % find resulting ISI penalty and noise amplification at 
    % the equalizer output
    eyeamp = sum(cheq);
    sig = max(cheq);
    desired = zeros(1,length(cheq));
    desired(find(cheq==max(cheq),1)) = sig;
    ISI = sum(abs(cheq-desired));
    if (sig>ISI) & -10*log10((sig-ISI)/eyeamp) < ISIpen_br(k)
        ISIpen_br(k) = -10*log10((sig-ISI)/eyeamp); 
        Namp_br(k) = sqrt(sum(c.^2));
        eqdelay_br(k) = eqdelay;
        phase_br(k) = phase;
    end
    end
    end
end

subplot(2,1,1)
plot(Ntaps./2,ISIpen_br,'bd-',Ntaps./2,ISIpen_fs,'ro-',Ntaps./2,ISIpen_fsdfe1,'mx-',Ntaps./2,ISIpen_fsdfe2,'gs-','linewidth',2,'markersize',10);
set(gca,'fontsize',14,'fontname','times','linewidth',2)
xlabel('FFE Span [UI]'); ylabel('ISI Penalty [dB]'); grid on
legend('Baud-Rate FFE, no DFE','T_b/2-spaced FFE, no DFE','T_b/2-spaced FFE, 1-tap DFE',''T_b/2-spaced FFE, 2-tap DFE')
subplot(2,1,2)
plot(Ntaps./2,20*log10(Namp_br),'b',Ntaps./2,20*log10(Namp_fs),'r',Ntaps./2,20*log10(Namp_fsdfe1),'m',Ntaps./2,20*log10(Namp_fsdfe2),'g','linewidth',2,'markersize',10);
