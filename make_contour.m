bins = -2:.02:2;
ndist = gausswin(length(bins),max(bins)/(nrms*sqrt(sum(ic.^2))));
[A,B] = eye_stats(hstep2,ic,2,jittrms,1.48e-12,Ts,OSR,T,bins,ndist);
%contour(1/OSR:1/OSR:2,mean([bins(1:end-1);bins(2:end)]),[B B],[1e-15 1e-15],'k');