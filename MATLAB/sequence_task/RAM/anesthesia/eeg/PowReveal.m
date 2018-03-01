% function pow = PowReveal(lfp)

for n = 90:97 %size(lfpCLDH1,1)
    hold all
%     fprintf('%02d%%',fix(n/size(lfp,1)*100));
    periodogram(lfp(n,:));
    pause
    clf
    
end