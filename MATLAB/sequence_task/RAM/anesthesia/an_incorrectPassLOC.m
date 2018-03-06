%temp
for i = 1:128
    hold all
    lfp = an_getlfp_ms_wrapper(i,incorrect,1500,-500,0);
    plotSEM(-499:1000,nanmean(lfp,1),SEM(lfp))
    xlabel('Time from Sound (ms)')
    ylabel ('voltage (mv)')
    title(i);
    set(gca,'ylim',[-100 100])
    plot([0 0],get(gca,'ylim'))
    pause
    clf
    
end