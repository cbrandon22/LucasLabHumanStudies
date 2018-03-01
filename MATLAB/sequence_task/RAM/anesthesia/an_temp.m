%temp
for i = 3
    hold all
    lfp = an_getlfp_ms_wrapper(i,soundEvents,1500,-500,0);
    plotSEM(-499:1000,nanmean(lfp,1),SEM(lfp))
    xlabel('Time from Sound (ms)')
    ylabel ('voltage (mv)')
    set(gca,'ylim',[-100 100])
    plot([0 0],get(gca,'ylim'))
    pause
    clf
    
end