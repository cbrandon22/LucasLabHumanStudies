%temp
for i = 1:128
    hold all
    lfp = an_getlfp_ms_wrapper(i,correct,1500,-500,0);
    plotSEM(-499:1000,nanmean(lfp,1),SEM(lfp))
    xlabel('Time from Sound (ms)')
    ylabel ('voltage (mv)')
    title(num2str(i));
    set(gca,'ylim',[-100 100])
    plot([0 0],get(gca,'ylim'))
    pause
    clf
    
end