function [moveWaitStat,moveMoveStat,waitWaitStat,winStats] = PLV_mm_mesPlot(moveDat,waitDat,t,sampleWin,sampleStep,cueInd,tit,iterations)
%%% Calculates hedges' g and plots with raw plv

%Create time bins
sampleBins = sampleWin:sampleStep:size(moveDat);
timeBins = ((sampleBins-0.5*sampleWin)/512)-2;

for tt = 1:length(sampleBins)
    %identify time windows
    sInd = (sampleBins(tt)-sampleWin+1):sampleBins(tt);
    binMoveMat(:,tt) = moveDat(sInd);
    binWaitMat(:,tt) = waitDat(sInd);
end

stats = mes(binMoveMat,binWaitMat,'hedgesg','nBoot',iterations);
preCueStats = mes(moveDat(1:cueInd),waitDat(1:cueInd),'hedgesg','nBoot',iterations);
postCueStats = mes(moveDat(cueInd+1:end),waitDat(cueInd+1:end),'hedgesg','nBoot',iterations);
moveStats = mes(moveDat(cueInd+1:end),moveDat(1:cueInd),'hedgesg','nBoot',iterations);
waitStats = mes(waitDat(cueInd+1:end),waitDat(1:cueInd),'hedgesg','nBoot',iterations);
moveWaitStat = [postCueStats.hedgesg,postCueStats.hedgesgCi(1),postCueStats.hedgesgCi(2)];
moveMoveStat = [moveStats.hedgesg,moveStats.hedgesgCi(2),moveStats.hedgesgCi(2)];
waitWaitStat = [waitStats.hedgesg,waitStats.hedgesgCi(2),waitStats.hedgesgCi(2)];
winStats = {stats.hedgesg,stats.hedgesgCi(2),stats.hedgesgCi(2),timeBins};

figure
hold on;
subplot(2,1,1)
hold on;
title(tit)
xlabel('Time (s)')
ylabel('hedges'' g')
plot([0 0], [-1.5 1.5],'b')
plot([-0.75 0.75],[moveStats.hedgesg moveStats.hedgesg],'g')
plot([-0.75 0.75],[waitStats.hedgesg waitStats.hedgesg],'r')
plot([-0.75 0.75],[moveStats.hedgesgCi(1) moveStats.hedgesgCi(1)],'g--')
plot([-0.75 0.75],[moveStats.hedgesgCi(2) moveStats.hedgesgCi(2)],'g--')
plot([-0.75 0.75],[waitStats.hedgesgCi(1) waitStats.hedgesgCi(1)],'r--')
plot([-0.75 0.75],[waitStats.hedgesgCi(2) waitStats.hedgesgCi(2)],'r--')
plot([-2 0],[preCueStats.hedgesg preCueStats.hedgesg],'b')
plot([0 5.25],[postCueStats.hedgesg postCueStats.hedgesg],'b')
plot([-2 0],[preCueStats.hedgesgCi(1) preCueStats.hedgesgCi(1)],'b--')
plot([-2 0],[preCueStats.hedgesgCi(2) preCueStats.hedgesgCi(2)],'b--')
plot([0 5.25],[postCueStats.hedgesgCi(1) postCueStats.hedgesgCi(1)],'b--')
plot([0 5.25],[postCueStats.hedgesgCi(2) postCueStats.hedgesgCi(2)],'b--')
scatter(timeBins,stats.hedgesg,'k','filled')
for hgind = 1:length(stats.hedgesgCi)
    plot([timeBins(hgind) timeBins(hgind)],stats.hedgesgCi(:,hgind),'k')
end
xlim([-2 5.25])
plot([-2 length(stats.hedgesg)], [0 0],'k')
subplot(2,1,2)
hold on;
title('normalized PLV')
plot(t,fgsmooth(moveDat,10),'g')
plot(t,fgsmooth(waitDat,10),'r')
xlim([-2 5.25])
plot([0 0], [-1 1],'b')

end