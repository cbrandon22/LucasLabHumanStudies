function [] = setFreqYTicks(ax,fBins)
yt = [[1:5:length(fBins)] length(fBins)];
if any(diff(yt)==0)
    yt = [1:5:length(fBins)];
end
set(ax,'xtick',yt,...
'xticklabel',round(fBins(yt)))