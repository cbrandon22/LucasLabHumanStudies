statDir = '/Volumes/Lucas_ECoG_Data/scratch/PHASE/stats/2sec';
cd([statDir '/broadTheta250_broadTheta300']);
load('HUP091_clus1_broadTheta250_moveBoot.mat');
load('/Volumes/Lucas_ECoG_Data/scratch/PHASE/k2/broadTheta250_broadTheta300/HUP091_clusPLV_PW_clus1.mat')
offset=-2;
dat = mean(PW_ThetaMove,1);
dat = dat(1:end-512);
startInd = length(dat)-length(shuffles)+1;
dat = dat(startInd:end);

for i=1:10
    shuffDat = dat(shuffles(i,:));
    startInd = length(dat)-length(shuffles)+1;
    t=(1:length(shuffDat))/512;
    t=t+offset;
    hold on
    title(['HUP091 broadTheta shuffled PLV' num2str(i)])
    plot(t,shuffDat,'g')
    plot(t,dat,'r')
    xlim([-2 6.25])
    ylim([0 1])
    plot([0 0], [-7 7], 'b-', 'linewidth', 2)
    plot([-3 7], [0 0], 'k--', 'linewidth', 2)
    plot([0 5], [-1 -1], 'k-', 'linewidth' ,2)
    xlabel('Time (s)')
    ylabel('Shuffled PLV')
    keyboard;
    close;
end