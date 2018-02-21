function CCDT_behavioral2
% Quick look at behavioral data for CCDT
subj = 'HUP143';
sessNum = []; % Pick specific session(s). Default combines all sessions
allRTs = CCDTcatBehavioral(subj,sessNum);
saveFigs = 1;
dirs = le_dirs('CCDT');
saveFigsDir = fullfile(dirs.scratch,'figs',subj,'behavioral');
overallTitle = strcat(subj,' CCDT Behavioral Summary');

%% Organize RTs into subgroups
%if strcmp(subj,'HUP133')
%    filtRTs = allRTs((allRTs(:,1)+allRTs(:,3)>300 & allRTs(:,1)~=1),:); % Filter out very early presses & non-response
%else
filtRTs = allRTs((allRTs(:,1)+allRTs(:,3)>300 & allRTs(:,1)~=1000),:);
%end
posRTs = filtRTs((filtRTs(:,1)>0),:); % Reaction times after CC
guessRTs = filtRTs((filtRTs(:,1)<150),:); % Early/subphysiological presses
fastRTs = posRTs((posRTs(:,1)<250),:); % Fast RTs
slowRTs = posRTs((posRTs(:,1)>750),:); % Slow RTs

%% Scatter plot
rts = filtRTs(:,1); % select subgroup to plot
tit = strcat(subj,' Filtered Reaction Times');

figure('Name',overallTitle,'NumberTitle','off','position',[800 50 800 800]);
subplot(3,2,[1 2]);
scatter([1:length(rts)],rts);
%text(length(rts)/3,max(rts)+500,overallTitle,'FontSize',18)
hold all;
pcoef = polyfit([1:length(rts)]',rts,1);
plot(polyval(pcoef,[1:length(rts)]'));
plot([0 length(rts)],[250 250],'g');
plot([0 length(rts)],[150 150],'g--');
plot([0 length(rts)],[0 0],'k');
title(tit);xlabel('trial number');ylabel('reaction time (ms)');
legend('RT','best-fit','Fast','SubPhys','Location','eastoutside')

%% Bar plot of all RTs
rts = posRTs(:,1);
tit = strcat(subj,' Reaction Time Distribution');
subplot(3,2,3);
[N,edges] = histcounts(rts);
x = edges(1:end-1) + .5*(unique(diff(edges)));
H = bar(x,N,'hist');
set(H,'facecolor',[.3 .3 .3]);
title(tit);xlabel('reaction time (ms)');ylabel('count');

%% Average RT per 9-target loop
rts = filtRTs(:,1); % select subgroup to plot
tit = strcat(subj,' Average RT per Loop');

goodTrial = 0;
loopRTs=[];
loopAvg=[];
for i=1:length(rts)
    if rts(i)>0 && rts(i)<1000
        loopRTs = [loopRTs;rts(i)];
        goodTrial = goodTrial+1;
    end
    if goodTrial == 9;
        loopAvg = [loopAvg; nanmean(loopRTs)];
        goodTrial = 0;
        loopRTs = [];
    end
end
subplot(3,2,4);
scatter(1:length(loopAvg),loopAvg);
hold all;
pcoef = polyfit([1:length(loopAvg)]',loopAvg,1);
plot(polyval(pcoef,[1:length(loopAvg)]'));
title(tit); xlabel('Loop Number (9-target set)'); ylabel('Average Reaction Time (ms)');

%% Average RT by delay length
rts = posRTs(:,[1 3]); % select subgroup to plot
tit = strcat(subj,' Average RT by Delay Length');
[~,ind] = sort(posRTs(:,3));
rts = rts(ind,:);
shortD = rts(rts(:,2)==500);
midD = rts(rts(:,2)==1000);
longD = rts(rts(:,2)==1500);
if isempty(midD)
    y1 = [mean(shortD) mean(longD)];
    y2 = [std(shortD) std(longD)];
    subplot(3,2,5);
    H = bar(1:2,y1,'hist'); hold all;
    e = errorbar(1:2,y1,y2,'.');
    set(H,'facecolor',[.3 .3 .3]);
    e.Color = 'k';
    set(gca,'XTickLabel',{'500 ms','1500 ms'});
else
    y1 = [mean(shortD) mean(midD) mean(longD)];
    y2 = [std(shortD) std(midD) std(longD)];
    subplot(3,2,5);
    H = bar(1:3,y1,'hist'); hold all;
    e = errorbar(1:3,y1,y2,'.');
    set(H,'facecolor',[.3 .3 .3]);
    e.Color = 'k';
    set(gca,'XTickLabel',{'500 ms','1000 ms','1500 ms'});
end
title(tit); xlabel('Delay Length');ylabel('Average RT');

%% Number of guesses per 9-target loop
rts = filtRTs(:,1); % select subgroup to plot
tit = strcat(subj,' Number of Guesses per Loop');
maxRT = 150; % set subphysiological ms cutoff
minRT = min(rts); % set cutoff for 'too early' guesses (note filtRTs scales 'too early' by delay time)
sessGuess = [];
goodTrial = 0;
loopGuess = 0;
for i=1:length(rts)
    if rts(i)<maxRT %&& allRTs(i,1) >minRT
        loopGuess = loopGuess+1;
    end
    if rts(i)>0 && rts(i)<1000
        goodTrial = goodTrial+1;
    end
    if goodTrial ==9
        sessGuess = [sessGuess;loopGuess];
        goodTrial = 0;
        loopGuess = 0;
    end
end
subplot(3,2,6);
H = bar(sessGuess,'hist');
set(H,'facecolor',[.3 .3 .3]);
title(tit); xlabel('Loop Number (9-target set)'); ylabel('Number of Guesses');

if saveFigs
    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
    cd(saveFigsDir);
    print(gcf,overallTitle,'-dpng');
end
end