function CCDT_behavioral
% Quick look at behavioral data for CCDT
subj = 'HUP133';
sessNum = [0]; % Pick specific session(s). Default combines all sessions
allRTs = CCDTcatBehavioral(subj,sessNum);
saveFigs = 1;
dirs = le_dirs('CCDT');
saveFigsDir = fullfile(dirs.scratch,'figs/CCDT',subj,'behavioral');

%% Organize RTs into subgroups
if strcmp(subj,'HUP133')
    filtRTs = allRTs((allRTs(:,1)+allRTs(:,3)>300 & allRTs(:,1)~=1),:); % Filter out very early presses & non-response
else
    filtRTs = allRTs((allRTs(:,1)+allRTs(:,3)>300 & allRTs(:,1)~=1000),:);
end
posRTs = filtRTs((filtRTs(:,1)>0),:); % Reaction times after CC
guessRTs = filtRTs((filtRTs(:,1)<150),:); % Early/subphysiological presses
fastRTs = posRTs((posRTs(:,1)<250),:); % Fast RTs
slowRTs = posRTs((posRTs(:,1)>750),:); % Slow RTs

%% Scatter plot
rts = filtRTs(:,1); % select subgroup to plot
tit = strcat(subj,' Filtered Reaction Times');

figure;
scatter([1:length(rts)],rts);
hold all;
pcoef = polyfit([1:length(rts)]',rts,1);
plot(polyval(pcoef,[1:length(rts)]'));
title(tit);xlabel('trial number');ylabel('reaction time (ms)');

if saveFigs
    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
    cd(saveFigsDir);
    print(gcf,'allRTscatter','-dpng');close
end
%% Bar plot of all RTs
rts = posRTs(:,1);
tit = strcat(subj,' Reaction Time Distribution');
figure;
[N,edges] = histcounts(rts);
x = edges(1:end-1) + .5*(unique(diff(edges)));
H = bar(x,N,'hist');
set(H,'facecolor',[.3 .3 .3]);
title(tit);xlabel('reaction time (ms)');ylabel('count');

if saveFigs
    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
    cd(saveFigsDir);
    print(gcf,'allRTbar','-dpng');close
end
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
figure;
scatter(1:length(loopAvg),loopAvg);
hold all;
pcoef = polyfit([1:length(loopAvg)]',loopAvg,1);
plot(polyval(pcoef,[1:length(loopAvg)]'));
title(tit); xlabel('Loop Number (9-target set)'); ylabel('Average Reaction Time (ms)');

if saveFigs
    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
    cd(saveFigsDir);
    print(gcf,'9targAvgScatter','-dpng');close
end
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
figure;
H = bar(sessGuess,'hist');
set(H,'facecolor',[.3 .3 .3]);
title(tit); xlabel('Loop Number (9-target set)'); ylabel('Number of Guesses');

if saveFigs
    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
    cd(saveFigsDir);
    print(gcf,'9targGuesses','-dpng');close
end
%% Fast and slow RTs per 9-target loop
rts = filtRTs(:,1); % select subgroup to plot
tit = strcat(subj,' Fast and Slow RTs per Loop');
fastRT = [0 250]; % Set ranges for fast and slow classifications
slowRT = [250 1000];

sessFast = [];
sessSlow = [];
goodTrial = 0;
loopFast = 0;
loopSlow = 0;
for i=1:length(rts)
    if rts(i)>fastRT(1) && rts(i,1) <fastRT(2)
        loopFast = loopFast+1;
    elseif rts(i)>slowRT(1) && rts(i,1) <slowRT(2)
        loopSlow = loopSlow+1;
    end
    if rts(i)>0 && rts(i)<1000
        goodTrial = goodTrial+1;
    end
    if goodTrial == 9
        sessFast = [sessFast;loopFast];
        sessSlow = [sessSlow;loopSlow];
        goodTrial = 0;
        loopFast = 0;
        loopSlow = 0;
    end
end
figure;
H = bar([sessFast sessSlow],'hist');
set(H(1),'facecolor','g')
set(H(2),'facecolor','r')
title(tit); xlabel('Loop Number (9-target set)');ylabel('Count');
legend({'Fast','Slow'})

if saveFigs
    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
    cd(saveFigsDir);
    print(gcf,'9targFastVsSlow','-dpng');close
end
%% Average RT by delay length
rts = posRTs(:,[1 3]); % select subgroup to plot
tit = strcat(subj,' Average RT by Delay Length');
[~,ind] = sort(posRTs(:,3));
rts = rts(ind,:);
shortD = rts(rts(:,2)==500);
midD = rts(rts(:,2)==1000);
longD = rts(rts(:,2)==1500);
y1 = [mean(shortD) mean(midD) mean(longD)];
y2 = [std(shortD) std(midD) std(longD)];
H = bar(1:3,y1,'hist'); hold all;
e = errorbar(1:3,y1,y2,'.');
set(H,'facecolor',[.3 .3 .3]);
e.Color = 'k';
set(gca,'XTickLabel',{'500 ms','1000 ms','1500 ms'})
title(tit); xlabel('Delay Length');ylabel('Average RT');

if saveFigs
    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
    cd(saveFigsDir);
    print(gcf,'DelayBar','-dpng');close
end

end