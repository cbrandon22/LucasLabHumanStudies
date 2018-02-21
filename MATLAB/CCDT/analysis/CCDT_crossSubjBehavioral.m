subjList = {'HUP136','HUP140','HUP142','HUP143'};
%figure('Name','CCDT Behavioral Summary','NumberTitle','off');
subProws = 2;
subPcols = 3;
plotType = 'rtFastSlow';%'rtScatter' 'rtSets' 'rtFastSlow'
threshList = [];
allShort = [];
allLong = [];
fastShort = [];
slowShort = [];
fastLong = [];
slowLong = [];
for ii=1:length(subjList)
    subj = subjList{ii};
    sessNum = []; % Pick specific session(s). Default combines all sessions
    allRTs = CCDTcatBehavioral(subj,sessNum);
    saveFigs = 1;
    dirs = le_dirs('CCDT');
    saveFigsDir = fullfile(dirs.scratch,'figs','allSubj','behavioral');
    
    %% Organize RTs into subgroups
    filtRTs = allRTs((allRTs(:,1)+allRTs(:,3)>300 & allRTs(:,1)~=1000),:);
    posRTs = filtRTs((filtRTs(:,1)>0),:); % Reaction times after CC
    guessRTs = filtRTs((filtRTs(:,1)<150),:); % Early/subphysiological presses
    fastRTs = posRTs((posRTs(:,1)<250),:); % Fast RTs
    slowRTs = posRTs((posRTs(:,1)>750),:); % Slow RTs
    switch plotType
        case 'rtScatter'
            %% Scatter plot
            rts = filtRTs(:,1); % select subgroup to plot
            tit = strcat(subj,' All Reaction Times');
            
            subplot(subProws,subPcols,ii);
            scatter([1:length(rts)],rts);
            hold all;
            pcoef = polyfit([1:length(rts)]',rts,1);
            plot(polyval(pcoef,[1:length(rts)]'));
            plot([0 length(rts)],[250 250],'g');
            plot([0 length(rts)],[150 150],'g--');
            plot([0 length(rts)],[0 0],'k');
            title(tit);xlabel('trial number');ylabel('reaction time (ms)');
            legend('RT','best-fit','Fast','SubPhys','Location','eastoutside')
            
        case 'rtSets'
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
            subplot(subProws,subPcols,ii);
            scatter(1:length(loopAvg),loopAvg);
            hold all;
            pcoef = polyfit([1:length(loopAvg)]',loopAvg,1);
            plot(polyval(pcoef,[1:length(loopAvg)]'));
            title(tit); xlabel('Loop Number (9-target set)'); ylabel('Average Reaction Time (ms)');
        case 'rtFastSlow'
            %% Average RT by delay length
            rts = posRTs(:,[1 3]); % select subgroup to plot
            [~,ind] = sort(posRTs(:,3));
            rts = rts(ind,:);
            shortD = rts(rts(:,2)==500);
            longD = rts(rts(:,2)==1500);
            allShort = [allShort;shortD];
            allLong = [allLong;longD];
            sortRts = sort(shortD);
            sortRtl = sort(longD);
            threshs = [sortRts(round(length(sortRts)/3)),sortRts(round(length(sortRts)/3)*2)];%fastest/slowest 3rd
            threshl = [sortRtl(round(length(sortRtl)/3)),sortRtl(round(length(sortRtl)/3)*2)];
            threshList = [threshList;threshs(1),threshs(2),threshl(1),threshl(2)];
            fastShort = [fastShort;sortRts(sortRts<threshs(1))];
            slowShort = [slowShort;sortRts(sortRts>threshs(2))];
            fastLong = [fastLong;sortRtl(sortRtl<threshl(1))];
            slowLong = [slowLong;sortRtl(sortRtl>threshl(2))];
            
    end
end
yy1 = [mean(fastShort) mean(slowShort) mean(fastLong) mean(slowLong)];
yy2 = [std(fastShort)/sqrt(length(fastShort)) std(slowShort)/sqrt(length(slowShort)) std(fastLong)/sqrt(length(fastLong)) std(slowLong)/sqrt(length(fastShort))];
tit = 'Average Fast/Slow RT by Delay Length';
figure;
H = bar(1:4,yy1,'hist'); hold all;
e = errorbar(1:4,yy1,yy2,'.');
set(H,'facecolor',[.3 .3 .3]);
e.Color = 'k';
set(gca,'XTickLabel',{'fast-500 ms','slow-500 ms','fast-1500 ms','slow-1500 ms'});
title(tit); xlabel('RT Subset');ylabel('Average RT');
sortRts = sort(allShort);
sortRtl = sort(allLong);
threshs = [sortRts(round(length(sortRts)/3)),sortRts(round(length(sortRts)/3)*2)];%fastest/slowest 3rd
threshl = [sortRtl(round(length(sortRtl)/3)),sortRtl(round(length(sortRtl)/3)*2)];

y1 = [mean(allShort) mean(allLong)];
y2 = [std(allShort) std(allLong)];
figure;
H = bar(1:2,y1,'hist'); hold all;
e = errorbar(1:2,y1,y2,'.');
set(H,'facecolor',[.3 .3 .3]);
e.Color = 'k';
set(gca,'XTickLabel',{'500 ms','1500 ms'});
title(tit); xlabel('RT Subset');ylabel('Average RT');



% bins = 0:50:1000
% n=histc(allShort,bins)
% n2=histc(allLong,bins)
% subplot(2,1,1)
% hold on
% H1=bar(bins(1:end-1),n(1:end-1))
% set(H1,'facecolor',[.3 .3 .3]);
% title('500ms delay')
% plot([threshs(1) threshs(1)],ylim,'r','LineWidth',2)
% plot([threshs(2) threshs(2)],ylim,'r','LineWidth',2)
% subplot(2,1,2)
% hold on
% H2=bar(bins(1:end-1),n2(1:end-1))
% set(H2,'facecolor',[.3 .3 .3]);
% title('1500ms delay')
% plot([threshl(1) threshl(1)],ylim,'r','LineWidth',2)
% plot([threshl(2) threshl(2)],ylim,'r','LineWidth',2)

keyboard
if saveFigs
    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
    cd(saveFigsDir);
    print(gcf,['CCDT Behavioral Summary- ' plotType] ,'-dpng');
end