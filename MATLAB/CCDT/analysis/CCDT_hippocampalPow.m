function CCDT_hippocampalPow
subjList = {'HUP136','HUP140','HUP142','HUP143','HUP145','HUP146'};
numElec = 3; % number of electrodes to include from hippocampal depths (most medial)
task = 'CCDT';
alignment = 'cue'; %cue,CC,response
comparison = 'RTsplit'; %delay, RT, time(early/late trials)Long or Short delay
saveFigs = 1;
configNum = 1;
powConfig = le_config_calcPow(configNum,task);
dirs = le_dirs('CCDT');
saveFigsDir = fullfile(dirs.scratch,'figs','allSubj','HippoPow');
trialCount = zeros(1,2); %Total trials,long trials
fRange = [4 8];
tRange = [1000 1500];

switch comparison
    case 'delay'
        allDelayPow = [];
        trialCount = zeros(1,2); %Total trials,long trials
    case 'RT'
        allFastPow = [];
        allSlowPow = [];
        allRTdiffPow = [];
        fastTrialCount = zeros(1,2); %Total trials,long trials
        slowTrialCount = zeros(1,2); %Total trials,long trials
    case 'RTsplit'
        longFastPow = [];
        shortFastPow =[];
        longSlowPow = [];
        shortSlowPow = [];
        diffTshort = [];
        diffTlong = [];
        longRTdiffPow = [];
        shortRTdiffPow = [];
        tfDistLongFast = [];
        tfDistLongSlow = [];
    case 'time'
        allEarlyPow = [];
        allLatePow = [];
        allTdiffPow = [];
        timeTall = [];
        timeTlate = [];
        timeTs = [];
        tfDistLongEarly = [];
        tfDistLongLate = [];
end
pow_cat = [];
for i=1:length(subjList)
    subj = subjList{i};
    [~,hSubs] = xlsread(fullfile(dirs.scratch,'CCDTsubjHippLeads.xlsx'),'A1:A20');
    [lia,locB] =ismember(subj,hSubs);
    if ~lia
        error('Subject not found in Hippocampal spreadsheet')
    end
    cells2read = ['B' num2str(locB) ':H' num2str(locB)];
    [~,hLabels] = xlsread(fullfile(dirs.scratch,'CCDTsubjHippLeads.xlsx'),cells2read);
    % open jacksheet
    jacFile = fullfile(dirs.data,'eeg',subj,'docs','jacksheet.txt');
    fid = fopen(jacFile,'r');
    JAC=textscan(fid,'%d%s','delimiter','\t');
    JAC{:,2} = strtrim(JAC{:,2});
    fclose(fid);
    channels = JAC{1};
    excludeChan = {'EKG1','EKG2','DC1'};
    [lia,locb] = ismember(excludeChan,JAC{2});
    if sum(lia)>0
        channels(locb(locb~=0)) = [];
    end
    hipChan =[];
    for j=1:length(JAC{2})
        [thisE,eNum] = strsplit(JAC{2}{j},'\d','DelimiterType','RegularExpression');
        [lia,locb] = ismember(thisE{1},hLabels);
        if lia && str2double(eNum)<=numElec
            hipChan = [hipChan JAC{1}(j)];
        end
    end
    evFile = fullfile(dirs.events,strcat(subj,'_events'));
    load(evFile,'events');
    switch comparison
        case 'delay'
            [shortD,longD] = le_events2retain(events,task,alignment,comparison);
            trialCount(1,1)=sum(shortD+longD);
            trialCount(1,2)=sum(longD);
        case 'RT'
            sortRt = sort(cell2mat({events.rt}));
            thresh = [sortRt(round(length(sortRt)/3)),sortRt(round(length(sortRt)/3)*2)];%fastest/slowest 3rd
            [fastRT_s,slowRT_s] = le_events2retain(events,task,alignment,'RTshort',thresh);
            [fastRT_l,slowRT_l] = le_events2retain(events,task,alignment,'RTlong',thresh);
        case 'RTsplit'
            delays = cell2mat({events.delay})==500;
            delayl = cell2mat({events.delay})==1500;
            allRTs = cell2mat({events.rt});
            sortRts = sort(allRTs(delays));
            sortRtl = sort(allRTs(delayl));
            
%             y1 = [mean(sortRt) mean(sortRtl)];
%             y2 = [std(sortRts) std(sortRtl)];
%             figure;
%             H = bar(1:2,y1,'hist'); hold all;
%             e = errorbar(1:2,y1,y2,'.');
%             set(H,'facecolor',[.3 .3 .3]);
%             e.Color = 'k';
%             set(gca,'XTickLabel',{'500 ms','1500 ms'});
%             rts = sortRts;
%             tit = strcat(subj,' Short Delay Reaction Time Distribution');
%             figure;
%             [N,edges] = histcounts(rts);
%             x = edges(1:end-1) + .5*(unique(diff(edges)));
%             H = bar(x,N,'hist');
%             set(H,'facecolor',[.3 .3 .3]);
%             title(tit);xlabel('reaction time (ms)');ylabel('count');
%             rts = sortRtl;
%             tit = strcat(subj,' Long Delay Reaction Time Distribution');
%             figure;
%             [N,edges] = histcounts(rts);
%             x = edges(1:end-1) + .5*(unique(diff(edges)));
%             H = bar(x,N,'hist');
%             set(H,'facecolor',[.3 .3 .3]);
%             title(tit);xlabel('reaction time (ms)');ylabel('count');
            
            threshs = [sortRts(round(length(sortRts)/3)),sortRts(round(length(sortRts)/3)*2)];%fastest/slowest 3rd
            threshl = [sortRtl(round(length(sortRtl)/3)),sortRtl(round(length(sortRtl)/3)*2)];%fastest/slowest 3rd
            [fastRT_s,slowRT_s] = le_events2retain(events,task,alignment,'RTshort',threshs);
            [fastRT_l,slowRT_l] = le_events2retain(events,task,alignment,'RTlong',threshl);
        case 'time'
            [earlyTrial_s,lateTrial_s] = le_events2retain(events,task,alignment,'timeShort');
            [earlyTrial_l,lateTrial_l] = le_events2retain(events,task,alignment,'timeLong');
    end
    % Load and concatenate power
    powDir = fullfile(dirs.scratch,'POWER',num2str(configNum),subj,strcat(task,'_events'));
    for j=1:length(hipChan)
        sessPowFile = fullfile(powDir,'Session_0',num2str(hipChan(j)));
        load(sessPowFile);
        pow_cat = pow;
        % Average across trials and within frequency band separately.
        freqLbl = [num2str(fRange(1)) ':' num2str(fRange(2))];
        [~,fInd_start] = min(abs(fRange(1) - powConfig.freQ));
        [~,fInd_end] = min(abs(fRange(2) - powConfig.freQ));
        fInd = fInd_start:fInd_end;
        [~,tInd_start] = ismember(tRange(1),round(powConfig.timeBins(:,1),-2));
        [~,tInd_end] = ismember(tRange(2),round(powConfig.timeBins(:,2),-2));
        tInd = tInd_start:tInd_end;
        
        switch comparison
            case 'delay'
                combD = or(shortD,longD);
                combD_pow_cat = pow_cat(:,1:23,combD);
                longD_pow_cat = pow_cat(:,24:35,longD);
                combD_pow_cat = nanmean(combD_pow_cat,3);
                longD_pow_cat = nanmean(longD_pow_cat,3);
                combD_pow_cat = cat(2,combD_pow_cat,longD_pow_cat);
                allDelayPow = cat(3,allDelayPow,combD_pow_cat);
                
            case 'RT'
                fastRTs = or(fastRT_s,fastRT_l);
                slowRTs = or(slowRT_s,slowRT_l);
                fast_pow_cat = pow_cat(:,1:23,fastRTs);
                fast_pow_catl = pow_cat(:,24:35,fastRT_l);
                slow_pow_cat = pow_cat(:,1:23,slowRTs);
                slow_pow_catl = pow_cat(:,24:35,slowRT_l);
                fast_pow_cat = nanmean(fast_pow_cat,3);
                fast_pow_catl = nanmean(fast_pow_catl,3);
                slow_pow_cat = nanmean(slow_pow_cat,3);
                slow_pow_catl = nanmean(slow_pow_catl,3);
                fast_pow_cat = cat(2,fast_pow_cat,fast_pow_catl);
                slow_pow_cat = cat(2,slow_pow_cat,slow_pow_catl);
                diff_pow_cat = fast_pow_cat - slow_pow_cat;
                allFastPow = cat(3,allFastPow,fast_pow_cat);
                allSlowPow = cat(3,allSlowPow,slow_pow_cat);
                allRTdiffPow = cat(3,allRTdiffPow,diff_pow_cat);
            case 'RTsplit'
                fastRT_pow_cats = pow_cat(:,1:24,fastRT_s);
                fastRT_pow_catl = pow_cat(:,1:35,fastRT_l);
                slowRT_pow_cats = pow_cat(:,1:24,slowRT_s);
                slowRT_pow_catl = pow_cat(:,1:35,slowRT_l);
                [~,p,~,stats] = ttest2(fastRT_pow_cats,slowRT_pow_cats,'dim',3);
                diffTshort = cat(3,diffTshort,stats.tstat);
                [~,p,~,stats] = ttest2(fastRT_pow_catl,slowRT_pow_catl,'dim',3);
                diffTlong = cat(3,diffTlong,stats.tstat);
                fastRT_pow_cats = nanmean(fastRT_pow_cats,3);
                fastRT_pow_catl = nanmean(fastRT_pow_catl,3);
                slowRT_pow_cats = nanmean(slowRT_pow_cats,3);
                slowRT_pow_catl = nanmean(slowRT_pow_catl,3);
                diff_pow_cats = fastRT_pow_cats - slowRT_pow_cats;
                diff_pow_catl = fastRT_pow_catl - slowRT_pow_catl;
                longFastPow = cat(3,longFastPow,fastRT_pow_catl);
                shortFastPow =cat(3,shortFastPow,fastRT_pow_cats);
                longSlowPow = cat(3,longSlowPow,slowRT_pow_catl);
                shortSlowPow = cat(3,shortSlowPow,slowRT_pow_cats);
                longRTdiffPow = cat(3,longRTdiffPow ,diff_pow_catl);
                shortRTdiffPow = cat(3,shortRTdiffPow,diff_pow_cats);
                tfDistLongFast = [tfDistLongFast; squeeze(nanmean(nanmean(pow_cat(fInd,tInd,fastRT_l),1),2))];
                tfDistLongSlow = [tfDistLongSlow; squeeze(nanmean(nanmean(pow_cat(fInd,tInd,slowRT_l),1),2))];
            case 'time'
                earlyTrials = or(earlyTrial_s,earlyTrial_l);
                lateTrials = or(lateTrial_s,lateTrial_l);
                early_pow_cat = pow_cat(:,1:23,earlyTrials);
                early_pow_catl = pow_cat(:,24:35,earlyTrial_l);
                late_pow_cat = pow_cat(:,1:23,lateTrials);
                late_pow_catl = pow_cat(:,24:35,lateTrial_l);
                tfDistLongEarly = [tfDistLongEarly; squeeze(nanmean(nanmean(pow_cat(fInd,tInd,earlyTrial_l),1),2))];
                tfDistLongLate = [tfDistLongLate; squeeze(nanmean(nanmean(pow_cat(fInd,tInd,lateTrial_l),1),2))];
                if(size(early_pow_cat,3)>size(late_pow_cat,3))
                    early_pow_cat = early_pow_cat(:,:,1:size(late_pow_cat,3));
                elseif(size(early_pow_cat,3)<size(late_pow_cat,3))
                    trimInd = size(late_pow_cat,3)-size(early_pow_cat,3);
                    late_pow_cat = late_pow_cat(:,:,trimInd+1:end);
                end
                if(size(early_pow_catl,3)>size(late_pow_catl,3))
                    early_pow_catl = early_pow_catl(:,:,1:size(late_pow_catl,3));
                elseif(size(early_pow_catl,3)<size(late_pow_catl,3))
                    trimInd = size(late_pow_catl,3)-size(early_pow_catl,3);
                    late_pow_catl = late_pow_catl(:,:,trimInd+1:end);
                end
                [~,p,~,stats] = ttest2(late_pow_cat,early_pow_cat,'dim',3);
                timeTall = stats.tstat;
                [~,p,~,stats] = ttest2(late_pow_catl,early_pow_catl,'dim',3);
                timeTlate = stats.tstat;
                combTs = cat(2,timeTall,timeTlate);
                timeTs = cat(3,timeTs,combTs);
                early_pow_cat = nanmean(early_pow_cat,3);
                early_pow_catl = nanmean(early_pow_catl,3);
                late_pow_cat = nanmean(late_pow_cat,3);
                late_pow_catl = nanmean(late_pow_catl,3);
                early_pow_cat = cat(2,early_pow_cat,early_pow_catl);
                late_pow_cat = cat(2,late_pow_cat,late_pow_catl);
                diff_pow_cat = late_pow_cat - early_pow_cat;
                allEarlyPow = cat(3,allEarlyPow,early_pow_cat);
                allLatePow = cat(3,allLatePow,late_pow_cat);
                allTdiffPow = cat(3,allTdiffPow,diff_pow_cat);
        end
    end
end
% Make plots
switch comparison
    case 'delay'
        cueTxt = 'cue';cueTxtOffset = 25;
        cc1Txt = 'color change 1';ccTxtOffset = 70;
        cc2Txt = 'color change 2';
        fBins = powConfig.freqBins;
        tBins = nanmean(powConfig.timeBins,2)';
        tBins = tBins(1:35);
        swagFig([1/8 1/4 3/4 1/2]); hold all
        imagesc(tBins,1:length(fBins),nanmean(allDelayPow,3));
        set(gca,'ydir','normal','CLim',[-3 3]);
        setFreqYTicks(gca,fBins)
        swagAxes(gca,10,'time (ms)','frequency (Hz)','Hippocampal Power')
        hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
        line([0,0],ylim,'color','w');
        line([500,500],ylim,'color','w');
        line([1500,1500],ylim,'color','w');
        texty = ylim;
        text(0-cueTxtOffset,texty(2)+.3,cueTxt);
        text(500-ccTxtOffset,texty(2)+.3,cc1Txt);
        text(1500-ccTxtOffset,texty(2)+.3,cc2Txt);
        keyboard
        
        if saveFigs
            if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
            cd(saveFigsDir);
            print(gcf,'all_delay','-dpng');close
        end
    case 'RT'
        cueTxt = 'cue';cueTxtOffset = 25;
        cc1Txt = 'color change 1';ccTxtOffset = 70;
        cc2Txt = 'color change 2';
        
        %fast RT's
        fBins = powConfig.freqBins;
        tBins = nanmean(powConfig.timeBins,2)';
        tBins = tBins(1:35);
        swagFig([1/8 1/4 3/4 1/2]); hold all
        imagesc(tBins,1:length(fBins),nanmean(allFastPow,3));
        set(gca,'ydir','normal','CLim',[-3 3]);
        setFreqYTicks(gca,fBins)
        swagAxes(gca,10,'time (ms)','frequency (Hz)','Hipp: Fast Trials')
        hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
        line([0,0],ylim,'color','w');
        line([500,500],ylim,'color','w');
        line([1500,1500],ylim,'color','w');
        texty = ylim;
        text(0-cueTxtOffset,texty(2)+.3,cueTxt);
        text(500-ccTxtOffset,texty(2)+.3,cc1Txt);
        text(1500-ccTxtOffset,texty(2)+.3,cc2Txt);
        
        keyboard
        if saveFigs
            if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
            cd(saveFigsDir);
            print(gcf,'all_fast','-dpng');close
        end
        
        %slow RTs
        fBins = powConfig.freqBins;
        tBins = nanmean(powConfig.timeBins,2)';
        tBins = tBins(1:35);
        swagFig([1/8 1/4 3/4 1/2]); hold all
        imagesc(tBins,1:length(fBins),nanmean(allSlowPow,3));
        set(gca,'ydir','normal','CLim',[-3 3]);
        setFreqYTicks(gca,fBins)
        swagAxes(gca,10,'time (ms)','frequency (Hz)','Hipp: Slow Trials')
        hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
        line([0,0],ylim,'color','w');
        line([500,500],ylim,'color','w');
        line([1500,1500],ylim,'color','w');
        texty = ylim;
        text(0-cueTxtOffset,texty(2)+1,cueTxt);
        text(500-ccTxtOffset,texty(2)+.3,cc1Txt);
        text(1500-ccTxtOffset,texty(2)+.3,cc2Txt);
        
        keyboard
        if saveFigs
            if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
            cd(saveFigsDir);
            print(gcf,'all_slow','-dpng');close
        end
        % Difference plot
        fBins = powConfig.freqBins;
        tBins = nanmean(powConfig.timeBins,2)';
        tBins = tBins(1:35);
        swagFig([1/8 1/4 3/4 1/2]); hold all
        imagesc(tBins,1:length(fBins),nanmean(allRTdiffPow,3));
        set(gca,'ydir','normal','CLim',[-3 3]);
        setFreqYTicks(gca,fBins)
        swagAxes(gca,10,'time (ms)','frequency (Hz)','Hipp: fast-slow')
        hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
        line([0,0],ylim,'color','w');
        line([500,500],ylim,'color','w');
        line([1500,1500],ylim,'color','w');
        texty = ylim;
        text(0-cueTxtOffset,texty(2)+1,cueTxt);
        text(500-ccTxtOffset,texty(2)+.3,cc1Txt);
        text(1500-ccTxtOffset,texty(2)+.3,cc2Txt);
        caxis([-0.75 0.75])
        
        keyboard
        if saveFigs
            if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
            cd(saveFigsDir);
            print(gcf,'RT_diff','-dpng');close
        end
    case 'RTsplit'
        cueTxt = 'cue';cueTxtOffset = 25;
        cc1Txt = 'color change';ccTxtOffset = 70;
        cc2Txt = 'color change';
        
        %fast RT's long delay
        fBins = powConfig.freqBins;
        tBins = nanmean(powConfig.timeBins,2)';
        tBins = tBins(1:35);
        swagFig([1/8 1/4 3/4 1/2]); hold all
        imagesc(tBins,1:length(fBins),nanmean(longFastPow,3));
        set(gca,'ydir','normal','CLim',[-3 3]);
        setFreqYTicks(gca,fBins)
        swagAxes(gca,10,'time (ms)','frequency (Hz)','Hipp: Fast Trials, long delay')
        hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
        line([0,0],ylim,'color','w');
        line([500,500],ylim,'color','w');
        line([1500,1500],ylim,'color','w');
        texty = ylim;
        text(0-cueTxtOffset,texty(2)+.3,cueTxt);
        text(1500-ccTxtOffset,texty(2)+.3,cc2Txt);
        
        keyboard
        if saveFigs
            if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
            cd(saveFigsDir);
            print(gcf,'long_fast','-dpng');close
        end
        
        %fast RT's short delay
        fBins = powConfig.freqBins;
        tBins = nanmean(powConfig.timeBins,2)';
        tBins = tBins(1:24);
        swagFig([1/8 1/4 3/4 1/2]); hold all
        imagesc(tBins,1:length(fBins),nanmean(shortFastPow,3));
        set(gca,'ydir','normal','CLim',[-3 3]);
        setFreqYTicks(gca,fBins)
        swagAxes(gca,10,'time (ms)','frequency (Hz)','Hipp: Fast Trials, short delay')
        hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
        line([0,0],ylim,'color','w');
        line([500,500],ylim,'color','w');
        texty = ylim;
        text(0-cueTxtOffset,texty(2)+.3,cueTxt);
        text(500-ccTxtOffset,texty(2)+.3,cc1Txt);
        
        keyboard
        if saveFigs
            if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
            cd(saveFigsDir);
            print(gcf,'short_fast','-dpng');close
        end
        
        %slow RT's long delay
        fBins = powConfig.freqBins;
        tBins = nanmean(powConfig.timeBins,2)';
        tBins = tBins(1:35);
        swagFig([1/8 1/4 3/4 1/2]); hold all
        imagesc(tBins,1:length(fBins),nanmean(longSlowPow,3));
        set(gca,'ydir','normal','CLim',[-3 3]);
        setFreqYTicks(gca,fBins)
        swagAxes(gca,10,'time (ms)','frequency (Hz)','Hipp: Slow Trials, long delay')
        hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
        line([0,0],ylim,'color','w');
        line([500,500],ylim,'color','w');
        line([1500,1500],ylim,'color','w');
        texty = ylim;
        text(0-cueTxtOffset,texty(2)+.3,cueTxt);
        text(1500-ccTxtOffset,texty(2)+.3,cc2Txt);
        
        keyboard
        if saveFigs
            if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
            cd(saveFigsDir);
            print(gcf,'long_slow','-dpng');close
        end
        
        %slow RT's short delay
        fBins = powConfig.freqBins;
        tBins = nanmean(powConfig.timeBins,2)';
        tBins = tBins(1:24);
        swagFig([1/8 1/4 3/4 1/2]); hold all
        imagesc(tBins,1:length(fBins),nanmean(shortSlowPow,3));
        set(gca,'ydir','normal','CLim',[-3 3]);
        setFreqYTicks(gca,fBins)
        swagAxes(gca,10,'time (ms)','frequency (Hz)','Hipp: Slow Trials, short delay')
        hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
        line([0,0],ylim,'color','w');
        line([500,500],ylim,'color','w');
        texty = ylim;
        text(0-cueTxtOffset,texty(2)+.3,cueTxt);
        text(500-ccTxtOffset,texty(2)+.3,cc1Txt);
        
        keyboard
        if saveFigs
            if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
            cd(saveFigsDir);
            print(gcf,'short_slow','-dpng');close
        end
        
        %Difference plot long delay
        fBins = powConfig.freqBins;
        tBins = nanmean(powConfig.timeBins,2)';
        tBins = tBins(1:35);
        swagFig([1/8 1/4 3/4 1/2]); hold all
        imagesc(tBins,1:length(fBins),nanmean(diffTlong,3));
        %imagesc(tBins,1:length(fBins),nanmean(longRTdiffPow,3));
        set(gca,'ydir','normal','CLim',[-3 3]);
        setFreqYTicks(gca,fBins)
        swagAxes(gca,10,'time (ms)','frequency (Hz)','Hipp: fast-slow, long delay')
        hc = colorbar; set(get(hc,'YLabel'),'String','t-stat'); colormap(jet);
        line([0,0],ylim,'color','w');
        line([500,500],ylim,'color','w');
        line([1500,1500],ylim,'color','w');
        texty = ylim;
        text(0-cueTxtOffset,texty(2)+.3,cueTxt);
        text(1500-ccTxtOffset,texty(2)+.3,cc2Txt);
        
        keyboard
        if saveFigs
            if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
            cd(saveFigsDir);
            print(gcf,'long_RT_diff','-dpng');close
        end
        
        %Difference plot short delay
        fBins = powConfig.freqBins;
        tBins = nanmean(powConfig.timeBins,2)';
        tBins = tBins(1:24);
        swagFig([1/8 1/4 3/4 1/2]); hold all
        imagesc(tBins,1:length(fBins),nanmean(diffTshort,3));
        %imagesc(tBins,1:length(fBins),nanmean(shortRTdiffPow,3));
        set(gca,'ydir','normal','CLim',[-3 3]);
        setFreqYTicks(gca,fBins)
        swagAxes(gca,10,'time (ms)','frequency (Hz)','Hipp: fast-slow, short delay')
        hc = colorbar; set(get(hc,'YLabel'),'String','t-stat'); colormap(jet);
        line([0,0],ylim,'color','w');
        line([500,500],ylim,'color','w');
        texty = ylim;
        text(0-cueTxtOffset,texty(2)+.3,cueTxt);
        text(500-ccTxtOffset,texty(2)+.3,cc1Txt);
        
        keyboard
        if saveFigs
            if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
            cd(saveFigsDir);
            print(gcf,'short_RT_diff','-dpng');close
        end
        
        % Collapsed t/f plot
        hold all;
%         plot(sort(tfDistLongFast));
%         plot(sort(tfDistLongSlow));
        [x1,c1]=hist(tfDistLongFast);
        [x2,c2]=hist(tfDistLongSlow);
        plot(c1,(x1/length(tfDistLongFast)));
        plot(c2,(x2/length(tfDistLongSlow)));
        keyboard
        mes(tfDistLongFast,tfDistLongSlow,'hedgesg')
    case 'time'
        cueTxt = 'cue';cueTxtOffset = 25;
        cc1Txt = 'color change 1';ccTxtOffset = 70;
        cc2Txt = 'color change 2';
        
        %early RT's
        fBins = powConfig.freqBins;
        tBins = nanmean(powConfig.timeBins,2)';
        tBins = tBins(1:35);
        swagFig([1/8 1/4 3/4 1/2]); hold all
        imagesc(tBins,1:length(fBins),nanmean(allEarlyPow,3));
        set(gca,'ydir','normal','CLim',[-3 3]);
        setFreqYTicks(gca,fBins)
        swagAxes(gca,10,'time (ms)','frequency (Hz)','Hipp: Early Trials')
        hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
        line([0,0],ylim,'color','w');
        line([500,500],ylim,'color','w');
        line([1500,1500],ylim,'color','w');
        texty = ylim;
        text(0-cueTxtOffset,texty(2)+.3,cueTxt);
        text(500-ccTxtOffset,texty(2)+.3,cc1Txt);
        text(1500-ccTxtOffset,texty(2)+.3,cc2Txt);
        
        keyboard
        if saveFigs
            if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
            cd(saveFigsDir);
            print(gcf,'early','-dpng');close
        end
        
        %late RTs
        fBins = powConfig.freqBins;
        tBins = nanmean(powConfig.timeBins,2)';
        tBins = tBins(1:35);
        swagFig([1/8 1/4 3/4 1/2]); hold all
        imagesc(tBins,1:length(fBins),nanmean(allLatePow,3));
        set(gca,'ydir','normal','CLim',[-3 3]);
        setFreqYTicks(gca,fBins)
        swagAxes(gca,10,'time (ms)','frequency (Hz)','Hipp: Late Trials')
        hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
        line([0,0],ylim,'color','w');
        line([500,500],ylim,'color','w');
        line([1500,1500],ylim,'color','w');
        texty = ylim;
        text(0-cueTxtOffset,texty(2)+.3,cueTxt);
        text(500-ccTxtOffset,texty(2)+.3,cc1Txt);
        text(1500-ccTxtOffset,texty(2)+.3,cc2Txt);
        
        keyboard
        if saveFigs
            if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
            cd(saveFigsDir);
            print(gcf,'late','-dpng');close
        end
        % Difference plot
        fBins = powConfig.freqBins;
        tBins = nanmean(powConfig.timeBins,2)';
        tBins = tBins(1:35);
        swagFig([1/8 1/4 3/4 1/2]); hold all
        imagesc(tBins,1:length(fBins),nanmean(timeTs(:,:,19:end),3));
        %imagesc(tBins,1:length(fBins),nanmean(allTdiffPow,3));
        set(gca,'ydir','normal','CLim',[-3 3]);
        setFreqYTicks(gca,fBins)
        swagAxes(gca,10,'time (ms)','frequency (Hz)','Hipp: late-early')
        hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
        line([0,0],ylim,'color','w');
        line([500,500],ylim,'color','w');
        line([1500,1500],ylim,'color','w');
        texty = ylim;
        text(0-cueTxtOffset,texty(2)+.3,cueTxt);
        text(500-ccTxtOffset,texty(2)+.3,cc1Txt);
        text(1500-ccTxtOffset,texty(2)+.3,cc2Txt);
        
        keyboard
        if saveFigs
            if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
            cd(saveFigsDir);
            print(gcf,'time_diff','-dpng');close
        end
        
        % Collapsed t/f plot
        hold all;
%         plot(sort(tfDistLongEarly));
%         plot(sort(tfDistLongLate));
        [x1,c1]=hist(tfDistLongEarly);
        [x2,c2]=hist(tfDistLongLate);
        h = plot(c1,(x1/length(tfDistLongEarly)));
        plot(c2,(x2/length(tfDistLongLate)));
        keyboard;
end
end

function [] = setFreqYTicks(ax,fBins)
yt = [[1:5:length(fBins)] length(fBins)];
if any(diff(yt)==0)
    yt = [1:5:length(fBins)];
end
set(ax,'ytick',yt,...
    'yticklabel',round(fBins(yt)))
end

function [] = setFreqXTicks(ax,fBins)
yt = [[1:5:length(fBins)] length(fBins)];
if any(diff(yt)==0)
    yt = [1:5:length(fBins)];
end
set(ax,'xtick',yt,...
    'xticklabel',round(fBins(yt)))
end