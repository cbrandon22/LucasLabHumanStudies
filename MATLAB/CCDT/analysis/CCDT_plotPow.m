function CCDT_plotPow
% Default plots freq x time for peri-event window for each item in
% comparison for each electrode and a summary plot of all electrodes
% averaged in fRange
% diffplot will plot t-stats for comparison in fRange. individualTs will
% save line plots for each comparison.

subj = 'HUP143';
task = 'CCDT';
sessList   = {'Session_0'};
alignment = 'cue'; %cue,CC,response
comparison = 'delay'; %delay, RTcombined, timeCombined(early/late trials), RTlong, RTshort, timeLong, timeShort
rawPlot = 0; % Not written yet, to plot raw voltages
plotTstats = 0; % make t-stat plot(s)
individualTs = 0; % make individual channel t-stat plots
fRange = [4 8]; %freqency band to collapse
freqName = '4-8 Hz'; %for plot titles
plotFBandTrials = 1; % collapse freq and plot each trial power for each electrode
saveFigs = 1;
configNum = 1;
dirs = le_dirs(task);
powConfig = le_config_calcPow(configNum,task);
saveFigsDir = fullfile(dirs.scratch,'figs',subj,'powPlots',num2str(configNum),strcat(alignment, '-',comparison));

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

if rawPlot
    %gete_wrapper
end

if plotTstats %Calculate t-stats and plot. 
    tStruct = struct();
    for chan=1:length(channels)
        if chan == 1
            [tStruct,config_pow] = le_calcTs(subj, sessList,channels(chan),task,alignment,comparison,fRange);
        else
            tStruct(chan) = le_calcTs(subj, sessList,channels(chan),task,alignment,comparison,fRange);
        end
        disp(chan)
    end
    % Single channel t-stats
    if individualTs
        le_plotTs(tStruct,config_pow,1);
    end
    % all electrode t-stats
    tMat_pow_cat = squeeze((cat(3,tStruct.tMat)))';
    tBins = nanmean(config_pow.timeBins,2)';
    fBins = [config_pow.freqBins];
    fBins = fBins(tStruct(1).fInd); % filter by freq range used for this tStruct
    
    %TMat Power
    cueTxt = 'cue';cueTxtOffset = 20;
    ccTxt = 'color change';ccTxtOffset = 75;
    tit = [subj ' ' freqName ' ' tStruct(1).retLbl1 ' vs ' tStruct(1).retLbl2 ' ',comparison];
    eLbl_list = {tStruct.eLbl};
    swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
    imagesc(1:length(tBins),1:length(eLbl_list),tMat_pow_cat);
    set(gca,'ydir','normal','clim',[-3 3])
    colormap jet; colorbar;
    swagAxes(gca,10,'Time (ms)','Electrodes',{'t(theta long - theta short)';tit});
    yt = [1:1:length(eLbl_list)];
    set(gca,'ytick',yt,'yticklabel',eLbl_list(yt),'yticklabelrotation',0)
    setFreqXTicks(gca,tBins)
    line([0,0],ylim,'color','w');
    line([500,500],ylim,'color','w');
    line([1000,1000],ylim,'color','w');
    line([1500,1500],ylim,'color','y');
    texty = ylim;
    text(0-cueTxtOffset,texty(2)+1,cueTxt);
    text(1500-ccTxtOffset,texty(2)+1,ccTxt);
    return
end

% Load and select events
evFile = fullfile(dirs.events,strcat(subj,'_events'));
load(evFile,'events');
switch comparison
    case 'delay'
        [shortD,longD] = le_events2retain(events,task,alignment,comparison);
    case 'RTcombined'
        sortRt = sort(cell2mat({events.rt}));
        thresh = [sortRt(round(length(sortRt)/3)),sortRt(round(length(sortRt)/3)*2)];%fastest/slowest 3rd
        [fastRT_s,slowRT_s] = le_events2retain(events,task,alignment,'RTshort',thresh);
        [fastRT_l,slowRT_l] = le_events2retain(events,task,alignment,'RTlong',thresh);
        if plotFBandTrials
            rtEvents = find(fastRT+slowRT);
            temp = cell2mat({events(fastRT|slowRT).rt});
            [~,sortedRTs] = sort(cell2mat({events(fastRT|slowRT).rt}));
            sortedEvInd = rtEvents(sortedRTs);
        end
    case 'RTshort'
        sortRt = sort(cell2mat({events.rt}));
        thresh = [sortRt(round(length(sortRt)/3)),sortRt(round(length(sortRt)/3)*2)];%fastest/slowest 3rd
        [fastRT,slowRT] = le_events2retain(events,task,alignment,comparison,thresh);
    case 'RTlong'
        sortRt = sort(cell2mat({events.rt}));
        thresh = [sortRt(round(length(sortRt)/3)),sortRt(round(length(sortRt)/3)*2)];%fastest/slowest 3rd
        [fastRT,slowRT] = le_events2retain(events,task,alignment,comparison,thresh);
    case 'timeShort'
        [earlyTrial,lateTrial] = le_events2retain(events,task,alignment,comparison);
    case 'timeLong'
        [earlyTrial,lateTrial] = le_events2retain(events,task,alignment,comparison);
    case 'timeCombined'
        [earlyTrial_s,lateTrial_s] = le_events2retain(events,task,alignment,'timeShort');
        [earlyTrial_l,lateTrial_l] = le_events2retain(events,task,alignment,'timeLong');
end

fRangeShort_pow_cat = [];
fRangeLong_pow_cat = [];
fRangeEarly_pow_cat = [];
fRangeLate_pow_cat = [];
fRangeFast_pow_cat = [];
fRangeSlow_pow_cat = [];
for chan = 1:length(channels)
    channel = num2str(channels(chan));
    % Load and concatenate power
    powDir = fullfile(dirs.scratch,'POWER',num2str(configNum),subj,strcat(task,'_events'));
    pow_cat = [];
    for i=1:length(sessList)
        sessPowFile = fullfile(powDir,sessList{i},channel);
        load(sessPowFile);
        pow_cat = cat(3,pow_cat,pow);
    end
    
    % Average across trials and within frequency band separately.
    switch comparison
        case 'delay'
            shortD_pow_cat = pow_cat(:,:,shortD);%Select events
            %midD_pow_cat = pow_cat(:,:,midD);
            longD_pow_cat = pow_cat(:,:,longD);
            combD_pow_cat = cat(3,shortD_pow_cat,longD_pow_cat);
            shortD_pow_cat = nanmean(shortD_pow_cat,3);%average events
            %midD_pow_cat = nanmean(midD_pow_cat,3);
            longD_pow_cat = nanmean(longD_pow_cat,3);
            combD_pow_cat = nanmean(combD_pow_cat,3);
            
            freqLbl = [num2str(fRange(1)) ':' num2str(fRange(2))];
            [~,fInd_start] = min(abs(fRange(1) - powConfig.freQ));
            [~,fInd_end] = min(abs(fRange(2) - powConfig.freQ));
            fInd = fInd_start:fInd_end;
            
            freqShort_pow_cat = pow_cat(fInd,:,shortD); %get short delay events, only in fRange
            freqShort_pow_cat = nanmean(freqShort_pow_cat,3);% Average trails
            freqShort_pow_cat = squeeze(nanmean(nanmean(freqShort_pow_cat,3),1));% Average within fRange
            fRangeShort_pow_cat = [fRangeShort_pow_cat; freqShort_pow_cat];

            freqLong_pow_cat = pow_cat(fInd,:,longD); %get long delay events, only in fRange
            freqLong_pow_cat = nanmean(freqLong_pow_cat,3);% Average trails
            freqLong_pow_cat = squeeze(nanmean(nanmean(freqLong_pow_cat,3),1));% Average within fRange
            fRangeLong_pow_cat = [fRangeLong_pow_cat; freqLong_pow_cat];

            % Make plots
            cueTxt = 'cue';cueTxtOffset = 20;
            ccTxt = 'color change';ccTxtOffset = 75;
            %short delay
            fBins = powConfig.freqBins;
            tBins = nanmean(powConfig.timeBins,2)';
            swagFig([1/8 1/4 3/4 1/2]); hold all
            imagesc(tBins,1:length(fBins),shortD_pow_cat);
            set(gca,'ydir','normal','CLim',[-3 3]);
            setFreqYTicks(gca,fBins)
            swagAxes(gca,10,'time (ms)','frequency (Hz)',['ch ' channel])
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
            line([500,500],ylim,'color','y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(500-ccTxtOffset,texty(2)+1,ccTxt);
            
            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_short'],'-dpng');close
            end
           
            
            %long delay
            fBins = powConfig.freqBins;
            tBins = nanmean(powConfig.timeBins,2)';
            swagFig([1/8 1/4 3/4 1/2]); hold all
            imagesc(tBins,1:length(fBins),longD_pow_cat);
            set(gca,'ydir','normal','CLim',[-3 3]);
            setFreqYTicks(gca,fBins)
            swagAxes(gca,10,'time (ms)','frequency (Hz)',['ch ' channel])
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
            line([500,500],ylim,'color','w');
            %line([1000,1000],ylim,'color','w');
            line([1500,1500],ylim,'color','y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(1500-ccTxtOffset,texty(2)+1,ccTxt);
            
            
            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_long'],'-dpng');close
            end
            
            %Combined delay
            cueTxt = 'cue';cueTxtOffset = 25;
            cc1Txt = 'color change 1';ccTxtOffset = 70;
            cc2Txt = 'color change 2';
            fBins = powConfig.freqBins;
            tBins = nanmean(powConfig.timeBins,2)';
            tBins = tBins(1:35);
            combPlot = [combD_pow_cat(:,1:23),longD_pow_cat(:,24:35)];
            swagFig([1/8 1/4 3/4 1/2]); hold all
            imagesc(tBins,1:length(fBins),combPlot);
            set(gca,'ydir','normal','CLim',[-3 3]);
            setFreqYTicks(gca,fBins)
            swagAxes(gca,10,'time (ms)','frequency (Hz)',['ch ' channel])
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
            line([500,500],ylim,'color','w');
            line([1500,1500],ylim,'color','w');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(500-ccTxtOffset,texty(2)+1,cc1Txt);
            text(1500-ccTxtOffset,texty(2)+1,cc2Txt);

            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_delay'],'-dpng');close
            end
        case 'RTshort'
            fast_pow_cat = pow_cat(:,:,fastRT);
            slow_pow_cat = pow_cat(:,:,slowRT);
            fast_pow_cat = nanmean(fast_pow_cat,3);
            slow_pow_cat = nanmean(slow_pow_cat,3);
            
            freqLbl = [num2str(fRange(1)) ':' num2str(fRange(2))];
            [~,fInd_start] = min(abs(fRange(1) - powConfig.freQ));
            [~,fInd_end] = min(abs(fRange(2) - powConfig.freQ));
            fInd = fInd_start:fInd_end;
            
            freqFast_pow_cat = pow_cat(fInd,:,fastRT); %get early trials, only in fRange
            freqFast_pow_cat = nanmean(freqFast_pow_cat,3);% Average trails
            freqFast_pow_cat = squeeze(nanmean(nanmean(freqFast_pow_cat,3),1));% Average within fRange
            fRangeFast_pow_cat = [fRangeFast_pow_cat; freqFast_pow_cat];
            
            freqSlow_pow_cat = pow_cat(fInd,:,slowRT); %get late trials, only in fRange
            freqSlow_pow_cat = nanmean(freqSlow_pow_cat,3);% Average trails
            freqSlow_pow_cat = squeeze(nanmean(nanmean(freqSlow_pow_cat,3),1));% Average within fRange
            fRangeSlow_pow_cat = [fRangeSlow_pow_cat; freqSlow_pow_cat];
            
            %make plots
            if strcmp(alignment,'cue')
                cueTxt = 'cue';cueTxtOffset = 20;
                ccTxt = 'color change';ccTxtOffset = 75;
            elseif strcmp(alignment,'response')
                cueTxt = 'response';cueTxtOffset = 40;
            end
            
            %fast trials
            fBins = powConfig.freqBins;
            tBins = nanmean(powConfig.timeBins,2)';
            swagFig([1/8 1/4 3/4 1/2]); hold all
            imagesc(tBins,1:length(fBins),fast_pow_cat);
            set(gca,'ydir','normal','CLim',[-3 3]);
            setFreqYTicks(gca,fBins)
            swagAxes(gca,10,'time (ms)','frequency (Hz)')
            yl = ylim;
            title(['ch ' channel ' fast'],'Position',[-1000,yl(2)],'fontsize',20)
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
            line([500,500],ylim,'color','y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(500-ccTxtOffset,texty(2)+1,ccTxt);
            
            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_fast'],'-dpng');close
            end
            
            %slow trials
            fBins = powConfig.freqBins;
            tBins = nanmean(powConfig.timeBins,2)';
            swagFig([1/8 1/4 3/4 1/2]); hold all
            imagesc(tBins,1:length(fBins),slow_pow_cat);
            set(gca,'ydir','normal','CLim',[-3 3]);
            setFreqYTicks(gca,fBins)
            swagAxes(gca,10,'time (ms)','frequency (Hz)')
            yl = ylim;
            title(['ch ' channel ' slow'],'Position',[-1000,yl(2)],'fontsize',20)
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
            line([500,500],ylim,'color','y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(500-ccTxtOffset,texty(2)+1,ccTxt);
            
            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_slow'],'-dpng');close
            end
            
            %Fast-slow difference plot
            fBins = powConfig.freqBins;
            tBins = nanmean(powConfig.timeBins,2)';
            swagFig([1/8 1/4 3/4 1/2]); hold all
            imagesc(tBins,1:length(fBins),fast_pow_cat-slow_pow_cat);
            set(gca,'ydir','normal','CLim',[-3 3]);
            setFreqYTicks(gca,fBins)
            swagAxes(gca,10,'time (ms)','frequency (Hz)')
            yl = ylim;
            title(['ch ' channel ' fast-slow'],'Position',[-1000,yl(2)],'fontsize',20)
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
            line([500,500],ylim,'color','y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(500-ccTxtOffset,texty(2)+1,ccTxt);
            caxis([-0.75 0.75])
            
            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_diff'],'-dpng');close
            end
            
        case 'RTlong'
            fast_pow_cat = pow_cat(:,:,fastRT);
            slow_pow_cat = pow_cat(:,:,slowRT);
            fast_pow_cat = nanmean(fast_pow_cat,3);
            slow_pow_cat = nanmean(slow_pow_cat,3);
            
            freqLbl = [num2str(fRange(1)) ':' num2str(fRange(2))];
            [~,fInd_start] = min(abs(fRange(1) - powConfig.freQ));
            [~,fInd_end] = min(abs(fRange(2) - powConfig.freQ));
            fInd = fInd_start:fInd_end;
            
            freqFast_pow_cat = pow_cat(fInd,:,fastRT); %get early trials, only in fRange
            freqFast_pow_cat = nanmean(freqFast_pow_cat,3);% Average trails
            freqFast_pow_cat = squeeze(nanmean(nanmean(freqFast_pow_cat,3),1));% Average within fRange
            fRangeFast_pow_cat = [fRangeFast_pow_cat; freqFast_pow_cat];
            
            freqSlow_pow_cat = pow_cat(fInd,:,slowRT); %get late trials, only in fRange
            freqSlow_pow_cat = nanmean(freqSlow_pow_cat,3);% Average trails
            freqSlow_pow_cat = squeeze(nanmean(nanmean(freqSlow_pow_cat,3),1));% Average within fRange
            fRangeSlow_pow_cat = [fRangeSlow_pow_cat; freqSlow_pow_cat];
            
            %make plots
            if strcmp(alignment,'cue')
                cueTxt = 'cue';cueTxtOffset = 20;
                ccTxt = 'color change';ccTxtOffset = 75;
            elseif strcmp(alignment,'response')
                cueTxt = 'response';cueTxtOffset = 40;
            end
            
            %fast trials
            fBins = powConfig.freqBins;
            tBins = nanmean(powConfig.timeBins,2)';
            swagFig([1/8 1/4 3/4 1/2]); hold all
            imagesc(tBins,1:length(fBins),fast_pow_cat);
            set(gca,'ydir','normal','CLim',[-3 3]);
            setFreqYTicks(gca,fBins)
            swagAxes(gca,10,'time (ms)','frequency (Hz)')
            yl = ylim;
            title(['ch ' channel ' fast'],'Position',[-1000,yl(2)],'fontsize',20)
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
            line([500,500],ylim,'color','w');
            line([1500,1500],ylim,'color','y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(1500-ccTxtOffset,texty(2)+1,ccTxt);
            
            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_fast'],'-dpng');close
            end
            
            %slow trials
            fBins = powConfig.freqBins;
            tBins = nanmean(powConfig.timeBins,2)';
            swagFig([1/8 1/4 3/4 1/2]); hold all
            imagesc(tBins,1:length(fBins),slow_pow_cat);
            set(gca,'ydir','normal','CLim',[-3 3]);
            setFreqYTicks(gca,fBins)
            swagAxes(gca,10,'time (ms)','frequency (Hz)')
            yl = ylim;
            title(['ch ' channel ' slow'],'Position',[-1000,yl(2)],'fontsize',20)
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
            line([500,500],ylim,'color','w');
            line([1500,1500],ylim,'color','y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(1500-ccTxtOffset,texty(2)+1,ccTxt);
            
            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_slow'],'-dpng');close
            end
            
            %Fast-slow difference plot
            fBins = powConfig.freqBins;
            tBins = nanmean(powConfig.timeBins,2)';
            swagFig([1/8 1/4 3/4 1/2]); hold all
            imagesc(tBins,1:length(fBins),fast_pow_cat-slow_pow_cat);
            set(gca,'ydir','normal','CLim',[-3 3]);
            setFreqYTicks(gca,fBins)
            swagAxes(gca,10,'time (ms)','frequency (Hz)')
            yl = ylim;
            title(['ch ' channel ' fast-slow'],'Position',[-1000,yl(2)],'fontsize',20)
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
            line([500,500],ylim,'color','w');
            line([1500,1500],ylim,'color','y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(1500-ccTxtOffset,texty(2)+1,ccTxt);
            caxis([-0.75 0.75])
            
            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_diff'],'-dpng');close
            end
            
        case 'RTcombined'
            if plotFBandTrials
                %select sorted trials
                rtSort_pow_cat = pow_cat(:,:,sortedEvInd);
                %Collapse freqency
                freqLbl = [num2str(fRange(1)) ':' num2str(fRange(2))];
                [~,fInd_start] = min(abs(fRange(1) - powConfig.freQ));
                [~,fInd_end] = min(abs(fRange(2) - powConfig.freQ));
                fInd = fInd_start:fInd_end;
                rtSort_pow_cat = squeeze(nanmean(rtSort_pow_cat(fInd,:,:)))';
                
                %make plot
                cueTxt = 'cue';cueTxtOffset = 20;
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(sortedRTs),rtSort_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                swagAxes(gca,10,'time (ms)','trial',['ch ' channel])
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                scatter(cell2mat({events(sortedEvInd).delay}),1:length(sortedRTs),'w','o');
                
                if saveFigs
                    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                    cd(saveFigsDir);
                    print(gcf,[channel '_rtSort'],'-dpng');close
                end
            else
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
                
                freqLbl = [num2str(fRange(1)) ':' num2str(fRange(2))];
                [~,fInd_start] = min(abs(fRange(1) - powConfig.freQ));
                [~,fInd_end] = min(abs(fRange(2) - powConfig.freQ));
                fInd = fInd_start:fInd_end;
                
                freqFast_pow_cat = fast_pow_cat(fInd,:); %get short delay events, only in fRange
                freqFast_pow_cat = nanmean(freqFast_pow_cat,3);% Average trails
                freqFast_pow_cat = squeeze(nanmean(nanmean(freqFast_pow_cat,3),1));% Average within fRange
                fRangeFast_pow_cat = [fRangeFast_pow_cat; freqFast_pow_cat];
                
                freqSlow_pow_cat = slow_pow_cat(fInd,:); %get long delay events, only in fRange
                freqSlow_pow_cat = nanmean(freqSlow_pow_cat,3);% Average trails
                freqSlow_pow_cat = squeeze(nanmean(nanmean(freqSlow_pow_cat,3),1));% Average within fRange
                fRangeSlow_pow_cat = [fRangeSlow_pow_cat; freqSlow_pow_cat];
                
                %make plots
                cueTxt = 'cue';cueTxtOffset = 20;
                %ccTxt = 'color change';ccTxtOffset = 75;
                
                %fast RT's
                fBins = powConfig.freqBins;
                tBins = nanmean(powConfig.timeBins,2)';
                tBins = tBins(1:35);
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(fBins),fast_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                setFreqYTicks(gca,fBins)
                swagAxes(gca,10,'time (ms)','frequency (Hz)',['ch ' channel ' fast'])
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');
                line([500,500],ylim,'color','w');
                %line([1000,1000],ylim,'color','y');
                line([1500,1500],ylim,'color','w');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                %text(1000-ccTxtOffset,texty(2)+1,ccTxt);
                
                if saveFigs
                    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                    cd(saveFigsDir);
                    print(gcf,[channel '_fast'],'-dpng');close
                end
                
                %slow RTs
                fBins = powConfig.freqBins;
                tBins = nanmean(powConfig.timeBins,2)';
                tBins = tBins(1:35);
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(fBins),slow_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                setFreqYTicks(gca,fBins)
                swagAxes(gca,10,'time (ms)','frequency (Hz)',['ch ' channel ' slow'])
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');
                line([500,500],ylim,'color','w');
                %line([1000,1000],ylim,'color','y');
                line([1500,1500],ylim,'color','w');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                %text(1000-ccTxtOffset,texty(2)+1,ccTxt);
                
                if saveFigs
                    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                    cd(saveFigsDir);
                    print(gcf,[channel '_slow'],'-dpng');close
                end
%                 % Difference plot
%                 fBins = powConfig.freqBins;
%                 tBins = nanmean(powConfig.timeBins,2)';
%                 tBins = tBins(1:35);
%                 swagFig([1/8 1/4 3/4 1/2]); hold all
%                 imagesc(tBins,1:length(fBins),diff_pow_cat);
%                 set(gca,'ydir','normal','CLim',[-3 3]);
%                 setFreqYTicks(gca,fBins)
%                 swagAxes(gca,10,'time (ms)','frequency (Hz)',['ch ' channel ' fast-slow'])
%                 hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
%                 line([0,0],ylim,'color','w');
%                 line([500,500],ylim,'color','w');
%                 %line([1000,1000],ylim,'color','y');
%                 line([1500,1500],ylim,'color','w');
%                 texty = ylim;
%                 text(0-cueTxtOffset,texty(2)+1,cueTxt);
%                 caxis([-0.75 0.75])
%                 %text(1000-ccTxtOffset,texty(2)+1,ccTxt);
%                 
%                 if saveFigs
%                     if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
%                     cd(saveFigsDir);
%                     print(gcf,[channel '_diff'],'-dpng');close
%                 end
            end
        case 'timeShort'
            early_pow_cat = pow_cat(:,:,earlyTrial);
            late_pow_cat = pow_cat(:,:,lateTrial);
            early_pow_cat = nanmean(early_pow_cat,3);
            late_pow_cat = nanmean(late_pow_cat,3);
            
            freqLbl = [num2str(fRange(1)) ':' num2str(fRange(2))];
            [~,fInd_start] = min(abs(fRange(1) - powConfig.freQ));
            [~,fInd_end] = min(abs(fRange(2) - powConfig.freQ));
            fInd = fInd_start:fInd_end;
            
            freqEarly_pow_cat = pow_cat(fInd,:,earlyTrial); %get early trials, only in fRange
            freqEarly_pow_cat = nanmean(freqEarly_pow_cat,3);% Average trails
            freqEarly_pow_cat = squeeze(nanmean(nanmean(freqEarly_pow_cat,3),1));% Average within fRange
            fRangeEarly_pow_cat = [fRangeEarly_pow_cat; freqEarly_pow_cat];
            
            freqLate_pow_cat = pow_cat(fInd,:,lateTrial); %get late trials, only in fRange
            freqLate_pow_cat = nanmean(freqLate_pow_cat,3);% Average trails
            freqLate_pow_cat = squeeze(nanmean(nanmean(freqLate_pow_cat,3),1));% Average within fRange
            fRangeLate_pow_cat = [fRangeLate_pow_cat; freqLate_pow_cat];
            
            %make plots
            if strcmp(alignment,'cue')
                cueTxt = 'cue';cueTxtOffset = 20;
                ccTxt = 'color change';ccTxtOffset = 75;
            elseif strcmp(alignment,'response')
                cueTxt = 'response';cueTxtOffset = 40;
            end
            
            %early trials
            fBins = powConfig.freqBins;
            tBins = nanmean(powConfig.timeBins,2)';
            swagFig([1/8 1/4 3/4 1/2]); hold all
            imagesc(tBins,1:length(fBins),early_pow_cat);
            set(gca,'ydir','normal','CLim',[-3 3]);
            setFreqYTicks(gca,fBins)
            swagAxes(gca,10,'time (ms)','frequency (Hz)')
            yl = ylim;
            title(['ch ' channel ' early'],'Position',[-1000,yl(2)],'fontsize',20)
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
            line([500,500],ylim,'color','y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(500-ccTxtOffset,texty(2)+1,ccTxt);
            
            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_early'],'-dpng');close
            end
            
            %Late trials
            fBins = powConfig.freqBins;
            tBins = nanmean(powConfig.timeBins,2)';
            swagFig([1/8 1/4 3/4 1/2]); hold all
            imagesc(tBins,1:length(fBins),late_pow_cat);
            set(gca,'ydir','normal','CLim',[-3 3]);
            setFreqYTicks(gca,fBins)
            swagAxes(gca,10,'time (ms)','frequency (Hz)')
            yl = ylim;
            title(['ch ' channel ' late'],'Position',[-1000,yl(2)],'fontsize',20)
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
            line([500,500],ylim,'color','y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(500-ccTxtOffset,texty(2)+1,ccTxt);
            
            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_late'],'-dpng');close
            end
            
            %Late-early difference plot
            fBins = powConfig.freqBins;
            tBins = nanmean(powConfig.timeBins,2)';
            swagFig([1/8 1/4 3/4 1/2]); hold all
            imagesc(tBins,1:length(fBins),late_pow_cat-early_pow_cat);
            set(gca,'ydir','normal','CLim',[-3 3]);
            setFreqYTicks(gca,fBins)
            swagAxes(gca,10,'time (ms)','frequency (Hz)')
            yl = ylim;
            title(['ch ' channel ' late-early'],'Position',[-1000,yl(2)],'fontsize',20)
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
            line([500,500],ylim,'color','y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(500-ccTxtOffset,texty(2)+1,ccTxt);
            caxis([-0.75 0.75])
            
            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_diff'],'-dpng');close
            end
            
        case 'timeLong'
            early_pow_cat = pow_cat(:,:,earlyTrial);
            late_pow_cat = pow_cat(:,:,lateTrial);
            early_pow_cat = nanmean(early_pow_cat,3);
            late_pow_cat = nanmean(late_pow_cat,3);
            
            freqLbl = [num2str(fRange(1)) ':' num2str(fRange(2))];
            [~,fInd_start] = min(abs(fRange(1) - powConfig.freQ));
            [~,fInd_end] = min(abs(fRange(2) - powConfig.freQ));
            fInd = fInd_start:fInd_end;
            
            freqEarly_pow_cat = pow_cat(fInd,:,earlyTrial); %get early trials, only in fRange
            freqEarly_pow_cat = nanmean(freqEarly_pow_cat,3);% Average trails
            freqEarly_pow_cat = squeeze(nanmean(nanmean(freqEarly_pow_cat,3),1));% Average within fRange
            fRangeEarly_pow_cat = [fRangeEarly_pow_cat; freqEarly_pow_cat];
            
            freqLate_pow_cat = pow_cat(fInd,:,lateTrial); %get late trials, only in fRange
            freqLate_pow_cat = nanmean(freqLate_pow_cat,3);% Average trails
            freqLate_pow_cat = squeeze(nanmean(nanmean(freqLate_pow_cat,3),1));% Average within fRange
            fRangeLate_pow_cat = [fRangeLate_pow_cat; freqLate_pow_cat];
            
            %make plots
            if strcmp(alignment,'cue')
                cueTxt = 'cue';cueTxtOffset = 20;
                ccTxt = 'color change';ccTxtOffset = 75;
            elseif strcmp(alignment,'response')
                cueTxt = 'response';cueTxtOffset = 40;
            end
            
            %early trials
            fBins = powConfig.freqBins;
            tBins = nanmean(powConfig.timeBins,2)';
            swagFig([1/8 1/4 3/4 1/2]); hold all
            imagesc(tBins,1:length(fBins),early_pow_cat);
            set(gca,'ydir','normal','CLim',[-3 3]);
            setFreqYTicks(gca,fBins)
            swagAxes(gca,10,'time (ms)','frequency (Hz)')
            yl = ylim;
            title(['ch ' channel ' early'],'Position',[-1000,yl(2)],'fontsize',20)
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
            line([500,500],ylim,'color','w');
            line([1500,1500],ylim,'color','y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(1500-ccTxtOffset,texty(2)+1,ccTxt);
            
            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_early'],'-dpng');close
            end
            
            %Late trials
            fBins = powConfig.freqBins;
            tBins = nanmean(powConfig.timeBins,2)';
            swagFig([1/8 1/4 3/4 1/2]); hold all
            imagesc(tBins,1:length(fBins),late_pow_cat);
            set(gca,'ydir','normal','CLim',[-3 3]);
            setFreqYTicks(gca,fBins)
            swagAxes(gca,10,'time (ms)','frequency (Hz)')
            yl = ylim;
            title(['ch ' channel ' late'],'Position',[-1000,yl(2)],'fontsize',20)
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
            line([500,500],ylim,'color','w');
            line([1500,1500],ylim,'color','y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(1500-ccTxtOffset,texty(2)+1,ccTxt);
            
            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_late'],'-dpng');close
            end
            
            %Late-early difference plot
            fBins = powConfig.freqBins;
            tBins = nanmean(powConfig.timeBins,2)';
            swagFig([1/8 1/4 3/4 1/2]); hold all
            imagesc(tBins,1:length(fBins),late_pow_cat-early_pow_cat);
            set(gca,'ydir','normal','CLim',[-3 3]);
            setFreqYTicks(gca,fBins)
            swagAxes(gca,10,'time (ms)','frequency (Hz)')
            yl = ylim;
            title(['ch ' channel ' late-early'],'Position',[-1000,yl(2)],'fontsize',20)
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
            line([500,500],ylim,'color','w');
            line([1500,1500],ylim,'color','y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(1500-ccTxtOffset,texty(2)+1,ccTxt);
            caxis([-0.75 0.75])
            
            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_diff'],'-dpng');close
            end
            
        case 'timeCombined'
            earlyTrials = or(earlyTrial_s,earlyTrial_l);
            lateTrials = or(lateTrial_s,lateTrial_l);
            early_pow_cat = pow_cat(:,1:23,earlyTrials);
            early_pow_catl = pow_cat(:,24:35,earlyTrial_l);
            late_pow_cat = pow_cat(:,1:23,lateTrials);
            late_pow_catl = pow_cat(:,24:35,lateTrial_l);
            early_pow_cat = nanmean(early_pow_cat,3);
            early_pow_catl = nanmean(early_pow_catl,3);
            late_pow_cat = nanmean(late_pow_cat,3);
            late_pow_catl = nanmean(late_pow_catl,3);
            early_pow_cat = cat(2,early_pow_cat,early_pow_catl);
            late_pow_cat = cat(2,late_pow_cat,late_pow_catl);
            diff_pow_cat = late_pow_cat - early_pow_cat;
            
            freqLbl = [num2str(fRange(1)) ':' num2str(fRange(2))];
            [~,fInd_start] = min(abs(fRange(1) - powConfig.freQ));
            [~,fInd_end] = min(abs(fRange(2) - powConfig.freQ));
            fInd = fInd_start:fInd_end;
            
            freqEarly_pow_cat = early_pow_cat(fInd,:); %get short delay events, only in fRange
            freqEarly_pow_cat = nanmean(freqEarly_pow_cat,3);% Average trails
            freqEarly_pow_cat = squeeze(nanmean(nanmean(freqEarly_pow_cat,3),1));% Average within fRange
            fRangeEarly_pow_cat = [fRangeEarly_pow_cat; freqEarly_pow_cat];

            freqLate_pow_cat = late_pow_cat(fInd,:); %get long delay events, only in fRange
            freqLate_pow_cat = nanmean(freqLate_pow_cat,3);% Average trails
            freqLate_pow_cat = squeeze(nanmean(nanmean(freqLate_pow_cat,3),1));% Average within fRange
            fRangeLate_pow_cat = [fRangeLate_pow_cat; freqLate_pow_cat];
            
            %make plots
            cueTxt = 'cue';cueTxtOffset = 25;
            cc1Txt = 'color change 1';ccTxtOffset = 70;
            cc2Txt = 'color change 2';
            
            %Early trials 
            fBins = powConfig.freqBins;
            tBins = nanmean(powConfig.timeBins,2)';
            tBins = tBins(1:35);
            swagFig([1/8 1/4 3/4 1/2]); hold all
            imagesc(tBins,1:length(fBins),early_pow_cat);
            set(gca,'ydir','normal','CLim',[-3 3]);
            setFreqYTicks(gca,fBins)
            swagAxes(gca,10,'time (ms)','frequency (Hz)',['ch ' channel ' early'])
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
            line([500,500],ylim,'color','w');
            %line([1000,1000],ylim,'color','y');
            line([1500,1500],ylim,'color','w');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(500-ccTxtOffset,texty(2)+1,cc1Txt);
            text(1500-ccTxtOffset,texty(2)+1,cc2Txt);
            
            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_early'],'-dpng');close
            end
            
            %late trials
            fBins = powConfig.freqBins;
            tBins = nanmean(powConfig.timeBins,2)';
            tBins = tBins(1:35);
            swagFig([1/8 1/4 3/4 1/2]); hold all
            imagesc(tBins,1:length(fBins),late_pow_cat);
            set(gca,'ydir','normal','CLim',[-3 3]);
            setFreqYTicks(gca,fBins)
            swagAxes(gca,10,'time (ms)','frequency (Hz)',['ch ' channel ' late'])
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
            line([500,500],ylim,'color','w');
            %line([1000,1000],ylim,'color','y');
            line([1500,1500],ylim,'color','w');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(500-ccTxtOffset,texty(2)+1,cc1Txt);
            text(1500-ccTxtOffset,texty(2)+1,cc2Txt);
            
            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_late'],'-dpng');close
            end
%             % Difference plot
%             fBins = powConfig.freqBins;
%             tBins = nanmean(powConfig.timeBins,2)';
%             tBins = tBins(1:35);
%             swagFig([1/8 1/4 3/4 1/2]); hold all
%             imagesc(tBins,1:length(fBins),diff_pow_cat);
%             set(gca,'ydir','normal','CLim',[-3 3]);
%             setFreqYTicks(gca,fBins)
%             swagAxes(gca,10,'time (ms)','frequency (Hz)',['ch ' channel ' late-early'])
%             hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
%             line([0,0],ylim,'color','w');
%             line([500,500],ylim,'color','w');
%             %line([1000,1000],ylim,'color','y');
%             line([1500,1500],ylim,'color','w');
%             texty = ylim;
%             text(0-cueTxtOffset,texty(2)+1,cueTxt);
%             text(500-ccTxtOffset,texty(2)+1,cc1Txt);
%             text(1500-ccTxtOffset,texty(2)+1,cc2Txt);
%             caxis([-0.75 0.75])
%             
%             if saveFigs
%                 if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
%                 cd(saveFigsDir);
%                 print(gcf,[channel '_diff'],'-dpng');close
%             end
    end
end
% Make plots
cueTxt = 'cue';cueTxtOffset = 20;
ccTxt = 'color change';ccTxtOffset = 75;

%All electrodes in frequency band
%set comparison to plot
switch comparison
    case 'delay'
        freq_pow_cat1 = fRangeShort_pow_cat;
        freq_pow_cat2 = fRangeLong_pow_cat;
        tit1 = [freqName 'Power Short Delay'];
        tit2 = [freqName 'Power Long Delay'];
        fileName1 = [freqName '_short'];
        fileName2 = [freqName '_long'];
    case 'RTshort'
        freq_pow_cat1 = fRangeFast_pow_cat;
        freq_pow_cat2 = fRangeSlow_pow_cat;
        tit1 = [freqName 'Power Fast Trials'];
        tit2 = [freqName 'Power Slow Trials'];
        tit3 = [freqName 'Power Fast - Slow Trials'];
        fileName1 = [freqName '_fast'];
        fileName2 = [freqName '_slow'];
    case 'RTlong'
        freq_pow_cat1 = fRangeFast_pow_cat;
        freq_pow_cat2 = fRangeSlow_pow_cat;
        tit1 = [freqName 'Power Fast Trials'];
        tit2 = [freqName 'Power Slow Trials'];
        tit3 = [freqName 'Power Fast - Slow Trials'];
        fileName1 = [freqName '_fast'];
        fileName2 = [freqName '_slow'];
    case 'RTcombined'
        freq_pow_cat1 = fRangeFast_pow_cat;
        freq_pow_cat2 = fRangeSlow_pow_cat;
        tit1 = [freqName 'Power Fast Trials'];
        tit2 = [freqName 'Power Slow Trials'];
        tit3 = [freqName 'Power Fast - Slow Trials'];
        fileName1 = [freqName '_fast'];
        fileName2 = [freqName '_slow'];
    case 'timeShort'
        freq_pow_cat1 = fRangeEarly_pow_cat;
        freq_pow_cat2 = fRangeLate_pow_cat;
        tit1 = [freqName 'Power Early Trials'];
        tit2 = [freqName 'Power Late Trials'];
        tit3 = [freqName 'Power Late - Early Trials'];
        fileName1 = [freqName '_early'];
        fileName2 = [freqName '_late'];
    case 'timeLong'
        freq_pow_cat1 = fRangeEarly_pow_cat;
        freq_pow_cat2 = fRangeLate_pow_cat;
        tit1 = [freqName 'Power Early Trials'];
        tit2 = [freqName 'Power Late Trials'];
        tit3 = [freqName 'Power Late - Early Trials'];
        fileName1 = [freqName '_early'];
        fileName2 = [freqName '_late'];
    case 'timeCombined'
        freq_pow_cat1 = fRangeEarly_pow_cat;
        freq_pow_cat2 = fRangeLate_pow_cat;
        tit1 = [freqName 'Power Early Trials'];
        tit2 = [freqName 'Power Late Trials'];
        tit3 = [freqName 'Power Late - Early Trials'];
        fileName1 = [freqName '_early'];
        fileName2 = [freqName '_late'];
end

%short delay/early trials
tBins = nanmean(powConfig.timeBins,2)';
swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
imagesc(tBins,1:length(channels),freq_pow_cat1);
set(gca,'ydir','normal','CLim',[-3 3]);
yt = [1:1:length(channels)];
set(gca,'ytick',yt,'yticklabel',channels,'yticklabelrotation',0)
swagAxes(gca,10,'time (ms)','electrode',tit1)
hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
line([0,0],ylim,'color','w');
%line([500,500],ylim,'color','y');% short delay
line([500,500],ylim,'color','w');
line([1000,1000],ylim,'color','w');
line([1500,1500],ylim,'color','y');
texty = ylim;
text(0-cueTxtOffset,texty(2)+1,cueTxt);
text(1500-ccTxtOffset,texty(2)+1,ccTxt);
%text(500-ccTxtOffset,texty(2)+1,ccTxt);% short delay

if saveFigs
    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
    cd(saveFigsDir);
    print(gcf,fileName1,'-dpng');close
end

%long delay/late trials
tBins = nanmean(powConfig.timeBins,2)';
swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
imagesc(tBins,1:length(channels),freq_pow_cat2);
set(gca,'ydir','normal','CLim',[-3 3]);
yt = [1:1:length(channels)];
set(gca,'ytick',yt,'yticklabel',channels,'yticklabelrotation',0)
swagAxes(gca,10,'time (ms)','electrode',tit2)
hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
line([0,0],ylim,'color','w');
%line([500,500],ylim,'color','y');% short delay
line([500,500],ylim,'color','w');
line([1000,1000],ylim,'color','w');
line([1500,1500],ylim,'color','y');
texty = ylim;
text(0-cueTxtOffset,texty(2)+1,cueTxt);
text(1500-ccTxtOffset,texty(2)+1,ccTxt);
%text(500-ccTxtOffset,texty(2)+1,ccTxt);% short delay

if saveFigs
    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
    cd(saveFigsDir);
    print(gcf,fileName2,'-dpng');close
end

% Subtract plots
if ~strcmp(comparison,'delay') && ~strcmp(comparison,'RTcombined')
    diff_pow_cat = freq_pow_cat2 - freq_pow_cat1;
    tBins = nanmean(powConfig.timeBins,2)';
    swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
    imagesc(tBins,1:length(channels),diff_pow_cat);
    set(gca,'ydir','normal','CLim',[-3 3]);
    yt = [1:1:length(channels)];
    set(gca,'ytick',yt,'yticklabel',channels,'yticklabelrotation',0)
    swagAxes(gca,10,'time (ms)','electrode',tit3)
    hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
    line([0,0],ylim,'color','w');
    %line([500,500],ylim,'color','y');% short delay
    line([500,500],ylim,'color','w');
    line([1000,1000],ylim,'color','w');
    line([1500,1500],ylim,'color','y');
    texty = ylim;
    text(0-cueTxtOffset,texty(2)+1,cueTxt);
    text(1500-ccTxtOffset,texty(2)+1,ccTxt);
    caxis([-0.75 0.75])
    %text(500-ccTxtOffset,texty(2)+1,ccTxt);% short delay
    
    if saveFigs
        if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
        cd(saveFigsDir);
        print(gcf,[freqName '_' comparison '_diff'],'-dpng');close
    end
end

% tBins = nanmean(powConfig.timeBins,2)';
% swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
% imagesc(tBins,1:length(channels),diff_theta_pow_cat);
% set(gca,'ydir','normal','CLim',[-3 3]);
% yt = [1:1:length(channels)];
% set(gca,'ytick',yt,'yticklabel',channels,'yticklabelrotation',0)
% swagAxes(gca,10,'time (ms)','electrode','Theta Power Short-Long Delay')
% hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
% line([0,0],ylim,'color','w');
% line([500,500],ylim,'color','w');
% line([1000,1000],ylim,'color','w');
% line([1500,1500],ylim,'color','y');
% texty = ylim;
% text(0-cueTxtOffset,texty(2)+1,cueTxt);
% text(1500-ccTxtOffset,texty(2)+1,ccTxt);
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