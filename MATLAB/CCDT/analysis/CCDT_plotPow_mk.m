function CCDT_plotPow_mk
% Default plots freq x time for peri-event window for each item in
% comparison for each electrode and a summary plot of all electrodes
% averaged in fRange
% diffplot will plot t-stats for comparison in fRange. individualTs will
% save line plots for each comparison.

subj = 'HUP136';
task = 'CCDT';
sessList   = {'Session_0'};
alignment = 'CC'; %cue,CC,response
comparison = 'RT1000'; %delay, RT, time(early/late trials)Long or Short delay %RTerr
%if using RT then change RT value in le_events2retain to the avg RT
rawPlot = 0; 
plotTstats = 0; % make t-stat plot(s)
individualTs = 0; % make individual channel t-stat plots
fRange = [70, 90]; %freqency band to collapse
freqName = '70, 90 Hz'; %for plot titles
plotFBandTrials = 0; % collapse freq and plot each trial power for each electrode
plotFastSlowRT = 1;
plotErrors = 0;
saveFigs = 1;
configNum = 1;
dirs = le_dirs;
powConfig = le_config_calcPow(configNum,task);
saveFigsDir = fullfile(dirs.scratch,'figs',task,subj,'powPlots',num2str(configNum),strcat(alignment, '-',comparison));
if saveFigs
   figDir = fullfile(dirs.scratch,'figs',task,subj,'powPlots',num2str(configNum),strcat(alignment, '-',comparison), freqName);
   if ~exist(figDir,'dir'),mkdir(figDir);end
end 
% open jacksheet
jacFile = fullfile(dirs.data,'eeg',task,subj,'docs','jacksheet.txt');
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
%     line([500,500],ylim,'color','w');
%     line([1000,1000],ylim,'color','w');
%     line([1500,1500],ylim,'color','y');
    texty = ylim;
    text(0-cueTxtOffset,texty(2)+1,cueTxt);
   % text(1500-ccTxtOffset,texty(2)+1,ccTxt);
                if saveFigs
                    saveFigsDir1 = fullfile(dirs.scratch,'figs',task,subj,'Session_0',strcat(comparison, '-'));
                    if ~exist(saveFigsDir1,'dir'),mkdir(saveFigsDir1);end
                    cd(saveFigsDir1);
                    print(gcf,['thetaLong-short'],'-dpng');close
                end
    return
end

% Load and select events
evFile = fullfile(dirs.events,task,strcat(subj,'_events'));
%evFile = fullfile(dirs.data,'eeg',task,subj, 'behavioral', sessList,'events');

load(evFile, 'events');
%load(evFile);

switch comparison
    case 'delay'
        [shortD,longD] = le_events2retain(events,task,alignment,comparison);
    case 'RT'
        [fastRT,slowRT] = le_events2retain(events,task,alignment,comparison);
        
        if plotFBandTrials
            rtEvents = find(fastRT+slowRT);
            temp = cell2mat({events(fastRT|slowRT).rt});
            [~,sortedRTs] = sort(cell2mat({events(fastRT|slowRT).rt}));
            sortedEvInd = rtEvents(sortedRTs);
        elseif plotFastSlowRT
            [fastRT,slowRT] = le_events2retain(events,task,alignment,comparison);
            rtEvents = find(fastRT+slowRT);
            temp = cell2mat({events(fastRT|slowRT).rt});
            [~,sortedRTs] = sort(cell2mat({events(fastRT|slowRT).rt}));
            sortedEvInd = rtEvents(sortedRTs);
        end 
    case 'RT1000'
        if plotFastSlowRT
            [fastRT,slowRT] = le_events2retain(events,task,alignment,comparison);
            rtEvents = find(fastRT+slowRT);
            temp = cell2mat({events(fastRT|slowRT).rt});
            [~,sortedRTs] = sort(cell2mat({events(fastRT|slowRT).rt}));
            sortedEvInd = rtEvents(sortedRTs);
        end
    case 'RTerr'
        if plotErrors
            [errRT,goodRT] = le_events2retain(events,task,alignment,comparison);
            rtEvents = find(errRT);
            temp = cell2mat({events(errRT).rt});
            [~,sortedRTs] = sort(cell2mat({events(errRT).rt}));
            sortedEvInd = rtEvents(sortedRTs);
        end
        
             
    case 'RT1500'
        [fastRT,slowRT] = le_events2retain(events,task,alignment,comparison);
        if plotFBandTrials
            rtEvents = find(fastRT+slowRT);
            temp = cell2mat({events(fastRT|slowRT).rt});
            [~,sortedRTs] = sort(cell2mat({events(fastRT|slowRT).rt}));
            sortedEvInd = rtEvents(sortedRTs);
        
        elseif plotFastSlowRT
            [fastRT,slowRT] = le_events2retain(events,task,alignment,comparison);
            rtEvents = find(fastRT+slowRT);
            temp = cell2mat({events(fastRT|slowRT).rt});
            [~,sortedRTs] = sort(cell2mat({events(fastRT|slowRT).rt}));
            sortedEvInd = rtEvents(sortedRTs);
        end 
    case 'RTerr500'
       if plotErrors
            [errorRT,goodRT] = le_events2retain(events,task,alignment,comparison);
            rtEvents = find(errorRT);
            temp = cell2mat({events(errorRT).rt});
            [~,sortedRTs] = sort(cell2mat({events(errorRT).rt}));
            sortedEvInd = rtEvents(sortedRTs);
        end
        
    case 'RT500'
        [fastRT,slowRT] = le_events2retain(events,task,alignment,comparison);
        if plotFBandTrials
            rtEvents = find(fastRT+slowRT);
            temp = cell2mat({events(fastRT|slowRT).rt});
            [~,sortedRTs] = sort(cell2mat({events(fastRT|slowRT).rt}));
            sortedEvInd = rtEvents(sortedRTs);
        elseif plotFastSlowRT
            [fastRT,slowRT] = le_events2retain(events,task,alignment,comparison);
            rtEvents = find(fastRT+slowRT);
            temp = cell2mat({events(fastRT|slowRT).rt});
            [~,sortedRTs] = sort(cell2mat({events(fastRT|slowRT).rt}));
            sortedEvInd = rtEvents(sortedRTs);
        end 
    case 'RTerr1500'
        if plotErrors
            [errorRT,goodRT] = le_events2retain(events,task,alignment,comparison);
            rtEvents = find(errorRT);
            temp = cell2mat({events(errorRT).rt});
            [~,sortedRTs] = sort(cell2mat({events(errorRT).rt}));
            sortedEvInd = rtEvents(sortedRTs);
        
        end 
    case 'timeShort'
        [earlyTrial,lateTrial] = le_events2retain(events,task,alignment,comparison);
    case 'timeLong'
        [earlyTrial,lateTrial] = le_events2retain(events,task,alignment,comparison);
end

fRangeShort_pow_cat = [];
fRangeLong_pow_cat = [];
fRangeEarly_pow_cat = [];
fRangeLate_pow_cat = [];
fRangeSlow_pow_cat = [];
fRangeFast_pow_cat = [];

for chan = 1:length(channels)
    channel = num2str(channels(chan));
    % Load and concatenate power
    powDir = fullfile(dirs.scratch,'POWER',num2str(configNum),subj,strcat(task,'_events'));
    pow_cat = [];
    sessPowFile = [];
    pow = [];
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
            shortD_pow_cat = nanmean(shortD_pow_cat,3);%average events
            %midD_pow_cat = nanmean(midD_pow_cat,3);
            longD_pow_cat = nanmean(longD_pow_cat,3);
            
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
           
            if strcmp (alignment, 'cue')
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
            swagAxes(gca,10,'time (ms)','frequency (Hz)')
            yl = ylim;
            title(['ch ' channel ' short delay'],'Position',[-1000,yl(2)],'fontsize',20)
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
            
            %mid delay
%             fBins = powConfig.freqBins;
%             tBins = nanmean(powConfig.timeBins,2)';
%             swagFig([1/8 1/4 3/4 1/2]); hold all
%             imagesc(tBins,1:length(fBins),midD_pow_cat);
%             set(gca,'ydir','normal','CLim',[-3 3]);
%             setFreqYTicks(gca,fBins)
%             swagAxes(gca,10,'time (ms)','frequency (Hz)',['ch ' channel])
%             hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
%             line([0,0],ylim,'color','w');
%             line([500,500],ylim,'color','w');
%             line([1000,1000],ylim,'color','y');
%             texty = ylim;
%             text(0-cueTxtOffset,texty(2)+1,cueTxt);
%             text(1000-ccTxtOffset,texty(2)+1,ccTxt);
%             
%             
%             if saveFigs
%                 if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
%                 cd(saveFigsDir);
%                 print(gcf,[channel '_mid'],'-dpng');close
%             end
            
            %long delay
            fBins = powConfig.freqBins;
            tBins = nanmean(powConfig.timeBins,2)';
            swagFig([1/8 1/4 3/4 1/2]); hold all
            imagesc(tBins,1:length(fBins),longD_pow_cat);
            set(gca,'ydir','normal','CLim',[-3 3]);
            setFreqYTicks(gca,fBins)
            swagAxes(gca,10,'time (ms)','frequency (Hz)')
            yl = ylim;
            title(['ch ' channel ' long delay'],'Position',[-1000,yl(2)],'fontsize',20)            
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
            %line([500,500],ylim,'color','w');
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
            
            elseif strcmp(alignment, 'CC')
                % Make plots
            cueTxt = 'color change';cueTxtOffset = 120;
            ccTxt = 'average response';ccTxtOffset = 90;
            %short delay
            fBins = powConfig.freqBins;
            tBins = nanmean(powConfig.timeBins,2)';
            swagFig([1/8 1/4 3/4 1/2]); hold all
            imagesc(tBins,1:length(fBins),shortD_pow_cat);
            set(gca,'ydir','normal','CLim',[-3 3]);
            setFreqYTicks(gca,fBins)
            swagAxes(gca,10,'time (ms)','frequency (Hz)')
            yl = ylim;
            title(['ch ' channel ' short delay'],'Position',[-1000,yl(2)],'fontsize',20)            
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
            line([525,525],ylim,'color','w');
            line ([375, 375], ylim, 'color', 'y');
            line ([725, 725], ylim, 'color', 'y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(540-ccTxtOffset,texty(2)+1,ccTxt);
            
            
            
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
            swagAxes(gca,10,'time (ms)','frequency (Hz)')
            yl = ylim;
            title(['ch ' channel ' long delay'],'Position',[-1000,yl(2)],'fontsize',20)            
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
            line([525,525],ylim,'color','w');
            line ([375, 375], ylim, 'color', 'y');
            line ([725, 725], ylim, 'color', 'y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(540-ccTxtOffset,texty(2)+1,ccTxt);
            
            
            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_long'],'-dpng');close
            end
                
            else strcmp(alignment, 'response')
            % Make plots
            cueTxt = 'response';cueTxtOffset = 110;
            %ccTxt = 'average response';ccTxtOffset = 90;
            %short delay
            fBins = powConfig.freqBins;
            tBins = nanmean(powConfig.timeBins,2)';
            swagFig([1/8 1/4 3/4 1/2]); hold all
            imagesc(tBins,1:length(fBins),shortD_pow_cat);
            set(gca,'ydir','normal','CLim',[-3 3]);
            setFreqYTicks(gca,fBins)
            swagAxes(gca,10,'time (ms)','frequency (Hz)')
            yl = ylim;
            title(['ch ' channel ' short delay'],'Position',[-1000,yl(2)],'fontsize',20)            
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
%             line([440,440],ylim,'color','w');
%             line ([325, 325], ylim, 'color', 'y');
%             line ([825, 825], ylim, 'color', 'y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            %text(440-ccTxtOffset,texty(2)+1,ccTxt);
            
            
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
            swagAxes(gca,10,'time (ms)','frequency (Hz)')
            yl = ylim;
            title(['ch ' channel ' long delay'],'Position',[-1000,yl(2)],'fontsize',20)            
            hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([0,0],ylim,'color','w');
%             line([440,440],ylim,'color','w');
%             line ([325, 325], ylim, 'color', 'y');
%             line ([825, 825], ylim, 'color', 'y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            %text(440-ccTxtOffset,texty(2)+1,ccTxt);
            
            
            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_long'],'-dpng');close
            end
             
            end 
            
        case 'RT'
            if plotFBandTrials
                %select sorted trials
                rtSort_pow_cat = pow_cat(:,:,sortedEvInd);
                %Collapse freqency
                freqLbl = [num2str(fRange(1)) ':' num2str(fRange(2))];
                [~,fInd_start] = min(abs(fRange(1) - powConfig.freQ));
                [~,fInd_end] = min(abs(fRange(2) - powConfig.freQ));
                fInd = fInd_start:fInd_end;
                rtSort_pow_cat = squeeze(nanmean(rtSort_pow_cat(fInd,:,:)))';
                
                if strcmp(alignment, 'cue')
                   cueTxt = 'cue';cueTxtOffset = 20;
                elseif strcmp(alignment, 'CC')
                    cueTxt = 'color change'; cueTxtOffset = 120;
                elseif strcmp (alignment, 'response')
                    cueTxt = 'response'; cueTxtOffset = 120;
                end 

                %make plot
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
                
            elseif plotFastSlowRT
                   %select sorted trials
                rtSort_pow_cat = pow_cat(:,:,sortedEvInd);
                %Collapse freqency
                freqLbl = [num2str(fRange(1)) ':' num2str(fRange(2))];
                [~,fInd_start] = min(abs(fRange(1) - powConfig.freQ));
                [~,fInd_end] = min(abs(fRange(2) - powConfig.freQ));
                fInd = fInd_start:fInd_end;
                rtSort_pow_cat = squeeze(nanmean(rtSort_pow_cat(fInd,:,:)))';
                
                if strcmp(alignment, 'cue')
                   cueTxt = 'cue';cueTxtOffset = 20;
                %make plot
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(sortedRTs),rtSort_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                swagAxes(gca,10,'time (ms)','trial',[freqName '_''ch ' channel ])
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                scatter(cell2mat({events(sortedEvInd).rt})+cell2mat({events(sortedEvInd).delay}),1:length(sortedRTs),'w','o');
                scatter(cell2mat({events(sortedEvInd).delay}),1:length(sortedRTs),'b','o');
                
                elseif strcmp(alignment, 'CC')
                    cueTxt = 'color change'; cueTxtOffset = 120;
                %make plot
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(sortedRTs),rtSort_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                swagAxes(gca,10,'time (ms)','trial',[freqName '_''ch ' channel ])
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                scatter(cell2mat({events(sortedEvInd).rt}),1:length(sortedRTs),'w','o');
                scatter(-1*cell2mat({events(sortedEvInd).delay}),1:length(sortedRTs),'b','o');

                elseif strcmp (alignment, 'response')
                    cueTxt = 'response'; cueTxtOffset = 120;
                end  
                if saveFigs
                    cd(figDir);
                    print(gcf,[channel '_' freqName '_FSRT'],'-dpng');close
                end
            else  
                fast_pow_cat = pow_cat(:,:,fastRT);
                slow_pow_cat = pow_cat(:,:,slowRT);
                fast_pow_cat = nanmean(fast_pow_cat,3);
                slow_pow_cat = nanmean(slow_pow_cat,3);
                freqLbl = [num2str(fRange(1)) ':' num2str(fRange(2))];
                [~,fInd_start] = min(abs(fRange(1) - powConfig.freQ));
                [~,fInd_end] = min(abs(fRange(2) - powConfig.freQ));
                fInd = fInd_start:fInd_end;
            
                freqSlow_pow_cat = pow_cat(fInd,:,slowRT); %get short delay events, only in fRange
                freqSlow_pow_cat = nanmean(freqSlow_pow_cat,3);% Average trails
                freqSlow_pow_cat = squeeze(nanmean(nanmean(freqSlow_pow_cat,3),1));% Average within fRange
                fRangeSlow_pow_cat = [fRangeSlow_pow_cat; freqSlow_pow_cat];
 
                freqFast_pow_cat = pow_cat(fInd,:,fastRT); %get long delay events, only in fRange
                freqFast_pow_cat = nanmean(freqFast_pow_cat,3);% Average trails
                freqFast_pow_cat = squeeze(nanmean(nanmean(freqFast_pow_cat,3),1));% Average within fRange
                fRangeFast_pow_cat = [fRangeFast_pow_cat; freqFast_pow_cat];
           
 
                if strcmp(alignment, 'cue')
                %make plots
                cueTxt = 'cue';cueTxtOffset = 20;
                cc1Txt = 'color change'; cc1TxtOffset = 75;
                cc2Txt = 'color change'; cc2TxtOffset = 50;
                
                %fast RT's
                fBins = powConfig.freqBins;
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(fBins),fast_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                setFreqYTicks(gca,fBins)
                swagAxes(gca,10,'time (ms)','frequency (Hz)',['ch ' channel])
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');
                line([500,500],ylim,'color','y');
                %line([1000,1000],ylim,'color','y');
                line([1500,1500],ylim,'color','y');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                text(1500-cc1TxtOffset,texty(2)+1,cc1Txt);
                text(500-cc2TxtOffset,texty(2)+1,cc2Txt);
                
                if saveFigs
                    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                    cd(saveFigsDir);
                    print(gcf,[channel '_fast'],'-dpng');close
                end
                
                %slow RTs
                fBins = powConfig.freqBins;
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(fBins),slow_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                setFreqYTicks(gca,fBins)
                swagAxes(gca,10,'time (ms)','frequency (Hz)',['ch ' channel])
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');
                line([500,500],ylim,'color','y');
                %line([1000,1000],ylim,'color','y');
                line([1500,1500],ylim,'color','y');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                text(1500-cc1TxtOffset,texty(2)+1,cc1Txt);
                text(500-cc2TxtOffset,texty(2)+1,cc2Txt);
                
                if saveFigs
                    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                    cd(saveFigsDir);
                    print(gcf,[channel '_slow'],'-dpng');close
                end
                
                elseif strcmp(alignment, 'response')
                fast_pow_cat = pow_cat(:,:,fastRT);
                slow_pow_cat = pow_cat(:,:,slowRT);
                fast_pow_cat = nanmean(fast_pow_cat,3);
                slow_pow_cat = nanmean(slow_pow_cat,3);
                
                %make plots
                cueTxt = 'response';cueTxtOffset = 20;
                
                %fast RT's
                fBins = powConfig.freqBins;
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(fBins),fast_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                setFreqYTicks(gca,fBins)
                swagAxes(gca,10,'time (ms)','frequency (Hz)')
                yl = ylim;
                title(['ch ' channel 'fast RT'],'Position',[-1000,yl(2)],'fontsize',20)                
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                
                if saveFigs
                    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                    cd(saveFigsDir);
                    print(gcf,[channel '_fast'],'-dpng');close
                end
                
                %slow RTs
                fBins = powConfig.freqBins;
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(fBins),slow_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                setFreqYTicks(gca,fBins)
                swagAxes(gca,10,'time (ms)','frequency (Hz)')
                yl = ylim;
                title(['ch ' channel 'slow RT'],'Position',[-1000,yl(2)],'fontsize',20)                
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                
                if saveFigs
                    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                    cd(saveFigsDir);
                    print(gcf,[channel '_slow'],'-dpng');close
                end
                end
            end
        case 'RT1000'
            plotFastSlowRT
                   %select sorted trials
                rtSort_pow_cat = pow_cat(:,:,sortedEvInd);
                %Collapse freqency
                freqLbl = [num2str(fRange(1)) ':' num2str(fRange(2))];
                [~,fInd_start] = min(abs(fRange(1) - powConfig.freQ));
                [~,fInd_end] = min(abs(fRange(2) - powConfig.freQ));
                fInd = fInd_start:fInd_end;
                rtSort_pow_cat = squeeze(nanmean(rtSort_pow_cat(fInd,:,:)))';
                
                if strcmp(alignment, 'cue')
                   cueTxt = 'cue';cueTxtOffset = 20;
                %make plot
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(sortedRTs),rtSort_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                swagAxes(gca,10,'time (ms)','trial',[freqName '_''ch ' channel ])
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                scatter(cell2mat({events(sortedEvInd).rt})+cell2mat({events(sortedEvInd).delay}),1:length(sortedRTs),'w','o');
                scatter(cell2mat({events(sortedEvInd).delay}),1:length(sortedRTs),'b','o');
                
                elseif strcmp(alignment, 'CC')
                    cueTxt = 'color change'; cueTxtOffset = 120;
                %make plot
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(sortedRTs),rtSort_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                swagAxes(gca,10,'time (ms)','trial',[freqName '_''ch ' channel ])
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                scatter(cell2mat({events(sortedEvInd).rt}),1:length(sortedRTs),'w','o');
                scatter(-1*cell2mat({events(sortedEvInd).delay}),1:length(sortedRTs),'b','o');

                elseif strcmp (alignment, 'response')
                    cueTxt = 'response'; cueTxtOffset = 120;
                end  
                if saveFigs
                    cd(figDir);
                    print(gcf,[channel '_' freqName '_FSRT'],'-dpng');close
                end
        case 'RTerr'
            if plotFastSlowRT
                   %select sorted trials
                rtSort_pow_cat = pow_cat(:,:,sortedEvInd);
                %Collapse freqency
                freqLbl = [num2str(fRange(1)) ':' num2str(fRange(2))];
                [~,fInd_start] = min(abs(fRange(1) - powConfig.freQ));
                [~,fInd_end] = min(abs(fRange(2) - powConfig.freQ));
                fInd = fInd_start:fInd_end;
                rtSort_pow_cat = squeeze(nanmean(rtSort_pow_cat(fInd,:,:)))';
                
                if strcmp(alignment, 'cue')
                   cueTxt = 'cue';cueTxtOffset = 20;
                elseif strcmp(alignment, 'CC')
                    cueTxt = 'color change'; cueTxtOffset = 120;
                elseif strcmp (alignment, 'response')
                    cueTxt = 'response'; cueTxtOffset = 120;
                end 

                %make plot
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(sortedRTs),rtSort_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                swagAxes(gca,10,'time (ms)','trial')
                yl = ylim;
                title([freqName '_''ch ' channel ],'Position',[-1000,yl(2)],'fontsize',20)            
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                scatter(cell2mat({events(sortedEvInd).rt})+cell2mat({events(sortedEvInd).delay}),1:length(sortedRTs),'w','o');
                scatter(cell2mat({events(sortedEvInd).delay}),1:length(sortedRTs),'b','o');
                if saveFigs
                    cd(figDir);
                    print(gcf,[channel '_' freqName 'FSRT'],'-dpng');close
                end
            end
        case 'RT1500'
             if plotFastSlowRT
                   %select sorted trials
                rtSort_pow_cat = pow_cat(:,:,sortedEvInd);
                %Collapse freqency
                freqLbl = [num2str(fRange(1)) ':' num2str(fRange(2))];
                [~,fInd_start] = min(abs(fRange(1) - powConfig.freQ));
                [~,fInd_end] = min(abs(fRange(2) - powConfig.freQ));
                fInd = fInd_start:fInd_end;
                rtSort_pow_cat = squeeze(nanmean(rtSort_pow_cat(fInd,:,:)))';
                
                if strcmp(alignment, 'cue')
                   cueTxt = 'cue';cueTxtOffset = 20;
                elseif strcmp(alignment, 'CC')
                    cueTxt = 'color change'; cueTxtOffset = 120;
                elseif strcmp (alignment, 'response')
                    cueTxt = 'response'; cueTxtOffset = 120;
                end 

                %make plot
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(sortedRTs),rtSort_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                swagAxes(gca,10,'time (ms)','trial',[freqName '_''ch ' channel])
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                if strcmp(alignment, 'cue')
                line([0,0],ylim,'color','w');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                scatter(cell2mat({events(sortedEvInd).rt})+cell2mat({events(sortedEvInd).delay}),1:length(sortedRTs),'w','o');
                scatter(cell2mat({events(sortedEvInd).delay}),1:length(sortedRTs),'b','o');
                  if saveFigs
                      cd(figDir);
                    print(gcf,[channel '_' freqName '_FSRT'],'-dpng');close
                  end
                
                elseif strcmp(alignment, 'CC')
                line([-1500,-1500],ylim,'color','w');
                texty = ylim;
                text(-1500, texty(2)+1, 'cue');
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                scatter(cell2mat({events(sortedEvInd).rt}),1:length(sortedRTs),'w','o');
                scatter(cell2mat({events(sortedEvInd).msoffset}),1:length(sortedRTs),'b','o');
                 
                
                if saveFigs
                    cd(figDir);
                    print(gcf,[channel '_' freqName '_FSRT'],'-dpng');close
                end
                end 
                
             else
                if strcmp(alignment, 'CC')
                fast_pow_cat = pow_cat(:,:,fastRT);
                slow_pow_cat = pow_cat(:,:,slowRT);
                fast_pow_cat = nanmean(fast_pow_cat,3);
                slow_pow_cat = nanmean(slow_pow_cat,3);
                
                %make plots
                cueTxt = 'color change';cueTxtOffset = 100;
                ccTxt = 'average response'; ccTxtOffset = 75;
              
                
                %fast RT's
                fBins = powConfig.freqBins;
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(fBins),fast_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                setFreqYTicks(gca,fBins)
                swagAxes(gca,10,'time (ms)','frequency (Hz)')
                yl = ylim;
                title(['ch ' channel 'fast RT'],'Position',[-1000,yl(2)],'fontsize',20)                
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
            line([525,525],ylim,'color','w');
            line ([375, 375], ylim, 'color', 'y');
            line ([725, 725], ylim, 'color', 'y');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                text(500-ccTxtOffset,texty(2)+1,ccTxt);
               
                
                if saveFigs
                    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                    cd(saveFigsDir);
                    print(gcf,[channel '_fast'],'-dpng');close
                end
                
                %slow RTs
                fBins = powConfig.freqBins;
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(fBins),slow_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                setFreqYTicks(gca,fBins)
                swagAxes(gca,10,'time (ms)','frequency (Hz)')
                yl = ylim;
                title(['ch ' channel 'slow RT' ],'Position',[-1000,yl(2)],'fontsize',20)                
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');

            line([525,525],ylim,'color','w');
            line ([375, 375], ylim, 'color', 'y');
            line ([725, 725], ylim, 'color', 'y');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                text(510-ccTxtOffset,texty(2)+1,ccTxt);
                
                if saveFigs
                    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                    cd(saveFigsDir);
                    print(gcf,[channel '_slow'],'-dpng');close
                end
                
            elseif strcmp (alignment, 'cue')
                fast_pow_cat = pow_cat(:,:,fastRT);
                slow_pow_cat = pow_cat(:,:,slowRT);
                fast_pow_cat = nanmean(fast_pow_cat,3);
                slow_pow_cat = nanmean(slow_pow_cat,3);
                
                %make plots
                cueTxt = 'cue';cueTxtOffset = 20;
                ccTxt = 'color change'; ccTxtOffset = 75;
              
                
                %fast RT's
                fBins = powConfig.freqBins;
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(fBins),fast_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                setFreqYTicks(gca,fBins)
                swagAxes(gca,10,'time (ms)','frequency (Hz)')
                yl = ylim;
                title(['ch ' channel 'fast RT'],'Position',[-1000,yl(2)],'fontsize',20)                
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');
                line([1500,1500],ylim,'color','y');
                
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                text(1500-ccTxtOffset,texty(2)+1,ccTxt);
               
                
                if saveFigs
                    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                    cd(saveFigsDir);
                    print(gcf,[channel '_fast'],'-dpng');close
                end
                
                %slow RTs
                fBins = powConfig.freqBins;
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(fBins),slow_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                setFreqYTicks(gca,fBins)
                swagAxes(gca,10,'time (ms)','frequency (Hz)')
                yl = ylim;
                title(['ch ' channel 'slow RT' ],'Position',[-1000,yl(2)],'fontsize',20)                
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');
                line([1500, 1500],ylim,'color','y');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                text(1500-ccTxtOffset,texty(2)+1,ccTxt);
                
                if saveFigs
                    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                    cd(saveFigsDir);
                    print(gcf,[channel '_slow'],'-dpng');close
                end
                
                elseif strcmp(alignment, 'response')
                fast_pow_cat = pow_cat(:,:,fastRT);
                slow_pow_cat = pow_cat(:,:,slowRT);
                fast_pow_cat = nanmean(fast_pow_cat,3);
                slow_pow_cat = nanmean(slow_pow_cat,3);
                
                %make plots
                cueTxt = 'response';cueTxtOffset = 20;

                %fast RT's
                fBins = powConfig.freqBins;
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(fBins),fast_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                setFreqYTicks(gca,fBins)
                swagAxes(gca,10,'time (ms)','frequency (Hz)')
                yl = ylim;
                title(['ch ' channel 'fast RT'],'Position',[-1000,yl(2)],'fontsize',20)                
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
               
                
                if saveFigs
                    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                    cd(saveFigsDir);
                    print(gcf,[channel '_fast'],'-dpng');close
                end
                
                %slow RTs
                fBins = powConfig.freqBins;
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(fBins),slow_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                setFreqYTicks(gca,fBins)
                swagAxes(gca,10,'time (ms)','frequency (Hz)')
                yl = ylim;
                title(['ch ' channel 'slow RT' ],'Position',[-1000,yl(2)],'fontsize',20)                
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                
                if saveFigs
                    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                    cd(saveFigsDir);
                    print(gcf,[channel '_slow'],'-dpng');close
                end
                end
            
            end
            
           
        case 'RT500'
            if plotFastSlowRT
                   %select sorted trials
                rtSort_pow_cat = pow_cat(:,:,sortedEvInd);
                %Collapse freqency
                freqLbl = [num2str(fRange(1)) ':' num2str(fRange(2))];
                [~,fInd_start] = min(abs(fRange(1) - powConfig.freQ));
                [~,fInd_end] = min(abs(fRange(2) - powConfig.freQ));
                fInd = fInd_start:fInd_end;
                rtSort_pow_cat = squeeze(nanmean(rtSort_pow_cat(fInd,:,:)))';
                
                if strcmp(alignment, 'cue')
                   cueTxt = 'cue';cueTxtOffset = 20;
                elseif strcmp(alignment, 'CC')
                    cueTxt = 'color change'; cueTxtOffset = 120;
                elseif strcmp (alignment, 'response')
                    cueTxt = 'response'; cueTxtOffset = 120;
                end 

                %make plot
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(sortedRTs),rtSort_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                swagAxes(gca,10,'time (ms)','trial',[freqName '-' 'ch ' channel ])
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                
                if strcmp(alignment, 'cue')
                line([0,0],ylim,'color','w');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                scatter(cell2mat({events(sortedEvInd).rt})+cell2mat({events(sortedEvInd).delay}),1:length(sortedRTs),'w','o');
                scatter(cell2mat({events(sortedEvInd).delay}),1:length(sortedRTs),'b','o');
                
                if saveFigs
                    cd(figDir);
                    print(gcf,[channel '_' freqName '_FSRT'],'-dpng');close
                end
                
                elseif strcmp(alignment, 'CC')
                line([-500,-500],ylim,'color','w');
                texty = ylim;
                text(-500, texty(2)+1, 'cue');
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                scatter(cell2mat({events(sortedEvInd).rt}),1:length(sortedRTs),'w','o');
                scatter(cell2mat({events(sortedEvInd).msoffset}),1:length(sortedRTs),'b','o');
                if saveFigs
                    cd(figDir);
                    print(gcf,[channel '_' freqName '_FSRT'],'-dpng');close
                end
                end
            
            else
             if strcmp(alignment, 'CC')
                fast_pow_cat = pow_cat(:,:,fastRT);
                slow_pow_cat = pow_cat(:,:,slowRT);
                fast_pow_cat = nanmean(fast_pow_cat,3);
                slow_pow_cat = nanmean(slow_pow_cat,3);
                
                %make plots
                cueTxt = 'color change';cueTxtOffset = 100;
                ccTxt = 'average response'; ccTxtOffset = 75;
              
                
                %fast RT's
                fBins = powConfig.freqBins;
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(fBins),fast_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                setFreqYTicks(gca,fBins)
                swagAxes(gca,10,'time (ms)','frequency (Hz)')
                yl = ylim;
                title(['ch ' channel 'fast RT'],'Position',[-1000,yl(2)],'fontsize',20)                
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');
            line([525,525],ylim,'color','w');
            line ([375, 375], ylim, 'color', 'y');
            line ([725, 725], ylim, 'color', 'y');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                text(525-ccTxtOffset,texty(2)+1,ccTxt);
                              
                
                if saveFigs
                    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                    cd(saveFigsDir);
                    print(gcf,[channel '_fast'],'-dpng');close
                end
                
                %slow RTs
                fBins = powConfig.freqBins;
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(fBins),slow_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                setFreqYTicks(gca,fBins)
                swagAxes(gca,10,'time (ms)','frequency (Hz)')
                yl = ylim;
                title(['ch ' channel 'slow RT'],'Position',[-1000,yl(2)],'fontsize',20)                
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');
            line([525,525],ylim,'color','w');
            line ([375, 375], ylim, 'color', 'y');
            line ([725, 725], ylim, 'color', 'y');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                text(525-ccTxtOffset,texty(2)+1,ccTxt);
                
                if saveFigs
                    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                    cd(saveFigsDir);
                    print(gcf,[channel '_slow'],'-dpng');close
                end
            elseif strcmp (alignment, 'cue')
                fast_pow_cat = pow_cat(:,:,fastRT);
                slow_pow_cat = pow_cat(:,:,slowRT);
                fast_pow_cat = nanmean(fast_pow_cat,3);
                slow_pow_cat = nanmean(slow_pow_cat,3);
                
                %make plots
                cueTxt = 'cue';cueTxtOffset = 20;
                ccTxt = 'color change'; ccTxtOffset = 50;
              
                
                %fast RT's
                fBins = powConfig.freqBins;
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(fBins),fast_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                setFreqYTicks(gca,fBins)
                swagAxes(gca,10,'time (ms)','frequency (Hz)')
                yl = ylim;
                title(['ch ' channel 'fast RT'],'Position',[-1000,yl(2)],'fontsize',20)                
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
                
                %slow RTs
                fBins = powConfig.freqBins;
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(fBins),slow_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                setFreqYTicks(gca,fBins)
                swagAxes(gca,10,'time (ms)','frequency (Hz)')
                yl = ylim;
                title(['ch ' channel 'slow RT' ],'Position',[-1000,yl(2)],'fontsize',20)                
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');
                line([500, 500],ylim,'color','y');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                text(500-ccTxtOffset,texty(2)+1,ccTxt);
                
                if saveFigs
                    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                    cd(saveFigsDir);
                    print(gcf,[channel '_slow'],'-dpng');close
                end
                
                elseif strcmp(alignment, 'response')
                fast_pow_cat = pow_cat(:,:,fastRT);
                slow_pow_cat = pow_cat(:,:,slowRT);
                fast_pow_cat = nanmean(fast_pow_cat,3);
                slow_pow_cat = nanmean(slow_pow_cat,3);
                
                %make plots
                cueTxt = 'response';cueTxtOffset = 20;
             
                %fast RT's
                fBins = powConfig.freqBins;
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(fBins),fast_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                setFreqYTicks(gca,fBins)
                swagAxes(gca,10,'time (ms)','frequency (Hz)')
                yl = ylim;
                title(['ch ' channel 'fast RT'],'Position',[-1000,yl(2)],'fontsize',20)                
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                             
                if saveFigs
                    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                    cd(saveFigsDir);
                    print(gcf,[channel '_fast'],'-dpng');close
                end
                
                %slow RTs
                fBins = powConfig.freqBins;
                tBins = nanmean(powConfig.timeBins,2)';
                swagFig([1/8 1/4 3/4 1/2]); hold all
                imagesc(tBins,1:length(fBins),slow_pow_cat);
                set(gca,'ydir','normal','CLim',[-3 3]);
                setFreqYTicks(gca,fBins)
                swagAxes(gca,10,'time (ms)','frequency (Hz)')
                yl = ylim;
                title(['ch ' channel 'slow RT' ],'Position',[-1000,yl(2)],'fontsize',20)                
                hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
                line([0,0],ylim,'color','w');
                texty = ylim;
                text(0-cueTxtOffset,texty(2)+1,cueTxt);
                
                if saveFigs
                    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                    cd(saveFigsDir);
                    print(gcf,[channel '_slow'],'-dpng');close
                end
            
             end 
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
            
            elseif strcmp(alignment, 'CC')
                 cueTxt = 'color change'; cueTxtOffset = 120;
                 ccTxt = 'response average'; ccTxtOffset = 50;
                  
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
            line([525,525],ylim,'color','w');
            line ([375, 375], ylim, 'color', 'y');
            line ([725, 725], ylim, 'color', 'y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(525-ccTxtOffset,texty(2)+1,ccTxt);
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
            line([525,525],ylim,'color','w');
            line ([375, 375], ylim, 'color', 'y');
            line ([725, 725], ylim, 'color', 'y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(525-ccTxtOffset,texty(2)+1,ccTxt);
              if saveFigs
                 if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                 cd(saveFigsDir);
                 print(gcf,[channel '_late'],'-dpng');close
              end


       elseif strcmp(alignment, 'response')
            cueTxt = 'response '; cueTxtOffset = 120;

                  
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
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);

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
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);

              if saveFigs
                 if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                 cd(saveFigsDir);
                 print(gcf,[channel '_late'],'-dpng');close
              end
             
            end 
            
        case 'timeLong'
            early_pow_cat = [];
            late_pow_cat = [];
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
            %line([1000,1000],ylim,'color','w');
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
            %line([1000,1000],ylim,'color','w');
            line([1500,1500],ylim,'color','y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(1500-ccTxtOffset,texty(2)+1,ccTxt);
            
            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_late'],'-dpng');close
            end
            
            elseif strcmp(alignment, 'CC')
                 cueTxt = 'color change'; 
                 cueTxtOffset = 120;
                 ccTxt = 'average response'; 
                 ccTxtOffset = 50;
                  
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
            line([525,525],ylim,'color','w');
            line ([375, 375], ylim, 'color', 'y');
            line ([725, 725], ylim, 'color', 'y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(525-ccTxtOffset,texty(2)+1,ccTxt);
            
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
            line([525,525],ylim,'color','w');
            line ([375, 375], ylim, 'color', 'y');
            line ([725, 725], ylim, 'color', 'y');
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            text(525-ccTxtOffset,texty(2)+1,ccTxt);
            
             if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_late'],'-dpng');close
             end
             
        elseif strcmp(alignment, 'response')
            cueTxt = 'response'; cueTxtOffset = 120;
                  
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
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            
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
            texty = ylim;
            text(0-cueTxtOffset,texty(2)+1,cueTxt);
            
             if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,[channel '_late'],'-dpng');close
             end   
            end
    end
        
end
% Make plots
if strcmp(alignment, 'cue')
    cueTxt = 'cue';cueTxtOffset = 20;
    ccTxt = 'color change';ccTxtOffset = 75;
elseif strcmp(alignment, 'CC')
    cueTxt = 'color change'; cueTxtOffset = 120;
    ccTxt = 'average response'; ccTxtOffset = 90;
elseif strcmp(alignment, 'response')
    cueTxt = 'response'; cueTxtOffset = 120;
end 
    
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
    case 'timeShort'
        freq_pow_cat1 = fRangeEarly_pow_cat;
        freq_pow_cat2 = fRangeLate_pow_cat;
        tit1 = [freqName 'Power Early Trials'];
        tit2 = [freqName 'Power Late Trials'];
        fileName1 = [freqName '_early'];
        fileName2 = [freqName '_late'];
    case 'timeLong'
        freq_pow_cat1 = fRangeEarly_pow_cat;
        freq_pow_cat2 = fRangeLate_pow_cat;
        tit1 = [freqName 'Power Early Trials'];
        tit2 = [freqName 'Power Late Trials'];
        fileName1 = [freqName '_early'];
        fileName2 = [freqName '_late'];
     case 'RT' 
         freq_pow_cat1 = fRangeSlow_pow_cat;
         freq_pow_cat2 = fRangeFast_pow_cat;
         tit1 = [freqName 'Power Slow Response'];
         tit2 = [freqName 'Power Fast Response'];
         fileName1 = [freqName '-slowRT'];
         fileName2 = [freqName '-fastRT'];
      case 'RT500' 
         freq_pow_cat1 = fRangeSlow_pow_cat;
         freq_pow_cat2 = fRangeFast_pow_cat;
         tit1 = [freqName 'Power Slow Response'];
         tit2 = [freqName 'Power Fast Response'];
         fileName1 = [freqName '-slowRT'];
         fileName2 = [freqName '-fastRT'];  
        case 'RT1500' 
         freq_pow_cat1 = fRangeSlow_pow_cat;
         freq_pow_cat2 = fRangeFast_pow_cat;
         tit1 = [freqName 'Power Slow Response'];
         tit2 = [freqName 'Power Fast Response'];
         fileName1 = [freqName '-slowRT'];
         fileName2 = [freqName '-fastRT'];
end

if rawPlot
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
%line([1000,1000],ylim,'color','w');
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
%line([500,500],ylim,'color','w');
%line([1000,1000],ylim,'color','w');
%line([1500,1500],ylim,'color','y');
texty = ylim;
text(0-cueTxtOffset,texty(2)+1,cueTxt);
%text(1500-ccTxtOffset,texty(2)+1,ccTxt);
%text(500-ccTxtOffset,texty(2)+1,ccTxt);% short delay

if saveFigs
    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
    cd(saveFigsDir);
    print(gcf,fileName2,'-dpng');close
end

% Subtract plots
diff_pow_cat = freq_pow_cat2 - freq_pow_cat1;
tit = [freqName 'Power Late - Early Trials'];
tBins = nanmean(powConfig.timeBins,2)';
swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
imagesc(tBins,1:length(channels),diff_pow_cat);
set(gca,'ydir','normal','CLim',[-3 3]);
yt = [1:1:length(channels)];
set(gca,'ytick',yt,'yticklabel',channels,'yticklabelrotation',0)
swagAxes(gca,10,'time (ms)','electrode',tit)
hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
line([0,0],ylim,'color','w');
%line([500,500],ylim,'color','y');% short delay
%line([500,500],ylim,'color','w');
%line([1000,1000],ylim,'color','w');
%line([1500,1500],ylim,'color','y');
texty = ylim;
text(0-cueTxtOffset,texty(2)+1,cueTxt);
%text(1500-ccTxtOffset,texty(2)+1,ccTxt);
%text(500-ccTxtOffset,texty(2)+1,ccTxt);% short delay

if saveFigs
    if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
    cd(saveFigsDir);
    print(gcf,[freqName '_' comparison '_diff'],'-dpng');close
end

tBins = nanmean(powConfig.timeBins,2)';
swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
imagesc(tBins,1:length(channels),diff_theta_pow_cat);
set(gca,'ydir','normal','CLim',[-3 3]);
yt = [1:1:length(channels)];
set(gca,'ytick',yt,'yticklabel',channels,'yticklabelrotation',0)
swagAxes(gca,10,'time (ms)','electrode','Theta Power Short-Long Delay')
hc = colorbar; set(get(hc,'YLabel'),'String','power'); colormap(jet);
line([0,0],ylim,'color','w');
line([500,500],ylim,'color','w');
line([1000,1000],ylim,'color','w');
line([1500,1500],ylim,'color','y');
texty = ylim;
text(0-cueTxtOffset,texty(2)+1,cueTxt);
text(1500-ccTxtOffset,texty(2)+1,ccTxt);
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
end 