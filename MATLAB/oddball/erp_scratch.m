 %% Inputs 
clear;
% subjList = {'HUP142_i','HUP144_e','HUP145_e','HUP147_e','HUP148_e','HUP149_e','HUP150_i','HUP151_e',...
%     'HUP152_e','HUP153_i','HUP154_e','HUP155_i','HUP156_i','HUP157_e','HUP159_e','HUP165_i','HUP166_i'};
subjList = {'HUP149_e'};%,'HUP157_e','HUP159_e','HUP165_i','HUP166_i'};
for s=1:length(subjList)
    subj = subjList{s};
    %ddir = fullfile('D:\TNL_Data\oddball\eeg',subj,'processed');
    ddir = fullfile('/Volumes/HumanStudies/HumanStudies/oddball/eeg',subj,'processed');
    elec_coordinates_file = fullfile('/Volumes/HumanStudies/HumanStudies/localization',subj(1:end-2),'post/electrodenames_coordinates_native.csv');
    trialType1 = {'TARGETHF','TARGETLF','BACKGROUNDHF','BACKGROUNDLF'};
    trialType2 = {}; % leave empty to only select trialType1
    subtract_trialTypes = 0; %set to 1 to plot difference btw types
    leadAvg_reref = 1; %re-reference to lead average
    zscore = 1; %z-score electrodes to pre-trial baseline
    
    % Use these to check different points in data manipulation
    plotLeadSingleTrials = 0; % make single trial plots for each lead
    overlayReref = 1; % plot rereference result on top of raw signal
    plotLeadAvgTrials = 0; % plot the average across trials
    elecTimePlot = 1; % imagesc plot of single electrode over time
    zScoreDemoPlots = 0; %single channel z score subplot demonstration
    respDistPlots = 1; %plot distance from first responding electrode x response time
    savePlots = 1;
    saveDir = fullfile('/Volumes/HumanStudies/HumanStudies/oddball/scratch/erp',subj);
    if ~exist(saveDir,'dir'),mkdir(saveDir);end
    
    % Plot title/legend labels
    plot_title = [];
    trialType1_label = 'All Tones';
    trialType2_label = [];
    load([ddir '/sessInfo.mat']);
    if exist(fullfile(ddir,'manual_bad_channels.mat'),'file')==2
        load(fullfile(ddir,'manual_bad_channels.mat'),'bad_channels')
    end
    includeTrials = [1 size(trial_type,2)];
    if elecTimePlot
        load('/Volumes/HumanStudies/HumanStudies/oddball/scratch/dci/plotEventsStruct.mat');
        if ismember(subj,{plotEventsStruct.subj})
            [~, sInd] = ismember(subj,{plotEventsStruct.subj});
            nlxEventsToPlot = plotEventsStruct(sInd).events;
        else
            disp('Select neuralynx events to plot for this patient')
            keyboard; % Edit the events field and run the next 3 lines:
            plotEventsStruct(length(plotEventsStruct)+1).subj = subj;
            plotEventsStruct(length(plotEventsStruct)).events = [8,9,10,14,15,17,19,22];
            save('/Volumes/HumanStudies/HumanStudies/oddball/scratch/dci/plotEventsStruct.mat','plotEventsStruct');
        end
    end
    %keyboard % manually set and run include trials based on nlxEvents
    %manually edit bad channels and resave file
    %bad_channels = [bad_channels,25];
    %save(fullfile(ddir,'manual_bad_channels.mat'),'bad_channels');
    
    %% Setup
    % select trials
    includeTrials1 = ismember(trial_type,trialType1);
    if ~isempty(trialType2)
        includeTrials2 = ismember(trial_type,trialType2);
    end
    if includeTrials(1)~=1 %trim early trials
        includeTrials1(1:includeTrials(1)) = false;
        if ~isempty(trialType2)
            includeTrials2(1:includeTrials(1)) = false;
        end
    end
    if includeTrials(2)~=length(trial_type) %trim late trials
        includeTrials1(includeTrials(2):end) = false;
        if ~isempty(trialType2)
            includeTrials2(includeTrials(2):end) = false;
        end
    end
    
    % initialize erp vectors (chan x t)
    erp1 = ones(size(elecInfo,1)-length(bad_channels),size(t,2))*NaN;
    if leadAvg_reref
        rerp1 = ones(size(elecInfo,1)-length(bad_channels),size(t,2))*NaN;
    end
    if ~isempty(trialType2)
        erp2 = ones(size(elecInfo,1)-length(bad_channels),size(t,2))*NaN;
        if leadAvg_reref
            rerp2 = ones(size(elecInfo,1)-length(bad_channels),size(t,2))*NaN;
        end
    end
    
    %% load each lead, rereference, average across trials and store in erp
    goodChan = 0;
    erp_elecInfo = {};
    files  = dir(ddir);
    for f=1:length(files)
        if ~ismember(files(f).name,{'.','..','sessInfo.mat','.DS_Store'})
            load(fullfile(ddir,files(f).name));
            includeChan = ~ismember(cell2mat(lead_elecInfo(:,1)),bad_channels);
            if plotLeadSingleTrials
                skipLead = 0;
                for trial=includeTrials(1):includeTrials(2)
                    if skipLead==1,break;end
                    cShift = 1.1*max(max(abs(eeg(1,:,:))));
                    cpos = [];
                    bc=0;
                    figure('Position',[0 0 1920 1080])
                    hold on
                    lead_avg = nanmean(eeg(trial,:,includeChan),3);
                    reeg = bsxfun(@minus,eeg(trial,:,:),lead_avg);
                    for c=1:length(includeChan)
                        if ~includeChan(c),bc=bc+1;continue;end
                        plot(t,eeg(trial,:,c)+((c-1-bc)*cShift),'Color',[0, 0.4470, 0.7410]);
                        if overlayReref
                            plot(t,reeg(1,:,c)+((c-1-bc)*cShift),'Color',[0.8500, 0.3250, 0.0980]);
                        end
                        plot([t(1), t(end)],[((c-1-bc)*cShift), ((c-1-bc)*cShift)],'k','LineWidth',0.5)
                        cpos= [cpos, ((c-1)*cShift)];
                    end
                    legend('raw','reref')
                    plot([t(1), t(end)],[((c-bc)*cShift), ((c-bc)*cShift)],'k','LineWidth',0.5)
                    plot(t,lead_avg+((c-bc)*cShift),'Color',[0.4660, 0.6740, 0.1880]);
                    plot([0, 0],ylim,'k')
                    cpos = [cpos,((c-bc)*cShift)];
                    clbls = [lead_elecInfo(includeChan,3); 'Avg'];
                    set(gca,'YTick',cpos,'YTickLabels',clbls,'FontSize',20);
                    title([subj ' Trial: ' num2str(trial) ' (' trial_type{trial} ')'],'Interpreter', 'none');
                    xlabel('Time (ms)');
                    ylabel('Channel');
                    keyboard
                end
            end
            if leadAvg_reref && elecTimePlot
                lead_avg = nanmean(eeg(:,:,includeChan),3);
                reeg=bsxfun(@minus,eeg,lead_avg);
                rerp1(goodChan+1:goodChan+sum(includeChan),:) = squeeze(nanmean(reeg(includeTrials1,:,includeChan),1))';
                winPerSess = 50;
                chanPlotDir = [saveDir,'/channelTimePlots/' num2str(winPerSess) '/'];
                if ~exist(chanPlotDir,'dir')
                    mkdir(chanPlotDir)
                end
                trialsPerWin = floor(size(eeg,1)/winPerSess);
                excludeTrials = size(eeg,1)-(winPerSess*trialsPerWin); %used to exclude early trials rather than late
                [~,toneInd] = min(abs(t-0));
                for c=1:sum(includeChan)
                    winErp = ones(floor(sum(includeTrials1)/trialsPerWin),size(t,2))*NaN;
                    trialInd = max(find(includeTrials1,excludeTrials))+1; 
                    for trialWin = 1:floor(sum(includeTrials1)/trialsPerWin)
                        winTrials = [];
                        while length(winTrials)<trialsPerWin
                            if includeTrials1(trialInd)
                                winTrials = [winTrials,trialInd];
                            end
                            trialInd = trialInd+1;
                        end
                        winErp(trialWin,:) = squeeze(nanmean(reeg(winTrials,:,c),1))';
                    end
                    winErp(~any(~isnan(winErp), 2),:)=[];
                    winErpMean = mean(winErp(:,t<0),2);
                    winErpStd = std(winErp(:,t<0),0,2);
                    winErp = bsxfun(@minus,winErp,winErpMean);
                    winErp = bsxfun(@rdivide,winErp,winErpStd);
                    figure('Position',[0 0 960 1080]);
                    imagesc(t,[1:size(winErp,1)],winErp);
                    hold on
                    nlxEvHandles = [];
                    for ii=1:length(nlxEventsToPlot)
                        blkNum = ceil((nlxEvents(nlxEventsToPlot(ii)).approxTrial-excludeTrials)/trialsPerWin);
                        h=plot([t(1) 0], [blkNum blkNum]);
                        nlxEvHandles = [nlxEvHandles,h];
                    end
                    xlabel('Time (ms)');
                    ylabel('Trial Block');
                    cb = colorbar;
                    ylabel(cb, 'Voltage');
                    title([subj '   Channel ' lead_elecInfo{c,3} '   Trials per Row: ' num2str(trialsPerWin)],'Interpreter','none')
                    caxis([-8 8]);
                    nlxEvLabels = {nlxEvents(nlxEventsToPlot).type};
                    legend(nlxEvHandles,nlxEvLabels,'Location','northeast')
                    if savePlots
                        print(gcf,[chanPlotDir lead_elecInfo{c,3} '.png'],'-dpng');
                        close
                    end
                end
            end
            erp1(goodChan+1:goodChan+sum(includeChan),:) = squeeze(nanmean(eeg(includeTrials1,:,includeChan),1))';
            if ~isempty(trialType2)
                erp2(goodChan+1:goodChan+sum(includeChan),:) = squeeze(nanmean(eeg(includeTrials2,:,includeChan),1))';
                if leadAvg_reref
                    rerp2(goodChan+1:goodChan+sum(includeChan),:) = squeeze(nanmean(reeg(includeTrials2,:,includeChan),1))';
                end
            end
            if plotLeadAvgTrials
                cShift = 1.1*max(max(abs(erp1(:,:))));
                cpos = [];
                figure('Position',[0 0 1920 1080])
                hold on
                for c=goodChan+1:goodChan+sum(includeChan)
                    plot(t,erp1(c,:)+((c-1)*cShift),'Color',[0, 0.4470, 0.7410]);
                    if overlayReref
                        plot(t,rerp1(c,:)+((c-1)*cShift),'Color',[0.8500, 0.3250, 0.0980]);
                    end
                    plot([t(1), t(end)],[((c-1)*cShift), ((c-1)*cShift)],'k','LineWidth',0.5)
                    cpos= [cpos, ((c-1)*cShift)];
                end
                legend('raw','reref')
                plot([0, 0],ylim,'k')
                clbls = lead_elecInfo(includeChan,3);
                set(gca,'YTick',cpos,'YTickLabels',clbls,'FontSize',20);
                title([subj ' Average ' trialType1_label],'Interpreter', 'none');
                xlabel('Time (ms)');
                ylabel('Channel');
                keyboard
            end
            goodChan = goodChan+sum(includeChan);
            erp_elecInfo = [erp_elecInfo;lead_elecInfo(includeChan,:)];
        end
    end
    if subtract_trialTypes
        mmn=erp1-erp2;
    end
    if zscore
        if subtract_trialTypes
            mmnmean = mean(mmn(:,t<0),2);
            mmnstd = std(erp1(:,t<0),0,2);
            zmmn = (mmn-mmnmean)./mmnstd;
        end
        erp1mean = mean(erp1(:,t<0),2);
        erp1std = std(erp1(:,t<0),0,2);
        rerp1mean = mean(rerp1(:,t<0),2);
        rerp1std = std(rerp1(:,t<0),0,2);
        zerp1 = bsxfun(@minus,erp1,erp1mean);
        zerp1 = bsxfun(@rdivide,zerp1,erp1std);
        zrerp1 = bsxfun(@minus,rerp1,rerp1mean);
        zrerp1 = bsxfun(@rdivide,zrerp1,rerp1std);
        if zScoreDemoPlots
            for sChan=1:size(erp1,1)
                clbl = erp_elecInfo{sChan,3};
                figure('Position',[0 0 1920 1080])
                ax1=subplot(3,2,1);
                plot(t,erp1(sChan,:),'Color',[0, 0.4470, 0.7410]);
                hold on
                plot([t(1) t(sum(t<0))],[erp1mean(sChan) erp1mean(sChan)],'Color',[0.8500, 0.3250, 0.0980]);
                plot([t(1) t(sum(t<0))],[erp1mean(sChan)-erp1std(sChan) erp1mean(sChan)-erp1std(sChan)],'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--');
                plot([t(1) t(sum(t<0))],[erp1mean(sChan)+erp1std(sChan) erp1mean(sChan)+erp1std(sChan)],'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--');
                title(['Channel ' clbl ' Average ' trialType1_label ' With Pre-Cue Dist'])
                ylabel('Voltage')
                xlabel('Time');
                ax2=subplot(3,2,3);
                plot(t,zerp1(sChan,:))
                hold on
                plot([t(1) t(end)],[0 0],'k')
                plot([t(1) t(end)],[2 2],'k--')
                plot([t(1) t(end)],[-2 -2],'k--')
                xlabel('Time');
                ylabel('Z-Score');
                title('Normalized to Pre-Cue, Emphasizing 2 Standard deviations')
                ax3=subplot(3,2,5);
                imagesc(t,1,zerp1(sChan,:))
                xlabel('Time');
                set(gca,'YTick',1,'YTickLabel',clbl);
                title('Scaled Z-Score Colormap: caxis = [-6 6]');
                linkaxes([ax1,ax2,ax3],'x')
                
                ax4=subplot(3,2,2);
                plot(t,rerp1(sChan,:),'Color',[0, 0.4470, 0.7410]);
                hold on
                plot([t(1) t(sum(t<0))],[rerp1mean(sChan) rerp1mean(sChan)],'Color',[0.8500, 0.3250, 0.0980]);
                plot([t(1) t(sum(t<0))],[rerp1mean(sChan)-rerp1std(sChan) rerp1mean(sChan)-rerp1std(sChan)],'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--');
                plot([t(1) t(sum(t<0))],[rerp1mean(sChan)+rerp1std(sChan) rerp1mean(sChan)+rerp1std(sChan)],'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--');
                title(['Channel ' clbl ' Reref Average ' trialType1_label ' With Pre-Cue Dist'])
                ylabel('Voltage')
                xlabel('Time');
                ax5=subplot(3,2,4);
                plot(t,zrerp1(sChan,:))
                hold on
                plot([t(1) t(end)],[0 0],'k')
                plot([t(1) t(end)],[2 2],'k--')
                plot([t(1) t(end)],[-2 -2],'k--')
                title('Normalized to Pre-Cue, Emphasizing 2 Standard deviations')
                xlabel('Time');
                ylabel('Z-Score');
                ax6=subplot(3,2,6);
                imagesc(t,1,zrerp1(sChan,:))
                caxis([-6 6]);
                xlabel('Time');
                set(gca,'YTick',1,'YTickLabel',clbl);
                title('Scaled Z-Score Colormap: caxis = [-6 6]');
                linkaxes([ax4,ax5,ax6],'x')
                keyboard
            end
        end
        if ~isempty(trialType2)
            erp2mean = mean(erp2(:,t<0),2);
            erp2std = std(erp2(:,t<0),0,2);
            zerp2 = (erp2-erp2mean)./erp2std;
        end
    end
    %% Make Plots
    if subtract_trialTypes
        if zscore==0
            h=imagesc(mmn);
        else
            h=imagesc(zmmn);
        end
    else
        if zscore==0
            h=imagesc(erp1);
        else
            % no rereference
            figure('Position',[0 0 960 1080]);
            h=imagesc(t,1:size(erp_elecInfo,1),zerp1);
            title([subj ' Average Trial Voltage ' trialType1_label],'Interpreter','none');
            c = colorbar;
            if zscore == 0
                ylabel(c, 'Voltage')
            else
                ylabel(c, 'Normalized Voltage')
            end
            xlabel('Time');
            set(gca,'YTick',1:size(erp_elecInfo,1),'YTickLabel', erp_elecInfo(:,3))
            caxis([-8 8])
            if savePlots
                print(gcf,[saveDir,'/ERP.png'],'-dpng');
                close;
            end
            
            %use rereference
            figure('Position',[0 0 960 1080]);
            h=imagesc(t,1:size(erp_elecInfo,1),zrerp1);
            title([subj ' Average, Rereferenced Trial Voltage ' trialType1_label],'Interpreter','none');
            c = colorbar;
            if zscore == 0
                ylabel(c, 'Voltage')
            else
                ylabel(c, 'Normalized Voltage')
            end
            xlabel('Time');
            set(gca,'YTick',1:size(erp_elecInfo,1),'YTickLabel', erp_elecInfo(:,3))
            caxis([-8 8])
            if savePlots
                print(gcf,[saveDir,'/zERP.png'],'-dpng');
            end
        end
    end
    %% Sort by first change post-tone
    keyboard
    [~,toneInd] = min(abs(t-0));
    zThresh = 4; % find z score values above/below this threshold
    consecutiveThresh = 10; % number of significant samples in a row
    respInd = zeros(size(zrerp1,1),1);
    for chan=1:size(zrerp1,1)
        startSig = strfind(abs(zrerp1(chan,toneInd:end))>zThresh,ones(1,consecutiveThresh));
        if ~isempty(startSig)
            respInd(chan) = toneInd+startSig(1)-1;
        else
            respInd(chan) = size(zrerp1,2);
        end
    end
    [respInd,channelOrder] = sort(respInd);
    respSeq_elecInfo = erp_elecInfo(channelOrder,:);
    hold on
    plot([t(respInd(1)) t(respInd(1))], ylim)
    if savePlots,close;end
    
    %% reorder z-score erp plot by response time
    figure('Position',[0 0 960 1080]);
    h=imagesc(t,1:size(erp_elecInfo,1),zrerp1(channelOrder,:));
    title([subj ' Average, Rereferenced Trial Voltage ' trialType1_label],'Interpreter','none');
    c = colorbar;
    if zscore == 0
        ylabel(c, 'Voltage')
    else
        ylabel(c, 'Normalized Voltage')
    end
    xlabel('Time');
    set(gca,'YTick',1:size(respSeq_elecInfo,1),'YTickLabel', respSeq_elecInfo(:,3))
    caxis([-8 8])
    if savePlots
        print(gcf,[saveDir,'/zERP_respTime.png'],'-dpng');
        close;
    end
    
%% get centroids
    fid = fopen(elec_coordinates_file);
    if fid==-1
        disp(['Missing coordinates for ' subj]);
        continue
    end
    thisLine = fgetl(fid);
    centroidLbls = {};
    centroids = [];
    while ischar(thisLine)
        elecLine = strsplit(thisLine,',');
        centroidLbls = [centroidLbls;elecLine{1}];
        if isnan(str2double(elecLine{2}))
            centroids = [centroids; str2double(elecLine{3}),str2double(elecLine{4}),str2double(elecLine{5})];
        else
            centroids = [centroids; str2double(elecLine{2}),str2double(elecLine{3}),str2double(elecLine{4})];
        end
        thisLine = fgetl(fid);
    end
    fclose(fid);
    
    numSeeds = 3; %number of origins in each hemisphere to plot
    for hemi=1:2
        if hemi ==1
            hlbl = 'L';
        else
            hlbl = 'R';
        end
        for chan=1:size(respSeq_elecInfo,1)
            if strcmp(respSeq_elecInfo{chan,3}(1),hlbl) && respInd(chan)<size(zrerp1,2);
                hemiInd(chan)=true; %select only electodes in this hemisphere that respond
            else
                hemiInd(chan)=false;
            end
        end
        if sum(hemiInd) == 0, continue;end
        for i=1:numSeeds
            c=1;j=1;
            while j<=i && c<size(respSeq_elecInfo,1)
                if strcmp(respSeq_elecInfo{c,3}(1),hlbl); %check current channel hemisphere
                    if j==i,[~,cenInd]=ismember(respSeq_elecInfo{c,3},centroidLbls);end %if this is correct seed #, save centroid index
                    j=j+1;
                end
                c=c+1; % next channel
            end
            seed = centroids(cenInd,:); % centroid of first tone response
            d = ones(1,sum(hemiInd))*NaN;
            skipC =0;
            for chan=1:size(respSeq_elecInfo,1)
                if hemiInd(chan)
                    [~,cenInd]=ismember(respSeq_elecInfo{chan,3},centroidLbls);
                    d(chan-skipC)=sqrt((seed(1)-centroids(cenInd,1))^2+(seed(2)-centroids(cenInd,2))^2+(seed(3)-centroids(cenInd,3))^2);
                else
                    skipC = skipC+1;
                end
            end
            figure;
            scatter(d,t(respInd(hemiInd)));
            if hemi == 1
                title(['Left Hemisphere Tone Response. Origin Electrode: ' respSeq_elecInfo{c-1,3}]);
            else
                title(['Right Hemisphere Tone Response. Origin Electrode: ' respSeq_elecInfo{c-1,3}]);
            end
            xlabel('Distance From First Responding Electrode')
            ylabel('Tone Response latency (ms)')
            fit = polyfit(d,t(respInd(hemiInd)),1);
            yfit = polyval(fit,d);
            hold on
            h= plot(d,yfit,'-');
            yresid = t(respInd(hemiInd)) - yfit;
            SSresid = sum(yresid.^2);
            SStotal = (length(t(respInd(hemiInd)))-1) * var(t(respInd(hemiInd)));
            rsq = 1 - SSresid/SStotal;  %%% percentage of variance of Y explained by the linear model
            legend(h,['R^2: ' num2str(rsq)])
            if savePlots
                print(gcf,[saveDir,'/respDist_',respSeq_elecInfo{c-1,3},'.png'],'-dpng');
                close;
            end
        end
    end


% keyboard
%sample single channel
% figure
% plot(t,erp1(10,:))
% hold on
% plot(t,erp2(10,:))
% title(plot_title)
% legend(trialType1_label,trialType2_label)

%manually edit bad channels and resave file
%bad_channels = [];
%save(fullfile(ddir,'sessInfo.mat'),'elecInfo','ref','events','trial_type','trial_resp','t','srate','nlxEvents','bad_channels','chan_stats');