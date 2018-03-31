 %% Inputs 
clear;
subj = 'HUP149_e';
%ddir = fullfile('D:\TNL_Data\oddball\eeg',subj,'processed');
ddir = fullfile('/Volumes/HumanStudies/HumanStudies/oddball/eeg',subj,'processed');
trialType1 = {'TARGETHF','TARGETLF','BACKGROUNDHF','BACKGROUNDLF'};
trialType2 = {}; % leave empty to only select trialType1
subtract_trialTypes = 0; %set to 1 to plot difference btw types
leadAvg_reref = 1; %re-reference to lead average
zscore = 1; %z-score electrodes to pre-trial baseline

% Use these to check different points in data manipulation
plotLeadSingleTrials = 1; % make single trial plots for each lead
overlayReref = 1; % plot rereference result on top of raw signal
plotLeadAvgTrials = 1; % plot the average across trials
zScoreDemoPlots = 1; %single channel z score subplot demonstration

% Plot title/legend labels
plot_title = [];
trialType1_label = 'All Tones';
trialType2_label = [];
load([ddir '/sessInfo.mat']);
if exist(fullfile(ddir,'manual_bad_channels.mat'),'file')==2
    load(fullfile(ddir,'manual_bad_channels.mat'),'bad_channels')
end
includeTrials = [1 size(trial_type,2)];
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
    if ~ismember(files(f).name,{'.','..','sessInfo.mat'})
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
        if leadAvg_reref
            lead_avg = nanmean(eeg(:,:,includeChan),3);
            reeg=bsxfun(@minus,eeg,lead_avg);
            rerp1(goodChan+1:goodChan+sum(includeChan),:) = squeeze(nanmean(reeg(includeTrials1,:,includeChan),1))';
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