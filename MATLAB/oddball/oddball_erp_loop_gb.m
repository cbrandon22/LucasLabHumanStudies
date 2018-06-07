 %% Inputs 
clear;
% xlfile = fullfile('C:\Users\gregb\Documents\MATLAB\induction_awake');
% [~,~,xlcells] = xlsread(xlfile);
% for s = 1:length(xlcells)
% subj = xlcells{s,1};
% subjlist = {'HUP144_e','HUP145_e','HUP147_e','HUP148_e','HUP149_e','HUP151_e','HUP152_e','HUP153_i'};
subjlist = {'HUP156_i','HUP157_e','HUP159_e','HUP165_i','HUP166_i'};
for s=1:(length(subjlist))
subj = char(subjlist(s))
ddir = fullfile('D:\TNL_Data\oddball\eeg',subj,'processed');
trialType1 = {'TARGETLF'};
trialType2 = {'BACKGROUNDLF'}; % leave empty to only select trialType1
subtract_trialTypes = 1; %set to 1 to plot difference btw types
leadAvg_reref = 1; %re-reference to lead average
% Plot title/legend labels
plot_title = 'Target-background LF';
trialType1_label = 'Target';
trialType2_label = 'Background';
load([ddir '/sessInfo.mat']);
zscore = 1; %z-score electrodes to pre-trial baseline
saveFigsDir = fullfile('D:\figures',subj);
includeTrials = [1 size(trial_type,2)];

% includeTrials = strsplit(xlcells{s,2});
% includeTrials = [str2double(includeTrials{1}),str2double(includeTrials{2})];

%%
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
erp1 = ones(size(elecInfo,1)-length(bad_channels),size(t,2))*NaN;
if ~isempty(trialType2)
    erp2 = ones(size(elecInfo,1)-length(bad_channels),size(t,2))*NaN;
end
goodChan = 0;
erp_elecInfo = {};
files  = dir(ddir);
for f=1:length(files)
    if ~ismember(files(f).name,{'.','..','sessInfo.mat'})
        load(fullfile(ddir,files(f).name));
        includeChan = ~ismember(cell2mat(lead_elecInfo(:,1)),bad_channels);
        if leadAvg_reref
            lead_avg = nanmean(eeg(:,:,includeChan),3);
            eeg=bsxfun(@minus,eeg,lead_avg);
        end
        erp1(goodChan+1:goodChan+sum(includeChan),:) = squeeze(nanmean(eeg(includeTrials1,:,includeChan),1))';
        if ~isempty(trialType2)
            erp2(goodChan+1:goodChan+sum(includeChan),:) = squeeze(nanmean(eeg(includeTrials2,:,includeChan),1))';
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
    zerp1 = (erp1-erp1mean)./erp1std;
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
        caxis([-8,8])
    end
else
    if zscore==0
        h=imagesc(erp1);
    else
        h=imagesc(zerp1);
        caxis([-8,8])
    end
end
title([subj ' ' plot_title],'Interpreter','none');
c = colorbar;
if zscore == 0
    ylabel(c, 'Voltage')
else
    ylabel(c, 'Z-score')
    caxis([-8,8])
end
xlabel('Time');
ylabel('Channels');
xticklabels = -400:50:950;
xticks = linspace(1, size(erp1, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
set(gcf, 'Position', [100, 100, 1200, 700])
if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
cd(saveFigsDir);
print(gcf,[subj ' ' plot_title],'-dpng');
close;
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