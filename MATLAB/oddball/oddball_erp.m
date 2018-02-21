%% Inputs
subj = 'HUP159_e';
ddir = fullfile('/Volumes/HumanStudies/HumanStudies/oddball/eeg',subj,'processed');
trialType1 = {'TARGETHF','BACKGROUNDHF'};
trialType2 = {'TARGETLF','BACKGROUNDLF'}; % leave empty to only select trialType1
subtract_trialTypes = 0; %set to 1 to plot difference btw types
leadAvg_reref = 1; %re-reference to lead average
% Plot title/legend labels
plot_title = 'High vs Low Freq. ERPs';
trialType1_label = 'High Freq.';
trialType2_label = 'Low Freq.';

load([ddir '/sessInfo.mat']);
includeTrials = [725 1245];
%keyboard % manually set and run include trials based on nlxEvents

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
        erp_elecInfo = [erp_elecInfo;lead_elecInfo];
    end
end
%% Make Plots
keyboard
%sample single channel
figure
plot(t,erp1(10,:))
hold on
plot(t,erp2(10,:))
title(plot_title)
legend(trialType1_label,trialType2_label)

%manually edit bad channels and resave file
%bad_channels = [];
%save(fullfile(ddir,'sessInfo.mat'),'elecInfo','ref','events','trial_type','trial_resp','t','srate','nlxEvents','bad_channels','chan_stats');