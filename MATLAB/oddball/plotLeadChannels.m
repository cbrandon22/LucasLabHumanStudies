% Manually check processed neuralynx signal (uses an_processSubj files and
% oddball_processing files)
clear
subj = 'HUP155_i';
ddir = fullfile('/Volumes/HumanStudies/HumanStudies/oddball/eeg',subj,'processed');


%% Calculate and plot ERPs
load(fullfile(ddir,'sessInfo.mat'));
files = dir(ddir);
for i=1:length(files)
    if ~ismember(files(i).name,{'.','..','sessInfo.mat'})
        load(fullfile(ddir,files(i).name));
        leadChannels = cell2mat(lead_elecInfo(:,1));
        leadEEG = [];
        for j=1:length(leadChannels)
            leadEEG(j,:) = look(events(1).lfpfile,leadChannels(j),[],1)';
        end
        eegplot(leadEEG,'srate',srate,'title',files(i).name);
    end
end
keyboard
% Manually fix bad_channels and re-save sessInfo.mat if altered
bad_channels = [109,79];
save(fullfile(ddir,'sessInfo.mat'),'elecInfo','ref','events','trial_type','trial_resp','t','srate','nlxEvents','bad_channels','chan_stats');
for i=2:size(chan_stats,2)
    figure;
    hist(chan_stats(:,i),20)
end