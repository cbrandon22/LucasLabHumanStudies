% Manually check processed neuralynx signal (uses an_processSubj files and
% oddball_processing files)
clear
subj = 'HUP156_i';
ddir = fullfile('/Volumes/HumanStudies/HumanStudies/oddball/eeg',subj,'processed');


%% Calculate and plot ERPs
load(fullfile(ddir,'sessInfo.mat'));
files = dir(ddir);
for i=1:length(files)
    if ~ismember(files(i).name,{'.','..','.DS_Store','sessInfo.mat'})
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
bad_channels = [37,90];
save(fullfile(ddir,'sessInfo.mat'),'elecInfo','ref','events','trial_type','trial_resp','t','srate','nlxEvents','bad_channels','chan_stats');
%% Set dynamical criticality parameters
% possibly useful plots:
figure;
dat=look(events(1).lfpfile,5);%condensed single channel plot
[~,maxInd] = max(dat);
[~,minInd] = min(dat);
figure
plot(dat(maxInd-50000:maxInd+50000)); %plot samples around max value
figure
plot(dat(minInd-50:minInd+50));
% set and save parameters
parameters.ana = [5731201]; % anesthesia drug start/stop
parameters.coc = [15326453];
parameters.artifact = []; % Artifact in electrical signal
save(fullfile(ddir,'parameters.mat'),'parameters');