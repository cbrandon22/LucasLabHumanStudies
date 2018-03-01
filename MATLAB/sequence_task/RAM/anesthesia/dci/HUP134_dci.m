%% HUP121_e
% This script will perform preprocessing analyses to get parameters ideal
% for AR model and eigenvalue analysis
%
% This includes aligning to the task, filtering between 3 and 200 Hz
%
% This subject does not have enough data points for this analysis to be
% worth it
%
%% Initialize and load behaviorally relavant variables

% for saving and loading
ddir = '/Users/tnl/Desktop/C/data/eeg/'; % directory
subj = 'HUP134_e'; % subject
subjF = 'HUP134_e_06Feb17_1651'; % raw data name?
save_dir = ['/Users/tnl/Desktop/C/data/eeg/', subj, '/processed_data/dci'];
elecs = 1:128;

% about the data
srate = 1000; % sampling rate (Hz) fter being resampled
nch = numel(elecs); % number of channels

% for timing and conversion
ev_srate = eegparams('samplerate',[ddir subj filesep 'lfp.noreref' filesep 'params.txt']); %obtain sample rate from params text file
%timefix = sr*60; % conversion to minutes from time points?
iso = []; % get this from .nev file

load([ddir subj filesep 'behavioral' filesep 'session_0' filesep 'events.mat']); % get behavioral data
fixEEGlfpfile_MAC % use this when working from the MAC
end_pt = floor(1e6*(srate/ev_srate)); % length of data in ms, determined with trial and error. Find better way to do this later
events = events(1); % take only the first event to look at entire time sequence.

% filtering variables
hpf = 3;
lpf = 180;

%% Preprocess and concatenate

% load first electrode to get size of final matrix
data = zeros(nch, end_pt);
for i = 1:nch
    elec = elecs(i);
    % align to pulses
    % I want the entire tsk recording, so I have selected events to be the
    % first event
    DAT = an_getlfp_ms_wrapper(elec,events,end_pt);

    % high and low pass filters
    % order 4 is the standard
    [a,b] = butter(4, [hpf/(srate/2), lpf/(srate/2)]); % returns filter
    DAT = filtfilt(a, b, DAT); % does the filtering in both directions
    % this ws taken out since noise shouldn't effect the model
    %60 hz
    [a,b] = butter(4, [59/(srate/2) 61/(srate/2)], 'stop');
    DAT = filtfilt(a, b, DAT);
    %harmonics
    [a,b] = butter(4, [119/(srate/2) 121/(srate/2)], 'stop');
    DAT = filtfilt(a, b, DAT);
    [a,b] = butter(4, [179/(srate/2) 181/(srate/2)], 'stop');
    DAT = filtfilt(a, b, DAT);
    [a,b] = butter(4, [239/(srate/2) 241/(srate/2)], 'stop');
    DAT = filtfilt(a, b, DAT);
    
    % resample to seconds - decided against this, changed freq
    % realtionships
    % DAT = resample(DAT,1,srate);
    
    % save
    save([save_dir, '/', num2str(elec), '_filt'], 'DAT', '-v7.3');
    data(i,:) = DAT;
end

%% Look for noisey electrodes in data

eegplot(data, 'srate', srate);

%% Clean data

% Create elec structure
% add noisy or not recording
elec_info.all = elecs;
elec_info.bad = sort([17:22 39 53]);
% logical indexing for good elecs
idx = true(size(elecs));
idx(elec_info.bad) = 0;
elec_info.good = elecs(idx);

data_good = data(idx,:);

%% Look for noisey electrodes in data

eegplot(data_good, 'srate', srate);

%% FFT to look for line noise
% spectopo(data, frames, srate)
% some electrodes still show noise in the FFT. Going to leave it for now
spectopo(data_good, 0, srate)

%%
% if everything looks good, save final data matrix and remove bad
% electrodes

save([save_dir, '/data'], 'data_good', '-v7.3');
save([save_dir, '/elecs'], 'elec_info', '-v7.3');