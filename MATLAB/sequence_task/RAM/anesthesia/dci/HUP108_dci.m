%% HUP108_i
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
top_dir = '/Users/tnl/Desktop/C/data/eeg/'; % directory
subj = 'HUP108_i'; % subject
subjF = 'HUP108_i_06Nov15_0746'; % raw data name?
save_dir = [top_dir, subj, '/processed_data/dci'];
notes_dir = [top_dir, subj,  '/docs/'];
data_dir = [top_dir, subj, '/raw/2015-11-06_07-10-39'];

% about elecs
[~, labels, ~] = xlsread([notes_dir, 'electrodeMap.xlsx']); % not sure how cell naming scheme works, so read in everything then select what you want
% look for number of elecs labeled as 'grid' or 'strip' (not depth)
elecs = 1:69;

% about the data
nch = numel(elecs); % number of channels

% for timing and conversion
ev_srate = eegparams('samplerate',[top_dir subj filesep 'lfp.noreref' filesep 'params.txt']); %obtain sample rate from params text file
srate = 2048; % taken from header
%iso = []; % get this from .nev file

% events
load([top_dir, subj, '/behavioral/session_0/events.mat']); % get behavioral data

% filtering variables
hpf = 3;
lpf = []; % no low pass

%% Behavioral Data anlysis
% This will calculate the time window we wish to use, as well as the
% rection time percentage to data, which will be compared to the
% criticality index result

% load first electrode to get size of final matrix
file_name = [data_dir, '/CSC1.ncs'];
% use mode (last input) of 1 to get all points
[timestamps, ~, DAT, header] = getRawCSCData(file_name, 1, 2, 1);
data = zeros(nch, numel(DAT));

% load nev
load(['/Users/tnl/Dropbox/NEV/', subj, 'NEV.mat']);
% get sampling rate conversion factor: unclear what the sampling rate of
% lfpoffset is
samp_dur = NEV(end).lfpoffsetNEV - NEV(1).lfpoffsetNEV;
nl_srate = samp_dur/(size(data,2)/srate);

%Change of consciousness:
% find correct, incorrect, and pass (LOC) recall trials. Result is a subset
% of events
[correctR, incorrectR] = eventStrategize(events);
% change of consciousness, as indexed by the last correct trial
coc = (correctR(size(correctR,2)).lfpoffset/ev_srate); % in s
% add on timestamps from before the task started: coc +first pulse
%coc = floor(coc + (125250496/nl_srate))*srate; % HUP121
% coc = floor(coc + (1.927974984000000e09/nl_srate))*srate; % HUP 117
coc = floor(coc + (478944116/nl_srate))*srate; % HUP 108

% start of anasthesia:
% use the .NEV file to get from notes when anesthesia was administaered, or
% stopped (needs to be done manually, since notes are always recorded
% differently). Need to subtract timepoint for the start of recording
% start time - start of recording
%ana_start = 810058533-215667820; % determined by looking through NEV file
%ana_start = 158612032 - ; % HUP121i
%ana_start = 914527073 - 0; % HUP117i
ana_start = 903788448 - 0; %HUP 108
%ana_start = 302486416 - 45939498; % HUP 60
ana_start = floor(ana_start/nl_srate)*srate; % dont resample, get to Neuralynx srate

% save parameters for later
parameters.ana = ana_start;
parameters.coc = coc;
save([save_dir, '/parameters'], 'parameters')

%% Preprocess and concatenate

fprintf('\nElectrode');
for i = 1:nch
    fprintf('\n...%d', i);
    elec = elecs(i);
    % align to pulses - skipping alignment for now
    % I want the entire tsk recording, so I have selected events to be the
    % first event
    % DAT = an_getlfp_ms_wrapper(elec,events,end_pt);
    
    % load raw data
    file_name = [data_dir, '/CSC', num2str(elec), '.ncs'];
    % use mode (last input) of 1 to get all points
    [timestamps, ~, DAT, header] = getRawCSCData(file_name, 1, 2, 1);
    
    % filtfilt
    % high and low pass filters
    % order 4 is the standard
    [a,b] = butter(4, hpf/(srate/2), 'high'); % returns filter
    DAT = filtfilt(a, b, DAT); % does the filtering in both directions
    %   [a,b] = butter(4, lpf/(srate/2), 'low');
    %   DAT = filtfilt(a, b, DAT); % does the filtering in both directions
    % ideally, wouldn't do notch filters, but data is noisy
    %60 hz
    [a,b] = butter(4, [59/(srate/2) 61/(srate/2)], 'stop');
    DAT = filtfilt(a, b, DAT);
    %harmonics
    % [a,b] = butter(4, [119/(srate/2) 121/(srate/2)], 'stop');
    % DAT = filtfilt(a, b, DAT);
    
    
    %     % idealfilter
    %     % demean for idealfilter
    %      m = mean(DAT);
    %      DAT = DAT - m;
    %
    %     % make DAT a timeseries object
    %     DAT_ts = timeseries(DAT,(1:numel(DAT))./srate);
    %
    %     % high and low pass filters
    %     interval = [hpf, lpf];
    %     DAT_ts = idealfilter(DAT_ts, interval, 'pass'); % does the filtering in both directions
    %     % ideally, wouldn't do notch filters, but data is noisy
    %     %60 hz
    %     DAT_ts = idealfilter(DAT_ts, [59, 61], 'notch');
    %     %harmonics
    %     %DAT_ts = idealfilter(DAT_ts, [119, 121], 'notch');
    %     %DAT_ts = idealfilter(DAT_ts, [179, 181], 'notch');
    %
    %     % make array again
    %     DAT = DAT_ts.data;
    %
    %     % add mean back
    %     DAT = DAT + m;
    
    % save
    %save([save_dir, '/', num2str(elec), '_filt'], 'DAT', '-v7.3');
    data(i,:) = DAT;
end

%% Look for noisey electrodes in data

eegplot(data, 'srate', srate);

%% Clean data

% Create elec structure
% add noisy or not recording
elec_info.all = elecs;
elec_info.bad = sort([62 45]);
% logical indexing for good elecs
idx = true(size(elecs));
idx(elec_info.bad) = 0;
elec_info.good = elecs(idx);

data_good = create_CAR(subj, numel(elec_info.all), elec_info.bad, 'gdat', data);

data_good = data_good(idx,:);

% FFT to look for line noise
% spectopo(data, frames, srate)
% some electrodes still show noise in the FFT. Going to leave it for now
spectopo(data_good, 0, srate);


%%
% if everything looks good, save final data matrix and remove bad
% electrodes

save([save_dir, '/data'], 'data_good', '-v7.3');
save([save_dir, '/elecs'], 'elec_info', '-v7.3');
save([save_dir, '/srate'], 'srate');

%%
for i = 1: numel(NEV)
    nev_time(i) = NEV(i).timestamp;
end