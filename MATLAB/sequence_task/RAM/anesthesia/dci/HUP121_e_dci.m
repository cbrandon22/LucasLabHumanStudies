%% HUP121_i
% This script will perform preprocessing analyses to get parameters ideal
% for AR model and eigenvalue analysis
%
% This includes aligning to the task, filtering between 3 and 300 Hz
%
%
%% Initialize and load behaviorally relavant variables

% add neurolynx stuff to path
addpath(genpath('/Users/tnl/Desktop/MATLAB/releaseDec2015/'));
addpath(genpath('/Users/tnl/Desktop/MATLAB/RAM/anesthesia/old_scripts_stiso/'));

% for saving and loading
top_dir = '/Users/tnl/Desktop/C/data/eeg/'; % directory
subj = 'HUP121_e'; % subject
subjF = 'HUP121_e_01Jul16_0822'; % raw data name?
save_dir = [top_dir, subj, '/processed_data/dci'];
notes_dir = [top_dir, subj,  '/docs/'];
data_dir = [top_dir, subj, '/raw/2016-06-22_13-35-07'];

% about elecs
[~, labels, ~] = xlsread([notes_dir, 'electrodeMap(missing2).xlsx']); % not sure how cell naming scheme works, so read in everything then select what you want
% now taking depths
elecs = 1:115;

% about the data
nch = numel(elecs); % number of channels

% for timing and conversion
ev_srate = eegparams('samplerate',[top_dir subj filesep 'lfp.noreref' filesep 'params.txt']); %obtain sample rate from params text file
srate = 2048; % taken from header
%iso = []; % get this from .nev file

% events
load([top_dir, subj, '/behavioral/session_0/events.mat']); % get behavioral data
fixEEGlfpfile_MAC % use this when working from the MAC
events_last = events(1);
dur = events(end).lfpoffset/srate*1000;

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
[~, ~, DAT, ~] = getRawCSCData(file_name, 1, 2, 1);

% load nev
%load(['/Users/tnl/Dropbox/NEV/', subj, 'NEV.mat']);
% get neuralynx srate
% duration of exp in seconds: taken from header of NEV file
%nl_dur = ((8*60*60) + (36*60) + 26) - ((8*60*60) + (20*60) + 3);% srate = total samples/ total seconds
%nl_srate = (NEV(end).lfpoffsetNEV - NEV(1).lfpoffsetNEV)/nl_dur;

%Change of consciousness:
% find correct, incorrect, and pass (LOC) recall trials. Result is a subset
% of events
%[correctR, incorrectR] = eventStrategize(events);
% change of consciousness, as indexed by the last correct trial
%coc = (correctR(size(correctR,2)).lfpoffset/ev_srate); % in s
% add on timestamps from before the task started: coc +first pulse
coc = 1; 

% start of anasthesia:
% use the .NEV file to get from notes when anesthesia was administaered, or
% stopped (needs to be done manually, since notes are always recorded
% differently). Need to subtract timepoint for the start of recording
% start time - start of recording
%ana_start = 810058533-215667820; % determined by looking through NEV file
%ana_start = 255155420/nl_srate; %- 125250496; % HUP121i
%ana_start = 2.085477936000000e9 - 1.927974984000000e9; % HUP117, using propofol increase
%ana_start = 903788448 - 478944116; %HUP 108
%ana_start = 302486416 - 45939498; % HUP 60
%ana_start = floor(ana_start)*srate; % get to Neuralynx srate

% save parameters for later
%parameters.ana = ana_start;
parameters.coc = coc;
save([save_dir, '/parameters'], 'parameters')

% visualize coc and ana
figure(1)
clf;
plot(DAT)
hold on
plot([coc,coc], [min(DAT), max(DAT)], 'r')
%plot([ana_start,ana_start], [min(DAT), max(DAT)], 'k')

%% Preprocess and concatenate

close
% load first electrode to get size of final matrix
file_name = [data_dir, '/CSC1.ncs'];
% use mode (last input) of 1 to get all points
[~, ~, DAT, header] = getRawCSCData(file_name, 1, 2, 1);
data = zeros(nch, numel(DAT));
for i = 1:nch
    fprintf('\n%d', i)
    elec = elecs(i);
    % align to pulses - skipping alignment for now
    % I want the entire tsk recording, so I have selected events to be the
    % first event
    % DAT = an_getlfp_ms_wrapper(elec,events,end_pt);
    
    % load raw data
    file_name = [data_dir, '/CSC', num2str(elec), '.ncs'];
    % use mode (last input) of 1 to get all points
    [~, ~, DAT, header] = getRawCSCData(file_name, 1, 2, 1);

    % high and low pass filters
    % order 4 is the standard
    [a,b] = butter(4, hpf/(srate/2), 'high'); % returns filter
    DAT = filtfilt(a, b, DAT); % does the filtering in both directions
    % ideally, wouldn't do notch filters, but data is noisy
    %60 hz
    [a,b] = butter(4, [59/(srate/2) 61/(srate/2)], 'stop');
    DAT = filtfilt(a, b, DAT);
%     %harmonics
%     [a,b] = butter(4, [119/(srate/2) 121/(srate/2)], 'stop');
%     DAT = filtfilt(a, b, DAT);
%     [a,b] = butter(4, [179/(srate/2) 181/(srate/2)], 'stop');
%     DAT = filtfilt(a, b, DAT);
%     [a,b] = butter(4, [239/(srate/2) 241/(srate/2)], 'stop');
%     DAT = filtfilt(a, b, DAT);

    % save
    %save([save_dir, '/', num2str(elec), '_filt'], 'DAT', '-v7.3');
    data(i,:) = DAT;
end
%% Look for noisey electrodes in data

eegplot(data, 'srate', srate);

% save the srate
save([save_dir, '/srate'], 'srate');

%% Clean data

% Create elec structure
% add noisy or not recording
elec_info.all = elecs;
elec_info.bad = sort([40 42 68 83 97:115]);
% logical indexing for good elecs
idx = true(size(elecs));
idx(elec_info.bad) = 0;
elec_info.good = elecs(idx);

% CAR
data_good = create_CAR(subj, numel(elec_info.all), elec_info.bad, 'gdat', data);

data_good = data_good(idx,:);

% FFT to look for line noise
% spectopo(data, frames, srate)
% some electrodes still show noise in the FFT. Going to leave it for now
spectopo(data_good, 0, srate);
%%
data = data_good;
%%
% if everything looks good, common average reference and save data and
% electrodes

save([save_dir, '/data'], 'data_good', '-v7.3');
save([save_dir, '/elecs'], 'elec_info', '-v7.3');