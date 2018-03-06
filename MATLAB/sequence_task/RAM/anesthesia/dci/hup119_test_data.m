
%% Initialize and load behaviorally relavant variables

% for saving and loading
top_dir = '/Users/tnl/Desktop/C/data/eeg/'; % directory
subj = 'HUP119_i'; % subject
subjF = 'HUP119_i_24Jun16_0959'; % raw data name?
save_dir = [top_dir, subj, '/processed_data/dci'];
notes_dir = [top_dir, subj,  '/docs/'];
data_dir = [top_dir, subj, '/raw/2016-06-24_09-56-00'];

% about elecs
[~, labels, ~] = xlsread([notes_dir, 'electrodeMap.xlsx']); % not sure how cell naming scheme works, so read in everything then select what you want
% look for number of elecs labeled as 'grid' or 'strip' (not depth)
elecs = 17:32;
    
% about the data
nch = numel(elecs); % number of channels

% for timing and conversion
ev_srate = eegparams('samplerate',[top_dir subj filesep 'lfp.noreref' filesep 'params.txt']); %obtain sample rate from params text file
srate = round(ev_srate);
%iso = []; % get this from .nev file

% events - not relevant when loadaing raw
% events = [];
% load([ddir subj filesep 'behavioral' filesep 'session_0' filesep 'events.mat']); % get behavioral data
% fixEEGlfpfile_MAC % use this when working from the MAC
% events = events(1); % take only the first event to look at entire time sequence.



%% Preprocess and concatenate

HUP119.srate = 2048;
HUP119.preproc = [{'60 Hz line noise Butterworth filter'}, {'Common Average Referencing'}];
HUP119.layout = {'Two 8 electrode strips over IFG/inferior sensorimotor. Elec 1 is most superior and anterios, elec 16 is most inferior and posterior. First and last electrode on second strip were noisy, and removed'};
% load first electrode to get size of final matrix
file_name = [data_dir, '/CSC1.ncs'];
% use mode (last input) of 1 to get all points
[timestamps, DAT, header] = getRawCSCData(file_name, 1, 2, 1);
HUP119.data = zeros(nch, numel(DAT));
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
    [timestamps, DAT, header] = getRawCSCData(file_name, 1, 2, 1);

    % high and low pass filters
    % order 4 is the standard
    %[a,b] = butter(4, [hpf/(srate/2), lpf/(srate/2)]); % returns filter
    %DAT = filtfilt(a, b, DAT); % does the filtering in both directions
    % ideally, wouldn't do notch filters, but data is noisy
    %60 hz
    [a,b] = butter(4, [59/(srate/2) 61/(srate/2)], 'stop');
    DAT = filtfilt(a, b, DAT);
    [a,b] = butter(4, [119/(srate/2) 121/(srate/2)], 'stop');
    DAT = filtfilt(a, b, DAT);
    [a,b] = butter(4, [179/(srate/2) 181/(srate/2)], 'stop');
    DAT = filtfilt(a, b, DAT);

    % save
    %save([save_dir, '/', num2str(elec), '_filt'], 'DAT', '-v7.3');
    HUP119.data(i,:) = DAT;
end

%% Clean data
data = HUP119.data;
% Create elec structure
% add noisy or not recording
elec_info.all = elecs;
elec_info.bad = sort([9 16]);
% logical indexing for good elecs
idx = true(size(elecs));
idx(elec_info.bad) = 0;
elec_info.good = elecs(idx);

% cut out arftifact
data_good = data(:, 1:2693120);
data_good = data_good(:, end-60001:end); % get shorter time

% CAR
data_good = create_CAR(subj, numel(elec_info.all), elec_info.bad, 'gdat', data_good);

data_good = data_good(idx,:);

% FFT to look for line noise
% spectopo(data, frames, srate)
% some electrodes still show noise in the FFT. Going to leave it for now
spectopo(data_good, 0, srate)
%% SAve
HUP119.data = data_good;
save('/Users/tnl/Desktop/C/data/eeg/HUP119_i/processed_data/tcn_human_ecog', 'HUP119');

