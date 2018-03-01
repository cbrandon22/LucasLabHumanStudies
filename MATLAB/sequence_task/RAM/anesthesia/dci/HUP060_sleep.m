%% HUP060
%
% sleep data


%% define constants
pwd_file = '/Users/tnl/Desktop/MATLAB/RAM/jas_ieeglogin.bin';
username = 'jastiso';
dataset = 'HUP060_phaseIV';
subj = 'HUP060_i';
% open a sesssion: name of dataset, username, password file
session = IEEGSession(dataset, username, pwd_file);
srate = session.data.sampleRate;

labels = cell(size(session.data.channelLabels,1),1);
% get labels
for i = 1:size(session.data.channelLabels,1)
    labels{i} = session.data.rawChannels(1,i).label;
end

% filtering variables
% filtering variables
hpf = 3;
lpf = []; % no low pass

%% Get data off portal

% get annotations - returns markers in microseconds
anns = session.data.annLayer.getEvents(0);
% after manual servey of events, I decided I want to take data from the one
% marked drowsy? to REMish, including some slow
events_idx = [5 6 8];
events = struct();
for i = 1:numel(events_idx)
    idx = events_idx(i);
    events(i).ev_type = anns(idx).description;
    events(i).start_us = anns(idx).start;
    events(i).start_samp = (anns(idx).start/(10^6))./srate;
    events(i).stop_us = anns(idx).stop;
    events(i).stop_samp = (anns(idx).stop/(10^6))./srate;
end
dur_drowsy = events(1).stop_us - events(1).start_us;
dur_slow = events(2).stop_us - events(2).start_us;
dur_rem = events(3).stop_us - events(3).start_us;
pad = 10*10^6; % s padding on either side

% get values: time, elecs
data_drowsy = session.data.getvalues(events(1).start_us - pad, dur_drowsy + 2*pad, 1:size(session.data.channelLabels,1))';
data_slow = session.data.getvalues(events(2).start_us - pad, dur_slow + 2*pad, 1:size(session.data.channelLabels,1))';
data_rem = session.data.getvalues(events(3).start_us - pad, dur_rem + 2*pad, 1:size(session.data.channelLabels,1))';

%% preprocess

eegplot(data_drowsy, 'srate', srate);
eegplot(data_rem, 'srate', srate);
eegplot(data_slow, 'srate', srate);

% elecs
% Create elec structure
% add noisy or not recording
elec_info.all = 1:size(session.data.channelLabels,1);
elec_info.bad = sort([1:12 17:21 22 25 56:59]);
% logical indexing for good elecs
idx = true(size(elec_info.all));
idx(elec_info.bad) = 0;
elec_info.good = elec_info.all(idx);

% high and low pass filters
% order 4 is the standard
for i = 1:size(data_rem,1)
    [a,b] = butter(2, hpf/(srate/2), 'high'); % returns filter
    data_drowsy(i,:) = filtfilt(a, b, data_drowsy(i,:)); % does the filtering in both directions
    data_slow(i,:) = filtfilt(a, b, data_slow(i,:));
    data_rem(i,:) = filtfilt(a, b, data_rem(i,:));
    
    % ideally, wouldn't do notch filters, but data is noisy
    %60 hz & harmonics
    [a,b] = butter(2, [59/(srate/2) 61/(srate/2)], 'stop');
    data_drowsy(i,:) = filtfilt(a, b, data_drowsy(i,:)); % does the filtering in both directions
    data_slow(i,:) = filtfilt(a, b, data_slow(i,:));
    data_rem(i,:) = filtfilt(a, b, data_rem(i,:));
    
    [a,b] = butter(2, [179/(srate/2) 181/(srate/2)], 'stop');
    data_drowsy(i,:) = filtfilt(a, b, data_drowsy(i,:)); % does the filtering in both directions
    data_slow(i,:) = filtfilt(a, b, data_slow(i,:));
    data_rem(i,:) = filtfilt(a, b, data_rem(i,:));
end

% FFT to look for line noise
% spectopo(data, frames, srate)
% some electrodes still show noise in the FFT. Going to leave it for now
figure(1)
spectopo(data_drowsy, 0, srate);
figure(2)
spectopo(data_slow, 0, srate);
figure(3)
spectopo(data_rem, 0, srate);

% CAR
data_r_good = create_CAR(subj, numel(elec_info.all), elec_info.bad, 'gdat', data_rem);
data_s_good = create_CAR(subj, numel(elec_info.all), elec_info.bad, 'gdat', data_slow);
data_d_good = create_CAR(subj, numel(elec_info.all), elec_info.bad, 'gdat', data_drowsy);

% remove bad elecs
data_r_good = data_r_good(idx, :);
data_s_good = data_s_good(idx, :);
data_d_good = data_d_good(idx, :);

%% save things

save_dir = ['/Users/tnl/Desktop/C/data/eeg/', subj, '/processed_data/sleep/'];

save([save_dir, 'srate'], 'srate')
save([save_dir, 'events'], 'events')
save([save_dir, 'elecs'], 'elec_info')
save([save_dir, 'data_drowsy'], 'data_d_good')
save([save_dir, 'data_rem'], 'data_r_good')
save([save_dir, 'data_slow'], 'data_s_good')

