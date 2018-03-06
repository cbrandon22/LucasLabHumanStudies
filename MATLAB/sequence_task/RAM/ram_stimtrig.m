function ram_stimtrig(varargin)
% function ram_stimtrig(varargin)
%   Find indices of stim train onsets and save with other information about
%   the stimulus train.
%
%   DR 03/2015

% parameters
ddir = 'E:\'; % data directory
subj = 'fidel'; % subject
tank = 'fidel_VPLT_DT1_071015'; % data tank
block = 1; % data block
offset = 0.6; % train time adjustment to account for delay stim slew rate (ms)
sfreq = 60; % train frequency (Hz)
spuls = 20; % train number of pulses
fs = 24414.06; % sampling rate (Hz)
if nargin
    v2struct(varargin{1});
end

% find stim triggers (on ch25)
cd([ddir subj '\' tank '\Block-' num2str(block)]);
fid = fopen([tank '_Block-' num2str(block) '_xRAW_ch25.sev'],'r');
fseek(fid,40,'bof');
dat = fread(fid,inf,'single');
fclose(fid);
itrain = find(dat(1:end-1)<0.5 & dat(2:end)>0.5); % trigger for each pulse train
itrain = itrain + round(offset/1000*fs);

% save
save([ddir subj '\processed\' tank '_' num2str(block) '.mat'],'itrain','sfreq','spuls','fs','-mat');