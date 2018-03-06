function anesthBH
% function ram_freq_cage
%   Compute behavioral metrics from events file.
% 
%  


% parameters

ddir = 'E:\data\C\data\eeg\'; % directory
subj = 'HUP119_i'; % subject
subjF = 'HUP119_i_24Jun16_0959'
% elec = 128; % Electrode number
tank = num2str(elec); % Electrode number 
fs = 1000; % sampling rate (Hz)
% nch = 1; % number of channels
% etime = 23.5; % length of time in minutes
% cor = [3.92 15.233]; % correct range
% LWM = 10.519; % LWM (loss of memory) in minutes from events file (lfpoffset / samplerate / 60) samplerate from params.txt
% LOC = 15.756; % LOC in minutes from events file
% sr = 2.047999e+03;

sr = eegparams('samplerate',[ddir subj filesep 'lfp.noreref' filesep 'params.txt']); %obtain sample rate from params text file
timefix = sr*60;

events = [];
load([ddir subj filesep 'behavioral' filesep 'session_0' filesep 'events.mat']);
run fixEEGlfpfile.m;
[correctR incorrectR] = eventStrategize(events); % find correct, incorrect, and pass (LOC) recall trials
cor = [correctR(1).lfpoffset/timefix correctR(size(correctR,2)).lfpoffset/timefix];
LWM = incorrectR(1).lfpoffset/timefix;
LOC = correctR(size(correctR,2)).lfpoffset/timefix; % estmated based on final correct trial
for n=1:length(events)
    events(n).lfpoffset2 = events(n).lfpoffset/timefix;
end

ci = 1;
% ii = 1;
% sr = 2.047999e+03 %sample rate from 

for i = 1:length(events)
    if i+2<length(events) && strcmp(events(i).type,'RECALL') && strcmp(events(i+2).type,'RECALL');
        events(i).correct = strcmp(events(i).response,events(i).target);
        events(i+1).correct = strcmp(events(i+1).response,events(i+1).target);
        events(i+2).correct = strcmp(events(i+2).response,events(i+2).target);
        events(i).percentcorrect = (events(i).correct + events(i+2).correct + events(i+1).correct)/3*100;
        events(i+1).percentcorrect = events(i).percentcorrect;
        events(i+2).percentcorrect = events(i).percentcorrect;
        pc(ci) = events(i);
        pc(ci+1) = events(i+1);
        pc(ci+2) = events(i+2);
        xpc(ci) = events(i).lfpoffset2;
        xpc(ci+1) = events(i+1).lfpoffset2;
        xpc(ci+2) = events(i+2).lfpoffset2;
        ypc(ci) = events(i).percentcorrect;
        ypc(ci+1) = events(i+1).percentcorrect;
        ypc(ci+2) = events(i+2).percentcorrect;
        ci = ci + 3;   

    end
end
% events = events(1); % take only the first event to look at entire time sequence.
% events.lfpfile = [ddir subj filesep 'lfp.noreref' filesep subjF];