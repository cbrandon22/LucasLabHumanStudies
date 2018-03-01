% Events File export
ddir = '/Users/tnl/Desktop/C/data/eeg/HUP138_i/processed/'; % directory 'E:\data\C\data\eeg\' 'L:\C\data\'
subj = 'HUP132_i'; % subject 'HUP117_i'
%[subjInfo,dataDir] = an_subjInfo(subj);
cd([ddir subj filesep 'raw/2017-04-21_07-04-47']);  %L:\C\data\HUP117_i\raw\2016-06-08_17-02-24\

[timestamps, EventIDs, TTLs, Extras, EventStrings, Header] = Nlx2MatEV('Events.nev',[1 1 1 1 1],1,1,[]);

lfpDir = [ddir subj filesep 'lfp.noreref'];
cd(lfpDir);
fname = ls('*.001');
[~,fileExt,~] = fileparts(fname); % automatically obtain the correct filename for electrodes

% lfp_samplerate = getNlxJackSheetInfo(lfpDir,fileExt,'DOWNSAMPLE_RATE','all');

% nev_samplerate = getNlxJackSheetInfo(lfpDir,fileExt,'ORIG_SAMPLERATE','all');
% first_sample_time = getNlxJackSheetInfo(lfpDir,fileExt,'FIRSTSAMPLETIME','all');
% nlxEvents_times = [ones(length(timestamps),1) timestamps]*b;
    
   for k=1:length(timestamps) % timestamps of important annotated text (e.g. drug delivery, patient response to a question, etc)
      thisNlxEvTime = timestamps(k);
      NEV(k).timestamp = thisNlxEvTime;
%       NEV(k).lfpoffsetNEV = round((thisNlxEvTime-(first_sample_time/1000))*(lfp_samplerate/1000)); %convert to relative offset in the LFP file, basically an index of LFP sample - currently not working...
      NEV(k).lfpoffsetNEV = thisNlxEvTime-NEV(1).timestamp; % conversion to generate microseconds
      NEV(k).note = EventStrings(k);
      NEV(k).TTL = TTLs(k);
      NEV(k).EventID = EventIDs(k);
    end
NEV(1).HeaderInfo = Header;
%cd_mkdir([ddir filesep subj filesep 'processed']);
save(['/Users/tnl/Desktop/C/data/eeg/HUP138_i/', subj, 'NEV.mat'],'NEV');
n = find(EventIDs == 4);
t = timestamps(n)';
notes = EventStrings(n);
% (t(20)-t(1))/2.59769e6/60 % calculate times to roughly correlate times


% thisNlxEvtime = timestamp_from_NEV; % each clinically relevant timestamp from .nev files
% lfpoffset = round((thisNlxEvTime-(first_sample_time/1000))*(lfp_samplerate/1000)); % convernt to relative offset in the LFP file