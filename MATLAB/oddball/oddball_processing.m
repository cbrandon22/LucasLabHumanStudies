function oddball_processing(subj)
%This script takes a look at raw data quality from the neuralynx
%Plots eeg for all channels to find bad channels or bad time windows
%% Inputs
task = 'oddball';
dirs = le_dirs(task);
saveDir = fullfile(dirs.data,'eeg',subj,'processed');
filt_ln = 1; % filter 60Hz and harmonics to 300Hz
check_bad_channels = 0; % double check auto assignment of bad channels (all channels will be saved regardless)

%% load nlx events and select timepoints of interest
% nlxEvents are manually typed into the neuralynx during a recording.
% Although best practice is to mark similar timepoints in the same way, it
% is subject to variance and typos. We must open the struct and select
% events of interest (drug delivery, change of consciousness (coc), etc.

load(fullfile(dirs.data,'eeg',subj,'behavioral/session_0/nlxEvents.mat'),'nlxEvents');
load(fullfile(dirs.data,'eeg',subj,'behavioral/session_0/events.mat'),'events');

% behavioral analysis
subjSplit = strsplit(subj,'_');
if strcmp(task,'oddball') && ismember('GO',{events.type}) %check for task response data
    for i=1:length(events)
        if isempty(events(i).response)
            events(i).response = NaN;
        end
        if isempty(events(i).trial)
            events(i).trial= NaN;
        end
    end
    goCorrect = strcmp({events.type},'GO') & cell2mat({events.response})==1;
    goIncorrect = strcmp({events.type},'GO') & cell2mat({events.response})==0;
    winMSdur = 30000; %30 second sliding window
    winStep = 1000; %1 second step size
    winStart = events(1).mstime;
    winEnd = winStart+winMSdur;
    while winEnd < events(end).mstime
        events.mstime
    end
end

% keyboard % manually select events from nlxEvents
% ana = 'custom'; %select index of anesthesia start(induction)/stop(emergence)
% coc = 'custom'; %select index of behavioral change of consciousness 

%% Electrophysiological data organization
lfp_fileroot = events(1).lfpfile;
% load jacksheet
fid = fopen([lfp_fileroot,'.JacksheetAndParams.txt'],'r');
jacksheet=strsplit(fgetl(fid));
tline = fgetl(fid);
jline = 2;
while ischar(tline)
    jacksheet(jline,:) = strsplit(tline);
    tline=fgetl(fid);
    jline=jline+1;
end
fclose(fid);
jacksheet = jacksheet(:,1:end-1);
[~,srateInd] = ismember('DOWNSAMPLE_RATE',jacksheet(1,:));
chan_srate = str2double(jacksheet(2:end,srateInd));
if length(unique(chan_srate))>1
    error('sample rate is not consistent for all channels')
end
srate = round(unique(chan_srate));
[~,csc_nameInd] = ismember('channel_name',jacksheet(1,:));
[~,chan_numInd] = ismember('channel_num',jacksheet(1,:));

% Get electrode names
elecMapFn = fullfile(dirs.data,'eeg',subj,'docs/electrodeMap.xlsx');
[~,~,xlcells] = xlsread(elecMapFn);
emptyChannels = [];
elecInfo = {};
for i=1:length(xlcells)
    if isnan(xlcells{i,3})
        emptyChannels = [emptyChannels, i];
    elseif i<129
        [~,cInd] = ismember([xlcells{i,1} '.ncs'],jacksheet(:,csc_nameInd));
        elecInfo = [elecInfo;{str2double(jacksheet{cInd,chan_numInd})}, xlcells{i,[1 3]}];
    end
end
[~,rind]=ismember(xlcells{129,1},elecInfo(:,2));
ref = elecInfo(rind,:);
channels = cell2mat(elecInfo(:,1));
% add anatomical info if available
anat_fn = fullfile(dirs.data(1:strfind(dirs.data,task)-1),'localization',subjSplit{1},'post/electrodenames_native.csv');
if exist(anat_fn)==2
    fid = fopen(anat_fn);
    tline = fgetl(fid);
    anat_lbls = {};
    while ischar(tline)
        anat_lbls = [anat_lbls;strsplit(tline,',')];
        tline = fgetl(fid);
    end
    for i=1:length(anat_lbls);
        [lia,locb]=ismember(anat_lbls{i,1},{elecInfo{:,3}});
        if lia,elecInfo(locb,4) = {strtrim(anat_lbls{i,2})};end
    end
    missing_anat = false;
else
    missing_anat = true;
end
%create eeg matrix channel x trial ms x trial number
trial_start = -400;
trial_end = 950;
stim_ev_types = {'BACKGROUNDLF','BACKGROUNDHF','TARGETLF','TARGETHF','GO','NOGO'};
stim_events = ismember({events.type},stim_ev_types);
if sum(stim_events(end-2:end))>0 %last trial was cut short, ignore
    stim_events(end-2:end) = false;
end
trial_type = {events(stim_events).type};
trial_resp = cell2mat({events(stim_events).response});
rt=[];
for i=1:length(stim_events) %calculate reaction times
    if stim_events(i)
        rt = [rt,events(i+2).mstime-events(i).mstime];
    end
end
trial_resp = [trial_resp;rt];
for i=1:length(nlxEvents)
    [~, index]=min(abs(cell2mat({events.lfpoffset})-nlxEvents(i).lfpoffset));
    nlxEvents(i).approxTrial = events(index).trial;
end
leads = regexp({elecInfo{:,3}},'\D*','match');
leads = vertcat(leads{:});
lead_prefixes = unique(leads);
bad_channels = [];
chan_stats = [];
for l=1:length(lead_prefixes)
    prefix = lead_prefixes{l};
    disp(prefix);
    mat_fn = [subj '_' prefix '.mat'];
    lead_chs = channels(strcmp(prefix,leads));
    eeg = ones(sum(stim_events),ceil((trial_end-trial_start)/1000*srate),length(lead_chs))*NaN;
    for i=1:length(lead_chs)
        elecNum = lead_chs(i);
        eeg(:,:,i) = an_getlfp_ms_wrapper(elecNum,events(stim_events),trial_end-trial_start,trial_start);
        dat = look(lfp_fileroot,elecNum,[],1)';
        [M,SD,sk,k,med,zlow,zhi,tM,tSD,~,ksh,ksp,ksstat,kscv]=signalstat_cb(dat,0);
        chan_stats = [chan_stats; elecNum,M,SD,sk,k,med,zlow,zhi,tM,tSD,ksh,ksp,ksstat,kscv];
        if ksstat>0.1,bad_channels = [bad_channels,elecNum];end
    end
    if filt_ln
        for c=1:size(eeg,3)
            for ii=1:size(eeg,1)
                if sum(isnan(eeg(ii,:,c)))<5;
                    ln = mtmlinenoise(eeg(ii,:,c),3,srate,srate,60:60:300)';
                    eeg(ii,:,c) = eeg(ii,:,c)-ln;
                end
            end
        end
    end
    lead_elecInfo = elecInfo(strcmp(prefix,leads),:);
    t = linspace(trial_start,trial_end,ceil((trial_end-trial_start)/1000*srate));
    if ~(exist(saveDir)==7),mkdir(saveDir);end
    save(fullfile(saveDir,mat_fn),'eeg','lead_elecInfo','-v7.3');
end
save(fullfile(saveDir,'sessInfo.mat'),'elecInfo','jacksheet','ref','events','trial_type','trial_resp','t','srate','nlxEvents','bad_channels','chan_stats');
disp('done')
if missing_anat,disp('Anatomical info excluded');end
if check_bad_channels
%     all_eeg = [];
%     for i=1:length(channels)
%         all_eeg(i,:) = look(lfp_fileroot,channels(i),[],1)';
%     end
%     eegplot(all_eeg,'srate',srate,'title','All Channels','dispchans',64);
%     disp('Manually Edit bad_channels, then run final 2 lines.');
%     keyboard;
%     bad_channels = [];
%     save(fullfile(saveDir,'sessInfo.mat'),'elecInfo','ref','events','trial_type','trial_resp','t','srate','nlxEvents','bad_channels','chan_stats');

    
    bad_eeg = [];
    for i=1:length(bad_channels)
        bad_eeg(i,:) = look(lfp_fileroot,bad_channels(i),[],1)';
    end
    eegplot(bad_eeg(7,:),'srate',srate,'title','Bad Channels');
    good_channels = channels(~ismember(channels,bad_channels));
    good_eeg = [];
    for i=1:length(good_channels)
        good_eeg(i,:) = look(lfp_fileroot,good_channels(i),[],1)';
    end
    eegplot(good_eeg,'srate',srate,'title','Good Channels');
end