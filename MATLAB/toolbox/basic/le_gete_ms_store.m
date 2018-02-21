%This script will wrap around le_calcPow.m
subj       = 'HUP111';
sessList   = [];% (optional; if empty, it loops through all sessions seperately all events); 
                %Note: must be a cell array of strings
task       = 'motormap'; %the label of the task dir in /data/events
evLabel    = [];
configNum  = [];


% set dirs
dirs = le_dirs;

% load anatStruct (must go to subject's tal dir to get that anatStruct)
cd(fullfile(dirs.data,'eeg',task,subj,'tal'))
[anatStruct] = le_centroids2Anatstruct(subj);


% This is the main move - we are running le_calcPow
for i=1:length(anatStruct)
    elecNum = anatStruct(i).elecNum;


%% set dirs and parse inputs
    dirs = le_dirs;
    if ~exist('evLabel','var') || isempty(evLabel)
        evLabel = 'events';
    end
    if ~exist('configNum','var') || isempty(configNum)
        configNum = 1;
    end

    %% load config
    config = le_config_calcPow(configNum,task);
    
    
    %load events
    [events] = ec_loadEvents(subj,sessList,task,evLabel);
    unique([events.session])
    %% if ~exist sessList, create it here
    if ~exist('sessList','var') || isempty(sessList)
        sessList = unique({events.session});
    end

    %% loop through sessList
    for i = 1:length(sessList)
        sess = sessList{i};
        %cd_mkdir(fullfile(e_msDir,sess))
        %e_msFile = [num2str(elecNum(1)) '-' num2str(elecNum(2)) '.mat'];

            %filter events
            sessEv = events(strcmp({events.session},sess));
          
    end
    %this function calculates power for a particular session
    %gete_ms_wrapper
    [EEG,resampleFreq] = gete_ms_wrapper(elecNum,events, ...
                   config.durationMS,config.offsetMS,config.bufferMS,...
                   config.filtfreq, config.filttype, config.filtorder);
    %disp(resampleFreq)
    %e_msDir = fullfile(dirs.scratch,'e_ms',num2str(config.configID),subj,[task '_' evLabel]);
    %cd_mkdir(e_msDir);
end

