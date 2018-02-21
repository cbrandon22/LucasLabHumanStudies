%This script will wrap around le_calcPow.m
subj       = 'HUP001';
sessList   = [];% (optional; if empty, it loops through all sessions seperately all events); 
                %Note: must be a cell array of strings
task       = 'motormap'; %the label of the task dir in /data/events
evLabel    = [];
configNum  = [];


% set dirs
dirs = le_dirs;

% load anatStruct (must go to subject's tal dir to get that anatStruct)
cd(fullfile(dirs.data,'eeg',subj,'tal'))
anatStruct = subj_anatStruct; 


% This is the main move - we are running le_calcPow
for i=1:length(anatStruct)
    le_calcPhase(subj,sessList,anatStruct(i).elecNum,task,evLabel,configNum);
    disp(i)
end