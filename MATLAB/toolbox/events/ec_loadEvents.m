function[events] = ec_loadEvents(subj,sessList,task,evLabel)
%% [events] = ec_loadEvents(subj,sess,task,evLabel);
%This function loads events for a particular subject from a particular
%task. 
%Inputs 
%subj  = 'UP034'; 
%sessList = {'0'}; % (optional; if empty, it loads all events); Note: must be a
               % cell array of strings
%task = 'pyFR'; %the label of the task dir in /data/events
%evLabel = 'math' (optional: default is 'events', but allows you to load
          %math or recog from pyFR if needed)
          
%Outputs:
%events

%Written by AGR 09-23-2013

%% set dirs and parse inputs
dirs = le_dirs(task);
if ~exist('evLabel','var') || isempty(evLabel)
    evLabel = 'events';
end

%% go to ev dir and load events
if exist(fullfile(dirs.events),'dir')
    cd(fullfile(dirs.events))
    evFile = [subj '_' evLabel '.mat'];
    load(evFile)
else
    ev = [];
    sessDir = fullfile('/','data','eeg',subj,'behavioral',task);
    d = dir(['/data/eeg/' subj '/behavioral/' task '/session*']);
    for foo=1:length(d)
        cd(fullfile(sessDir,d(foo).name));
        if exist('events.mat','file')
            load('events.mat')
            ev = [ev events];
        end
    end
    events = ev;
end

%% conv session field to char
if ~ischar(events(1).session)
     for i = 1:length(events); events(i).session = num2str(events(i).session); end
end
%% filter by session if needed
if exist('sessList','var') && ~isempty(sessList)

    sess = {events.session};

    % identify which events to return
    retain = false(1,length(events));
    for i = 1:length(sessList)
       retain = retain|strcmp(sess,sessList{i}); 
    end
else
    retain = true(1,length(events));
end

%% filter events
events = events(retain);
