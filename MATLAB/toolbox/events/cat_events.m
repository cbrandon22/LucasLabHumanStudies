%This script can be used to concatonate events from a particular subject
%and save under /data/events. It requires the user to enter the paths to
%each session. This allows to concatonate sessions in a particular folder even
%even when they are across different "laptop" directories. Also, cats math
%and recog events for pyFR
%Written by AGR 09-13-2013 (ashwinramayya@gmail.com)


% USER INPUTS
subj = 'HUP168';
task = 'CCDT'; % options: 'CCDT'; 'oddball'; 'netLearn';
sessList = {'Session_0'}; % 'Session_1','Session_2','Session_3'}; % Must be in order
dirs = le_dirs(task);
hadProblems = true; % true or false (no quotations)
description = 'no bad leads info'; % if there is a problem, describe it here
sessDirBase = fullfile(dirs.data,'eeg',subj,'behavioral');
sessDirList = {fullfile(sessDirBase,sessList{1})};
for i=2:length(sessList)
    sessDirList(i) = {fullfile(sessDirBase,sessList{i})};
end

evDir = dirs.events;
taskDir = evDir;

ev = [];

if strcmp(task,'pyFR') || strcmp(task,'RAM_FR1') || strcmp(task,'RAM_FR2')
    mathEv = [];
    recogEv = [];
end

for i = 1:length(sessDirList)
    % load events
    cd(sessDirList{i})
    load('events.mat')
    ev = [ev events];

    if strcmp(task,'RAM_FR1') || strcmp(task,'RAM_FR2')
        if exist('MATH_events.mat','file')
            load('MATH_events.mat')
            mathEv = [mathEv events];
        end
    end
    
    if strcmp(task,'pyFR')
        if exist('MATH_events.mat','file')
            load('MATH_events.mat')
            mathEv = [mathEv events];
        end
        if exist('RECOG_events.mat','file')
            load('RECOG_events.mat')
            recogEv = [recogEv events];
        end
    end
end

%save events
cd(taskDir)

%save expInfo
save([subj '_expinfo'],'hadProblems','description')


% check if events already exist
if exist([subj '_events.mat'],'file')
    s = input('WARNING: _events.mat file already exists, Overwrite??? (y/n)','s');
    if strcmp(s,'y')
        wrFlag = 1;
    elseif strcmp(s,'n')
        wrFlag = 0;
    end
else
    wrFlag = 1;
end

if wrFlag == 1
    % save  events;
    events = ev;
    if ~isempty(events)
    save([subj '_events'],'events');

    end

    if strcmp(task,'pyFR') || strcmp(task,'RAM_FR1') || strcmp(task,'RAM_FR2')
        %save math events
        events = mathEv;
        if ~isempty(events)
        save([subj '_math'],'events');
        end

        %save recog events
        events = recogEv;
        if ~isempty(events)
        save([subj '_recog'],'events');
        end
    end
end