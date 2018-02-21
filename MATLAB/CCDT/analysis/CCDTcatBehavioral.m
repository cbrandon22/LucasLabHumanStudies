function [allRTs, allSyncs] = CCDTcatBehavioral(subj,sessNum)
% Concatenates behavioral data for specified subject and sessions 
% (default = all sessions)
dirs = le_dirs('CCDT');
dataDir = fullfile(dirs.data,'eeg',subj,'behavioral');

allRTs = [];
allSyncs = [];
if isempty(sessNum)
    sessNum = 0;
    while exist(fullfile(dataDir,['Session_' num2str(sessNum)]),'dir')==7
        sessDir = (fullfile(dataDir,['Session_' num2str(sessNum)]));
        cd(sessDir);
        if exist('sessRTs.mat','file') == 2 % Not broken into sets, add session
            load('sessRTs.mat');
            load('syncTimes.mat');
            sessRTs(:,4) = sessNum;
            allRTs = [allRTs; sessRTs];
            allSyncs = [allSyncs; syncTimes];
        else % Concatenate sets within session
            setNum = 1;
            while exist(fullfile(sessDir,['set' num2str(setNum)]),'dir')==7
                cd(fullfile(sessDir,['set' num2str(setNum)]));
                load('sessRTs.mat');
                load('syncTimes.mat');
                sessRTs(:,4) = sessNum;
                allRTs = [allRTs; sessRTs];
                allSyncs = [allSyncs; syncTimes];
                setNum = setNum + 1;
            end
        end
        sessNum = sessNum + 1;
    end
else
    if iscell(sessNum) % Convert to vector of doubles
        sessNumStr = '';
        for i=1:length(sessNum)
            sess = strsplit(sessNum{i},'_');
            sessNumStr = [sessNumStr ' ' sess{2}];
        end
        sessNum = sscanf(sessNumStr,'%d');
    end
    for i=1:length(sessNum)
        sessDir = fullfile(dataDir,['Session_' num2str(sessNum(i))]);
        cd(sessDir);
        if exist('sessRTs.mat','file') == 2 % Not broken into sets, add session
            load('sessRTs.mat');
            load('syncTimes.mat');
            sessRTs(:,4) = sessNum;
            allRTs = [allRTs; sessRTs];
            allSyncs = [allSyncs; syncTimes];
        else % Concatenate sets within session
            setNum = 1;
            while exist(fullfile(sessDir,['set' num2str(setNum)]),'dir')==7
                cd(fullfile(sessDir,['set' num2str(setNum)]));
                load('sessRTs.mat');
                load('syncTimes.mat');
                sessRTs(:,4) = sessNum;
                allRTs = [allRTs; sessRTs];
                allSyncs = [allSyncs; syncTimes];
                setNum = setNum + 1;
            end
        end
    end
end
end