function allRTs = CCDT_loadRTs(subj,sessNum)
%%% Loads reaction times for given sessions, default = all sessions.
%%% Returns allRTs (RT, response time, fixDur, sessNum, setNum)
dirs = le_dirs;
dataDir = fullfile(dirs.data,'eeg','CCDT',subj,'behavioral');

allRTs = [];
if isempty(sessNum)
    sessNum = 0;
    while exist(fullfile(dataDir,['Session_' num2str(sessNum)]),'dir')==7
        sessDir = (fullfile(dataDir,['Session_' num2str(sessNum)]));
        cd(sessDir);
        setNum = 1;
        while exist(fullfile(sessDir,['set' num2str(setNum)]),'dir')==7
            cd(fullfile(sessDir,['set' num2str(setNum)]));
            load('sessRTs.mat');
            sessRTs(:,4) = sessNum;
            sessRTs(:,5) = setNum;
            allRTs = [allRTs; sessRTs];
            setNum = setNum + 1;
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
        setNum = 1;
        while exist(fullfile(sessDir,['set' num2str(setNum)]),'dir')==7
            cd(fullfile(sessDir,['set' num2str(setNum)]));
            load('sessRTs.mat');
            sessRTs(:,4) = sessNum(i);
            sessRTs(:,5) = setNum;
            allRTs = [allRTs; sessRTs];
            setNum = setNum + 1;
        end
    end
end
end