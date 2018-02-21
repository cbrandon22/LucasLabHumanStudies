%This script will wrap around le_calcPow.m
subjList = {'HUP149_e'};
for jj=1:length(subjList)
    subj       = subjList{jj};
    sessList   = {};% (optional; if empty, it loops through all sessions seperately all events);
    %Note: must be a cell array of strings
    task       = 'oddball'; %the label of the task dir in /data/events
    evLabel    = 'events';
    configNum  = 1;
    excludeChan = {};%{'DC1'};
    
    
    % set dirs
    dirs = le_dirs(task);
    
    % load anatStruct (must go to subject's tal dir to get that anatStruct)
    % cd(fullfile(dirs.data,'eeg',task,subj,'tal'))
    % [~,anatStruct] = le_centroids2Anatstruct(subj);
    
    % open jacksheet
    if strcmp(task,'oddball')
        jacFile = fullfile(dirs.data,'eeg',subj,'docs','electrodeMap.xlsx');
        [~,~,xlcells] = xlsread(jacFile);
        emptyChannels = [];
        lbls = {};
        for i=1:length(xlcells)
            if isnan(xlcells{i,3})
                emptyChannels = [emptyChannels, i];
            elseif i<129
                lbls = [lbls;{i}, xlcells{i,3}];
            end
        end
        JAC{1} = cell2mat(lbls(:,1));
        JAC{2} = lbls(:,2);
    else
        jacFile = fullfile(dirs.data,'eeg',subj,'docs','jacksheet.txt');
        fid = fopen(jacFile,'r');
        JAC=textscan(fid,'%d%s','delimiter','\t');
        JAC{:,2} = strtrim(JAC{:,2});
        fclose(fid);
    end
   
    [lia,locb] = ismember(excludeChan,JAC{2});
    if sum(lia)>0
        JAC{1}(locb(locb~=0)) = [];
        JAC{2}(locb(locb~=0)) = [];
    end
    
    config = le_config_calcPow(configNum,task);
    
    % This is the main move - we are running le_calcPow
    for i=1:length(JAC{:,1})
        disp([subj ': ' num2str(i) '/' num2str(length(JAC{:,1}))])
        le_calcPow(subj,sessList,JAC{1}(i),task,evLabel,configNum);
    end
    powConfigFile = fullfile(dirs.scratch,'POWER',num2str(config.configID),subj,[task '_' evLabel],'powConfig.mat');
    save(powConfigFile,'config')
end