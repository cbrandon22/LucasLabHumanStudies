function [tStruct,config_pow] = le_calcTs(subj, sessList, elecNum,task,alignment,comparison,fRange,collapseFreqFlag,powConfigNum);
%% [tStruct] = le_calcTs(subj, sessList, elecNum,task,comparison,powConfigNum)
%This function computes t-statistics comparing power distributions across
%particular task conditions.  Must specify a subject, electrode number, task, and comparison
%(allows you to select between various comparisons of interest for a
%particular task). Optional inputs include sessList (particular list of sessions), 
% fRange (frequency of interest (default = [70 200]) and powConfigNum (determines how power is calculated)
% INPUTS

if ~exist('sessList','var') || isempty(sessList)
    sessList = [];
    sessLbl = 'all';
else
    sessLbl = [sessList{:}];
end
if ~exist('fRange','var') || isempty(fRange)
    fRange = [70 200];
end
if ~exist('collapseFreqFlag','var') || isempty(collapseFreqFlag)
    collapseFreqFlag = true;
end
if ~exist('powConfigNum','var') || isempty(powConfigNum)
    powConfigNum = 1;
end

%% dirs and config
dirs = le_dirs(task);
config_task = le_config_task(task,alignment,comparison);

%% load all events
[events] = le_loadEvents(subj,sessList,task,config_task.evLabel);

%% load all power
[subjPow,config_pow] = le_calcPow(subj,sessList,elecNum,task,config_task.evLabel,powConfigNum);

%% error check: has there power for each event?
if length(events)~=size(subjPow,3)
    error('There is a mismatch between the number of events loaded and the number of trials in the power matrix')
end

%% get retainIdx (logical index selecting the appropriate events to be used in the comparison)
[retain1,retain2] = le_events2retain(events,task,alignment,comparison);


%% compute t-stats
if collapseFreqFlag
    freqLbl = [num2str(fRange(1)) ':' num2str(fRange(2))];
else
    freqLbl = ['TF-' num2str(fRange(1)) ':' num2str(fRange(2))];
end
    
saveDir = fullfile(dirs.scratch,'tStats',freqLbl,...
    subj,config_task.task,comparison,sessLbl);
if size(elecNum,2)>1
    saveFile = [num2str(elecNum(1)) '-' num2str(elecNum(2)) '.mat'];
else
    saveFile = [num2str(elecNum) '.mat'];
end
cd_mkdir(saveDir)

if exist([saveFile '.mat'],'file')
    load([saveFile '.mat'])
else
        % Identifying Information
        tStruct.subj = subj;
        tStruct.task = config_task.task;
        tStruct.comparison = comparison;
        tStruct.retLbl1 = config_task.retLbl1;
        tStruct.retLbl2 = config_task.retLbl2;
        tStruct.sessLbl = sessLbl;
        if size(elecNum,2)>1
            tStruct.eLbl = [num2str(elecNum(1)) '-' num2str(elecNum(2))];
        else
            tStruct.eLbl = [num2str(elecNum)];
        end
    
        %compute tstats from retain
        [tStruct.tMat,tStruct.pMat,tStruct.pow1,tStruct.pow2,...
            tStruct.fInd, tStruct.tInd] = le_pow2tstat(subjPow,retain1,retain2,fRange,config_task.tRange,collapseFreqFlag,config_pow);
end

cd(saveDir)
save(saveFile,'tStruct');



function [outStruct] = le_config_task(task,alignment,comparison)
%this config maintains comparison-specific paramters.
switch task
    case 'motormap'
        outStruct.task = 'motormap';%'motormap';%'trackball';%'prob_sel2';%
        outStruct.evLabel = 'events';
        outStruct.tRange = [-1000 5000];
        switch comparison
            case 'moveWait'
                outStruct.retLbl1 = 'move';
                outStruct.retLbl2 = 'wait';
            case 'leftWait'
                outStruct.retLbl1 = 'left';
                outStruct.retLbl2 = 'wait';
            case 'rightWait'
                outStruct.retLbl1 = 'right';
                outStruct.retLbl2 = 'wait';
            case 'mouthWait'
                outStruct.retLbl1 = 'mouth';
                outStruct.retLbl2 = 'wait';
            case 'moveInstruct'
                outStruct.retLbl1 = 'move';
                outStruct.retLbl2 = 'instruct';
            case 'handWait'
                outStruct.retLbl1 = 'hand';
                outStruct.retLbl2 = 'wait';
            case 'instructWait'
                outStruct.retLbl1 = 'instruct';
                outStruct.retLbl2 = 'wait';
        end
    case 'CCDT'
        outStruct.task = 'CCDT';
        outStruct.evLabel = 'events';
        outStruct.tRange = [-2000 2500];
        switch alignment
            case 'cue'
                switch comparison
                    case 'delay'
                        outStruct.retLbl1 = 'Short';
                        %outStruct.retLbl2 = 'Mid';
                        outStruct.retLbl2 = 'Long';
                    case 'RT'
                        outStruct.retLbl1 = 'Fast';
                        outStruct.retLbl2 = 'Slow';
                    case 'time'
                        outStruct.retLbl1 = 'Early';
                        outStruct.retLbl2 = 'Late';
                end
            case 'CC'
                outStruct.retLbl1 = strcmp({events.type},'CC');
            case 'response'
                outStruct.retLbl1 = strcmp({events.type},'RESPONSE');
        end
 end
    

