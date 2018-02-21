function [phaseStruct,config_phase] = le_calcPhaseReset(subj, sessList, elecNum,task,comparison,fRange,...
    tRange,collapseFreqFlag,collapseTimeFlag,phaseConfigNum);
   
%% [tStruct] = le_calcPhaseReset(subj, sessList, elecNum,task,comparison,phaseConfigNum)
%This function computes phase reset testing the non-uniformity of phase distributions across
%particular task conditions using circularANOVA.  Must specify a subjec, electrode number, task, and comparison
%(allows you to select between various comparisons of interest for a
%particular task). 
%Optional inputs include:
%sessList (particular list of sessions), 
%fRange (frequency of interest (default = []; loops through all frequencies)
%tRange (default all time bins in config_phase)
%phaseConfigNum (default =1 )

if ~exist('sessList','var') || isempty(sessList)
    sessList = [];
    sessLbl = 'all';
else
    sessLbl = [sessList{:}];
end
if ~exist('fRange','var') || isempty(fRange)
    fRange = [];
end
if ~exist('tRange','var') || isempty(tRange)
    tRange = [];
end
if ~exist('collapseFreqFlag','var') || isempty(collapseFreqFlag)
    collapseFreqFlag = true;
end
if ~exist('collapseTimeFlag','var') || isempty(collapseTimeFlag)
    collapseTimeFlag = true;
end
if ~exist('powConfigNum','var') || isempty(phaseConfigNum)
    phaseConfigNum = 1;
end

%% dirs and config
dirs = le_dirs;
config_task = le_config_task(task,comparison);

%% load all events
[events] = le_loadEvents(subj,sessList,task,config_task.evLabel);

%% load all power
[subjPhase,config_phase] = le_calcPhase(subj,sessList,elecNum,task,config_task.evLabel,phaseConfigNum);

%% error check: has there power for each event?
if length(events)~=size(subjPhase,3)
    error('There is a mismatch between the number of events loaded and the number of trials in the power matrix')
end

%% get retainIdx (logical index selecting the appropriate events to be used in the comparison)
[retain1,retain2] = le_events2retain(events,task,comparison);


%% compute t-stats
if collapseFreqFlag
    freqLbl = [num2str(fRange(1)) ':' num2str(fRange(2))];
else
    freqLbl = ['TF-' num2str(fRange(1)) ':' num2str(fRange(2))];
end
    
saveDir = fullfile(dirs.scratch,'phaseReset',freqLbl,...
    subj,config_task.task,comparison,sessLbl);

saveFile = [num2str(elecNum(1)) '-' num2str(elecNum(2)) '.mat'];
cd_mkdir(saveDir)

if exist([saveFile '.mat'],'file')
    load([saveFile '.mat'])
else
        % Identifying Information
        phaseStruct.subj = subj;
        phaseStruct.task = config_task.task;
        phaseStruct.comparison = comparison;
        phaseStruct.retLbl1 = config_task.retLbl1;
        phaseStruct.retLbl2 = config_task.retLbl2;
        phaseStruct.sessLbl = sessLbl;
        phaseStruct.eLbl = [num2str(elecNum(1)) '-' num2str(elecNum(2))];
    
        %compute tstats from retain
        [phaseStruct.zMat,phaseStruct.pMat,phaseStruct.fMat] = le_phaseReset_local(subjPhase,retain1,retain2,fRange,...
        tRange,collapseFreqFlag,collapseTimeFlag,config_phase);
end

cd(saveDir)
save(saveFile,'phaseStruct');

function[zMat,pMat,fMat] = le_phaseReset_local(phaseMat,retain1,retain2,fRange,tRange,collapseFreqFlag,collapseTimeFlag,config);
%computes phase reset using circ_ktest
%outputs

% A and B are the phase matrixes for each condition
A = phaseMat(:,:,retain1);
B = phaseMat(:,:,retain2);

% parse tRange, collapseTimeFlag
if isempty(tRange)
    tRange = config.timeBins;
else
    tbins = nanmean(config.timeBins,2)';
    [~,tInd_start] = min(abs(tRange(1) - tbins));
    [~,tInd_end] = min(abs(tRange(2) - tbins));
    tInd = tInd_start:tInd_end;
    tRange = config.timeBins(tInd,:);
end
% This flag collapses it into one time window
if collapseFreqFlag 
    fRange = [tRange(1,1) tRange(size(tRange,1),2)];
end

% parse fRange, collapseFreqFlag
% if fRange has multiple rows (suggesting freq bands, it wont edit it)
if isempty(fRange)
    fRange = [config.freQ' config.freQ'];
elseif size(fRange,1) == 1
    [~,fInd_start] = min(abs(fRange(1) - config.freQ));
    [~,fInd_end] = min(abs(fRange(2) - config.freQ));
    fInd = fInd_start:fInd_end;
    fRange = [config.freQ(fInd)' config.freQ(fInd)'];
end

% This flag collapses it into one frequency range
if collapseFreqFlag 
    fRange = [fRange(1,1) fRange(size(fRange,1),2)];
end

% initialize (loop through f and time bins) 
fMat = nan(size(fRange,1) ,size(tRange,1));
pMat = nan(size(fRange,1) ,size(tRange,1));
zMat = nan(size(fRange,1) ,size(tRange,1));

% loop through freq
for f = 1:size(fRange,1) 
    for t = 1:size(tRange,1) 
        %get freq range (forbinning)
        [~,fInd_start] = min(abs(fRange(f,1) - config.freQ));
        [~,fInd_end] = min(abs(fRange(f,2) - config.freQ));
        fInd = fInd_start:fInd_end;
        
        % calculate time epoch of interest
        tInd = [tRange(t,1):tRange(t,2)] - config.offsetMS;
  
        % collapse phase dist using circ_mean
        phase_dist_A = shiftdim(circ_mean(circ_mean(A(fInd,tInd,:),[],1),[],2),1)';
        phase_dist_B = shiftdim(circ_mean(circ_mean(B(fInd,tInd,:),[],1),[],2),1)';

        % do circ_ktest
        [pMat(f,t),fMat(f,t),zMat(f,t)]  = circ_ktest(phase_dist_A,phase_dist_B);
    end
    disp(f)
end




function [outStruct] = le_config_task(task,comparison)
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
        end
end
% 
% % calc rbar
% rbar1 = compute_rbar(phaseMat(:,:,retain1));
% rbar2 = compute_rbar(phaseMat(:,:,retain2));

% function [rbar] = compute_rbar(phase_dist)
% % This function will compute rbar or the degree of non-uniformity in a
% % phase distribution. This can be phase values from a single electrode 
% % or phase differences between two electrodes
% % The phase_dist matrix must be in the following format:
% % rows = F 
% % columns = T
% % 3rd dim = trial
% rbar = sqrt(((sum(cos(phase_dist),3)).^2) + ((sum(sin(phase_dist),3)).^2))/size(phase_dist,3);