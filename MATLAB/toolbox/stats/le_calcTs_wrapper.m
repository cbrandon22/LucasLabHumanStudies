function[tStruct,config_pow] = le_calcTs_wrapper(subj,sessList,task,alignment,comparison,fRange,collapseFreqFlag,powConfigNum,eLbl_list)

%This wrapper calls le_calcTs.m to return tStruct and power configuration.
%Inputs
%subj           = 'HUP001';
%sessList       = [];% (optional; if empty, it loops through all sessions seperately all events); 
                    %Note: must be a cell array of strings
%task           = 'motormap'; %the label of the task dir in /data/events
%comparison     = 'moveWait';
%fRange         = [70 200];
%powConfigNum   = [];
%elecLblList    = cell arry of strings listing of elecLbls to plot, if empty, loads all tStructs (optional)
%collapseFreqFlag = if TRUE, it will collapse across the frequency range
                    %(default = TRUE)


% set dirs
dirs = le_dirs(task);

% load anatStruct (must go to subject's tal dir to get that anatStruct)
[anatStruct] = le_centroids2Anatstruct(subj);
tStruct = struct();

% filter for elecLbl_list (elecs of interest)
if exist('eLbl_list','var') && ~isempty(eLbl_list)
    idx_to_keep = false(1,length(anatStruct));
    for e = 1:length(eLbl_list)
        idx_to_keep(strcmp({anatStruct.eLbl},eLbl_list{e})) = true;     
    end
    
    %filter anatstruct and tStruct
    anatStruct = anatStruct(idx_to_keep);
end


% This is the main move - we are running le_calcTs
disp(['Calculating Session t-stats... ' subj])
for i=1:length(anatStruct)
    if i == 1
        [tStruct,config_pow] = le_calcTs(subj, sessList,anatStruct(i).elecNum,task,alignment,comparison,fRange,collapseFreqFlag,powConfigNum);
    else
        tStruct(i) = le_calcTs(subj, sessList,anatStruct(i).elecNum,task,alignment,comparison,fRange,collapseFreqFlag,powConfigNum);
    end
end


