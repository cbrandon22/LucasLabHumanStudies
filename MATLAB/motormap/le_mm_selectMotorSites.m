function [eLbl_list_motor,eLbl_list_nonmotor] = le_mm_selectMotorSites (subj)
% This function returns a list of motor and non-motor sites by performing
% FDR-correction of distribution of p-values (move to cluster based
% correction in the future?)

% Inputs: 
% subj  =                'HUP001'


% 1) 'tStruct_hfa' for HFA band (70-200) for 'moveWait' comparison for all electrodes. This is
% created by the report generator pipeline
% ('le_calcPow' and wrapper,'le_calcTs'(uses pow2stat) and wrapper).

% 2) 'config_pow' associated with this tStruct

% Load both of these variables using:
% 
% Outputs:
% elecLbl_list_motor     elecLbls of motor site electrodes
% elecLbl_list_nonmotor  elecLbls of non motor site electrodes
% AGR 05/2016

%Initialize vars
config_mm = le_mm_config;

% set dirs
dirs = le_dirs('motormap');

% lGenerate anatStruct
[anatStruct] = le_centroids2Anatstruct(subj);

% load tStruct and configPow
[tStruct_hfa,config_pow] = le_calcTs_wrapper(subj,[],'motormap',[],config_mm.comparison,...
    config_mm.hfa_range,config_mm.hfa_collapseFreqFlag,config_mm.powConfigNum);

% concatonate p values from each elec to a single matrix
if size(tStruct_hfa(1).pMat,1)>1
    error ('tStruct needs to be computed for a single HFA band...pMat needs to be a row vec for each elec')
end
pMat = cat(1,tStruct_hfa.pMat);

% % first, identify post-movement time bins
% tIdx_postMove = find(mean(config_pow.timeBins,2) > 0,1);
% 
% % clean pMat by removing prestim values
% pMat = pMat(:,tIdx_postMove:end);
% apply FDR correction on all p-values
% pMat_sig_idx = pMat <= 0.001;
pID = fdr(pMat,config_mm.fdr_p_thresh);
if isempty(pID)
    eLbl_list_motor = {};
    eLbl_list_nonmotor = {anatStruct.eLbl};
    return;
end
pMat_sig_idx = pMat <= max(pID);   %0.001;
% identify motor sites as those that have a sig p - value during the
% post-movement interval

% first, identify post-movement time bins
tIdx_postMove = find(mean(config_pow.timeBins,2) > 0,1);

% clean pMat by removing prestim values
pMat = pMat(:,tIdx_postMove:end);
pMat_sig_idx = pMat_sig_idx(:,tIdx_postMove:end);

% identify motor electrodes
motor_idx = sum(pMat_sig_idx,2)>config_mm.num_sig_bins_thresh;


% anat
eLbl_list_motor = {anatStruct(motor_idx).eLbl};
eLbl_list_nonmotor = {anatStruct(~motor_idx).eLbl};