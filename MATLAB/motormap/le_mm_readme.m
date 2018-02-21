% This script is a readme for motor map analyses. The general goal of these
% analyses is to perform hase-based connectivity analyses based on
% pre-selected motor a HFA sites so as to understand connectivity patterns
% associated with "motor sites".There are two parts of these
% analyses: 1) Scripts to perform HFA-based power analyses on all
% electrodes and return a list of motor and non-motor sites, 2) perform
% pairwise low-freq power and phase-based analyses on these electrode pairs (motor - motor;
% non-motor - non-motor; motor - non-motor)
% AGR 05/2016


% Config file for these analyses
% config_mm = le_mm_config;

% For each subject, 
%subj = 'HUP111';

% (X)  Select Motor Sites 

%  [eLbl_list_motor,eLbl_list_nonmotor] = le_mm_selectMotorSites (subj)

    % This function returns a list of motor and non-motor sites by performing
    % FDR-correction of distribution of p-values (move to cluster based
    % correction in the future?)
    
    % Inputs: subj
    % Outputs: anatStruct for motor and non motor sites
    
    % Comments:
    % it loads the following variables using: 
    % [tStruct_hfa,config_pow_hfa] = le_calcTs_wrapper('HUP001',[],'motormap','moveWait',[70 200],1);

    % 1) 'tStruct_hfa' for HFA band (70-200) for 'moveWait' comparison for all electrodes. This is
    % created by the report generator pipeline
    % ('le_calcPow' and wrapper,'le_calcTs'(uses pow2stat) and wrapper).

    % 2) 'config_pow_hfa' associated with this tStruct

    % Have these variables calculated using parallel computing before
    % running this function
    
    % It then selects significant motor sites using FDR correction


% (X) Compute low frequency time freq plots for motor electrodes and
% categorize the various movement-related spectral patterns

% le_mm_cat_tfPowPhase(subj,eLbl_list_motor,'motor')

     
    % ( ) Assess for TF phase reset (using powConfigNum = 2). Use PLV, perform single -site phase
    % distribution analyses. Should see a desynchronization. 
    
    % ( ) Assess for TF phase reset (using powConfigNum = 2). 
    % Perform phase locking through all non-motor sites. Hyp - should see
    % desynchronization with all brain regions
    
    % ( ) Assess for PAC - motor site HFA amp vs TF phase at non-motor
    % sites

