%This wrapper function runs through the functions that are necessary to
%process raw behavioral and ecog data, and generate a report of the task.
%The steps are as follows:
%STAGE 1: Post-processing data

%1) Make behavioral events (task-specific, e.g. extractMOTORMAPevents.m, buildCCDTevents)

%2) process raw eeg file  (edf_anon and edf_split_wrapper)

%3) align the eeg data with the behavioral data, and update the events structure (alignTool)

%4) Run cat_events.m to concanate the Events


%STAGE 2: Electrophyiology Report Generation
%To Do:
%( ) write description of subj_anatStruct, le_calcPow_wrapper,  functions
%

%1) set dirs appropriately
%%%%% le_dirs.m

%2) set up subject's anatomical structure (list of bipolar pairs, w/ CT
%coords, MRI coords etc). Can be done automatic or manually.
%%%%% subj_anatStruct.m (data_structure in the subject's tal directory)
% - manually enter electrode pairs and labels

%3) for each bipolar pair, calculate ERPs and power for all events (config file is contained 
%within the power generation code)
%%%%% le_calcPow.m and le_calcPow_wrapper.m (to loop over a list)

%4) for each bipolar pair, for a given comparison, calculate t-stats for an
%freq of interest (default high gamma), and comparison of interest
%%%%% le_calcTs.m and le_calcTs_wrapper.m 

%5) Plot these data for each electrode
%%%%% le_plotTs.m

%6) Select motor sites for low frequency power and phase analysis. Open
%le_mm_readme.m for a walkthrough of the process.
%%%%% le_mm_selectMotorSites.m and le_mm_cat_tfPowPhase.m

%STAGE 3: Display results of t-test on the brain
%( ) use voxTool and freeSurfer to generate localizations







