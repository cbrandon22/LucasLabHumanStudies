%% Batch Run Function

% input files for batch function
dir = '/media/vpbuch/VIVEK/motormap/';
load('/media/vpbuch/VIVEK/motormap/allSubj.mat')
subj=allSubj.name(1:11);
resampleFreq=512;

clear i
for i=1:length(subj)
    load([dir subj{i} '_eegdatnew.mat'])
    load([dir 'anat_data/' subj{i} '_talLocs_database_monopol.mat'])
    disp(subj{i})
    
    % RUN function of interest
    try
        [plvAll, plvSourceSpace, plvBrodmannSpace, t] = vbrun_PLV(dat, ind, talStruct, resampleFreq, subj{i});
        disp('done')
    catch
        EX
    end
    clear dat ind talStruct plv* t
end