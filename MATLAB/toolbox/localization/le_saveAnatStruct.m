cd /Volumes/Lucas_ECoG_Data/scratch/anatStructs

subjList = {'HUP084','HUP087','HUP089','HUP090','HUP091','HUP092','HUP111'};
for i=1:length(subjList)
    subj = subjList{i};
    [anatStruct_bipol,anatStruct_mono] = le_centroids2Anatstruct(subj)
    bpfilename = strcat(subj,'_anatStruct_bipol.mat');
    monofilename = strcat(subj,'_anatStruct_mono.mat');
    save(bpfilename,'anatStruct_bipol');
    save(monofilename,'anatStruct_mono');
end