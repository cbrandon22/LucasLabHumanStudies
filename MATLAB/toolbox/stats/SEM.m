%This function calculates SEM for a group of vectors
function[SEM_groupvar]=SEM(groupvar)

SEM_groupvar=nanstd(groupvar,[],1)/sqrt(size(groupvar,1));

