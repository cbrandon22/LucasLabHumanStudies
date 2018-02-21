dirs = le_dirs;
cd (dirs.scratch);
load('moveWaitClusterStruct.mat');
e1 = cell(length(tStruct_lf),1);
e2 = cell(length(tStruct_lf),1);
elec1 = cell(length(tStruct_lf),1);
elec2 = cell(length(tStruct_lf),1);
for i=1:length(tStruct_lf)
    eLbl = strsplit(tStruct_lf(i).eLbl,'-');
    e1(i) = eLbl(1);
    e2(i) = eLbl(2);
    if ~exist('subj','var') | ~strcmp(subj,{tStruct_lf(i).subj})
        subj = tStruct_lf(i).subj;
        [~,monoStruct] = le_centroids2Anatstruct(subj);
        elecNums = cell(length(monoStruct),1);
        for ii=1:length(monoStruct)
            elecNums{ii} = num2str(monoStruct(ii).elecNum);
        end
    end
    [~,locb1] = ismember(e1(i),elecNums);
    [~,locb2] = ismember(e2(i),elecNums);
    elec1{i} = monoStruct(locb1).elecLbl;
    elec2{i} = monoStruct(locb2).elecLbl;
end
fclose('all');