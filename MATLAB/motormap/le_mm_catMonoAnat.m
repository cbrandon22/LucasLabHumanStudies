subjList = {'HUP084','HUP087','HUP091','HUP092','HUP111'};%'HUP089','HUP090'
dirs = le_dirs('motormap');
config_mm = le_mm_config;
allMonoAnat = {};
motorRegions = {'R-PRE','R-Post','L-PRE','L-POST'};
for i=1:length(subjList)
    subj = subjList{i};
    thisSubjAnat = {};
    [~,mono_anatStruct] = le_centroids2Anatstruct(subj);
    for j=1:length(mono_anatStruct)
        if ismember({mono_anatStruct(j).anatAbbr},motorRegions)
            label = 1;
        else
            label=2;
        end
        thisSubjAnat(j,1)={subj};
        thisSubjAnat(j,2)={mono_anatStruct(j).X};
        thisSubjAnat(j,3)={mono_anatStruct(j).Y};
        thisSubjAnat(j,4)={mono_anatStruct(j).Z};
        thisSubjAnat(j,5)={label};
    end
    allMonoAnat=[allMonoAnat;thisSubjAnat];
end