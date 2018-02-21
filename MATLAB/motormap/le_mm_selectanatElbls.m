function [eLbl_list] = le_mm_selectanatElbls(subj)
% Select bipolars within anatomical regions to plot, regardless of
% activity

%select regions
regions= {'L-HIP','R-HIP','L-PHIP','R-PHIP','L-AMYG','R-AMYG'};
%Initialize vars
config_mm = le_mm_config;

% set dirs
dirs = le_dirs;

% lGenerate anatStruct
[anatStruct] = le_centroids2Anatstruct(subj);
eLbl_list = [];
for i=1:length(anatStruct)
   if ismember({anatStruct(i).anatAbbr},regions)
       eLbl_list = [eLbl_list,{anatStruct(i).elecLbl}];
   end
end
end
