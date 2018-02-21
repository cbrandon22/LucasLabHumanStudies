function [ anatStruct ] = le_loadAnatStruct( subj,eLbl_list )
%This function loads the anat_structure for each subject
% Inputs:
% subj.... 'HUP001'
% elecLbl_list .... optional; filters by a list of elecLbls


dirs = le_dirs;
cd(fullfile(dirs.data,'eeg',subj,'tal'))
anatStruct = subj_anatStruct; 

if exist('eLbl_list','var') && ~isempty(eLbl_list)
    idx_to_keep = false(1,length(anatStruct));
    for e = 1:length(eLbl_list)
        %if sum(~cellfun(@isempty,strcmp({anatStruct.eLbl},eLbl_list{e})))>1
        %    keyboard
        %end
        idx_to_keep(strcmp({anatStruct.eLbl},eLbl_list{e})) = true;
        
        %idx_to_keep(~cellfun(@isempty,strfind({anatStruct.eLbl},eLbl_list{e}))) = true;
    end
    
    %filter anatstruct and tStruct
    anatStruct = anatStruct(idx_to_keep);
end
    

