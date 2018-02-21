% Calculate tstats for power across all patients by anatomical label and
% make summary color plot

% User Inputs
subjList = {'HUP084','HUP087','HUP091','HUP092','HUP111'};
task     = 'motormap';

% set dirs
dirs = le_dirs;

anatList = {};

% Create anatList including all locations across patients
for i=1:length(subjList)
    subj = subjList{i};
    
    % select motor sites
    %[eLbl_list_motor,eLbl_list_nonmotor] = le_mm_selectMotorSites (subj)
    
    % load anatStruct
    [anatStruct] = le_centroids2Anatstruct(subj);
    
    % add additional anatomical areas to overall cell
    for j=1:length(anatStruct)
        if ~ismember({anatStruct(j).anatAbbr},anatList)
            anatList(length(anatList)+1)={anatStruct(j).anatAbbr};
        end
    end  
end
