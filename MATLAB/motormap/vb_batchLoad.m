function [batchFile]=vb_batchLoad(subj, loadfile, renamefile)
%% Load and store the same subject specific file into one structure to enable batch analysis

% Inputs
%        subj: cell array of subject specific names preceding loadfile
%        loadfile: string defining subject specific file of interest for batch analysis
%        renamefile: cell array of variables of interest for batch analysis
%    all subj specific files of interest should be in current directory


%    examples
%             subj=allSubj.name; 
%             loadfile = '_X.mat'
%             renamefile= {'var1' 'var2'}; 


% Output
%        batchFile: batchFile.name
%                   batchFile.renamefiles

% example command line: [plvBatch]=vb_batchLoad(allSubj.name, '_PLV.mat', {'plvSourceSpace'});


% VPB - vivek.buch@uphs.upenn.edu
% University of Pennsylvania
% Dept of Neurosurgery, PGY3

%%
clear i j xy
xy=0;
for i=1:length(subj)
    if exist([subj{i} loadfile])
        disp(subj{i})
        xy=xy+1;
        batchFile{xy}.name=subj{i}; 
        for j=1:length(renamefile)
            load([subj{i} loadfile], renamefile{j});
            newfile=([subj{i} '_' renamefile{j}]);
            disp(newfile)
            batchFile{xy}.(renamefile{j})=eval(renamefile{j}); 
            clear(renamefile{j})
        end
    else 
        disp([subj{i} ' does not have '  loadfile ' file, skipping...'])
        continue;
    end
end
    