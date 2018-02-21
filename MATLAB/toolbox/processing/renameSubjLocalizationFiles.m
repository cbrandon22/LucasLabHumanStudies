function changeLog = renameSubjLocalizationFiles%(baseDir,inSubj,outSubj)
% This script renames all files within baseDir and all sub-directories,
% replacing all instances of inSubj with outSubj.
% Intended to change RID### labels to HUP### labels from imaging data.
% Written 2/6/2018 - CB, cbrandon22@gmail.com

%Inputs:
baseDir = '/Volumes/HumanStudies/HumanStudies/localization/HUP154/post';
inSubj = 'RID146'; %old subject label
outSubj = 'HUP154';%new subject label

changeLog = [];% log changes just in case...
changeLog = renameSubjsInDir(baseDir,inSubj,outSubj,changeLog);


function changeLog = renameSubjsInDir(baseDir,inSubj,outSubj,changeLog)
% Actually change all file/folder names given directory
files = dir(baseDir);
for i=1:length(files)
    if ~ismember(files(i).name,{'.','..'}) % skip pwd and upper directory
        if files(i).isdir % If this is a folder, enter and restart
            newDir = fullfile(baseDir,files(i).name);
            tmp_changeLog = renameSubjsInDir(newDir,inSubj,outSubj,changeLog); % recursive call
            changeLog = tmp_changeLog;
        end
        start_idx = strfind(files(i).name,inSubj);
        if ~isempty(start_idx)
            infn = fullfile(baseDir,files(i).name);
            outfn = fullfile(baseDir,[files(i).name(1:start_idx-1) outSubj files(i).name(start_idx+length(inSubj):end)]);
            if ~isempty(files(i).date) % Alias files cause problems. May need to adjust if non-alias files have empty date strings
                changeLog = [changeLog; {infn},{outfn}];
                movefile(infn,outfn);
            else
                disp(['Skipping: ' files(i).name]);
            end
        end
    end
end