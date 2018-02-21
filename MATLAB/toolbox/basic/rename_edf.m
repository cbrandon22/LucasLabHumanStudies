% Batch rename edf files
% Changes subject names of all edf files within a folder
%
% Enter old subject label and new subject label
% Set path to desired directory
% Make sure subject label in dir function matches the oldSubj label

oldSubj = 'UP045';
newSubj = 'HUP091';
cd /Volumes/Lucas_ECoG_Data/data/eeg/HUP091/eeg.reref/OLD
files = dir('UP045*');
if isempty(files);
    disp('Error: No such files in set directory')
end

for i=1:length(files)
    [pathname,filename,extension] = fileparts(files(i).name);
    keepName = files(i).name(length(oldSubj)+1:end);
    newFileName = strcat(newSubj,keepName);
    movefile([filename extension], [newFileName]);
end