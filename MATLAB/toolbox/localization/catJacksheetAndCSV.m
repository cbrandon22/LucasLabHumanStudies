function catJacksheetAndCSV
task = 'CCDT';
subj = 'HUP143';
d = le_dirs(task);
saveFile=1;
saveFName = fullfile(d.data,'eeg',subj,'docs','jacksheetAnat.mat');
csvfile = 'electrodenames_native.csv';
csvfile = fullfile('/Volumes/HumanStudies/HumanStudies/localization',subj,'post',csvfile);
csvLabelCol = 1;
csvColumns2add = [2];

% open jacksheet
%jacFile = fullfile(d.data,'eeg',task,subj,'docs','jacksheet.txt');
jacFile = fullfile(d.data,'eeg',subj,'docs','jacksheet.txt');
fid = fopen(jacFile,'r');
JAC=textscan(fid,'%d%s','delimiter','\t');
JAC{:,2} = strtrim(JAC{:,2});
fclose(fid);
jacLabels = JAC{2};

fid=fopen(csvfile);
csvMat=textscan(fid,'%s%s','delimiter',',');
fclose(fid);
csvMat = [csvMat{1},csvMat{2}];

for c=1:size(jacLabels,1)
    [lia,locb]=ismember(jacLabels{c,1},{csvMat{:,csvLabelCol}});
    if ~lia,warning(sprintf('%s not found in csv', jacLabels{c,1}));end
    jacLabels{c,2:1+length(csvColumns2add)}=csvMat{c,csvColumns2add};
end
if saveFile
    save(saveFName,'jacLabels');
end