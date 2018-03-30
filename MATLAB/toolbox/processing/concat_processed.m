% General non-task related processing to combine split sessions, must run
% edf_split_wrapper on both parts to create eeg.noreref folder of channel
% folders first.
ddir = '/Volumes/HumanStudies/HumanStudies/sleep/eeg';
p1fn = fullfile(ddir,'HUP166_1.1/eeg.noreref/HUP166_1.1_20Mar18_2320');
p2fn = fullfile(ddir,'HUP166_1.2/eeg.noreref/HUP166_1.2_21Mar18_0009');
outdir = fullfile(ddir,'HUP166_1/processed');

% Read jacksheet
fid = fopen(fullfile(ddir, 'HUP166_1.1/eeg.noreref/jacksheet.txt'));
C = textscan(fid,'%s');
fclose(fid);
JAC{:,1} = C{1}(1:2:end-1);
JAC{:,2} = C{1}(2:2:end);
cNum = cellfun(@str2double, JAC{1});
clbl = JAC{2};
leadNames = regexp(clbl,'\D*','match');
leadNames = vertcat(leadNames{:});

elecInfo = {};
for i=1:length(cNum)
    elecInfo = [elecInfo;{cNum(i)},clbl{i}];
end

%Read parameters
fid = fopen(fullfile(ddir, 'HUP166_1.1/eeg.noreref/params.txt'));
C = textscan(fid,'%s');
fclose(fid);
srate= str2double(C{1}(2,1));
dat_gain = str2double(C{1}(6,1));

%Save session info
if ~exist(outdir,'dir'),mkdir(outdir);end
save(fullfile(outdir, 'sessInfo.mat'),'srate','elecInfo','dat_gain');

i=1;
while i<length(cNum)
    thisLead = leadNames(i);
    disp(thisLead{1});
    lead_info = [cNum(strcmp(leadNames,thisLead)),{clbl(strcmp(leadNames,thisLead))}];
    eeg = [];
    while strcmp(leadNames(i),thisLead)
        chanfile1 = sprintf('%s.%03i', p1fn,cNum(i));
        chanfile2 = sprintf('%s.%03i', p2fn,cNum(i));
        ceeg=[look(p1fn,cNum(i),[],1)',look(p2fn,cNum(i),[],1)'];
        eeg = [eeg;ceeg];
        i = i+1;
    end
    save(fullfile(outdir, [thisLead{1} '.mat']),'eeg','lead_info','-v7.3');
end