% General non-task related processing to combine split sessions, must run
% edf_split_wrapper on both parts to create eeg.noreref folder of channel
% folders first.
ddir = '/Volumes/HumanStudies/HumanStudies/sleep/eeg';
concat = 0; %combine 2 different eeg.noreref (useful for long recordings (~> 1hr:20min))
p1fn = fullfile(ddir,'HUP167_a2/eeg.noreref/HUP167_a2_07Apr18_1702');
%p2fn = fullfile(ddir,'HUP166_1.2/eeg.noreref/HUP166_1.2_21Mar18_0009');
outdir = fullfile(ddir,'HUP167_a2/processed');

% Read jacksheet
fid = fopen(fullfile(ddir, 'HUP167_a2/eeg.noreref/jacksheet.txt'));
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
fid = fopen(fullfile(ddir, 'HUP167_a2/eeg.noreref/params.txt'));
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
        if concat
            chanfile2 = sprintf('%s.%03i', p2fn,cNum(i));
            ceeg=[look(p1fn,cNum(i),[],1)',look(p2fn,cNum(i),[],1)'];
        else
            ceeg=look(p1fn,cNum(i),[],1)';
        end
        eeg = [eeg;ceeg];
        i = i+1;
    end
    save(fullfile(outdir, [thisLead{1} '.mat']),'eeg','lead_info','-v7.3');
end