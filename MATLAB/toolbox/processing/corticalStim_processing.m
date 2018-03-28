ddir = '/Volumes/HumanStudies/HumanStudies/corticalStim/eeg';
p1fn = fullfile(ddir,'HUP128_p1/eeg.noreref/HUP128_p1_21Dec16_1528');
p2fn = fullfile(ddir,'HUP128_p2/eeg.noreref/HUP128_p2_21Dec16_1630');

fid = fopen(fullfile(ddir, 'HUP128_p1/docs/jacksheet.txt'));
C = textscan(fid,'%s');
fclose(fid);
JAC{:,1} = C{1}(1:2:end-1);
JAC{:,2} = C{1}(2:2:end);

cNum = cellfun(@str2double, JAC{1});
clbl = JAC{2};
leadNames = regexp(clbl,'\D*','match');
leadNames = vertcat(leadNames{:});

i=1;
while i<length(cNum)
    thisLead = leadNames(i);
    disp(thisLead{1});
    eeg = [];
    while strcmp(leadNames(i),thisLead)
        chanfile1 = sprintf('%s.%03i', p1fn,cNum(i));
        chanfile2 = sprintf('%s.%03i', p2fn,cNum(i));
        ceeg=[look(p1fn,cNum(i),[],1)',look(p2fn,cNum(i),[],1)'];
        eeg = [eeg;ceeg];
        i = i+1;
    end
    save(fullfile(ddir,'HUP128/processed', [thisLead{1} '_dat.mat']),'eeg','JAC','-v7.3');
end