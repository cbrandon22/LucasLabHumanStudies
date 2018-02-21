task = 'CCDT';
subjList = {'HUP140'};
dirs = le_dirs(task);
cd(dirs.events)
for sub=1:length(subjList)
    load([subjList{sub} '_events']);
    for i=1:length(events)
        if ~isempty(events(i).eegfile)
            thisEEGfile = events(i).eegfile;
            if strcmp(task,'CCDT')
                processedSubDir = strfind(thisEEGfile,['/' subjList{sub} '/']);
                events(i).eegfile = fullfile(dirs.data,'eeg',thisEEGfile(processedSubDir(1):end));
            else
                processedSubDir = strfind(thisEEGfile,'eeg/');
                events(i).eegfile = fullfile(dirs.data,thisEEGfile(processedSubDir(1):end));
            end
        end
    end
    save([subjList{sub} '_events.mat']);
end
disp('done')