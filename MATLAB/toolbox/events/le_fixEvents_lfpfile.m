task = 'oddball';
subjList = {'HUP144_e'};
dirs = le_dirs(task);
for sub=1:length(subjList)
    cd(fullfile(dirs.data,'eeg',subjList{sub},'behavioral','session_0'))
    load('events.mat','events');
    for i=1:length(events)
        if ~isempty(events(i).lfpfile)
            thisEEGfile = events(i).lfpfile;
            if strcmp(task,'CCDT')
                processedSubDir = strfind(thisEEGfile,['/' subjList{sub} '/']);
                events(i).eegfile = fullfile(dirs.data,'eeg',thisEEGfile(processedSubDir(1):end));
            else
                processedSubDir = strfind(thisEEGfile,'eeg/');
                events(i).lfpfile = fullfile(dirs.data,thisEEGfile(processedSubDir(1):end));
            end
        end
    end
    save([subjList{sub} '_events.mat']);
end
disp('done')