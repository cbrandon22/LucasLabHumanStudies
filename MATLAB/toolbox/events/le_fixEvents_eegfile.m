task = 'CCDT';
subjList = {'HUP069','HUP133','HUP136','HUP139','HUP140',...
    'HUP142','HUP143','HUP145','HUP146','HUP150','HUP152','HUP153',...
    'HUP154','HUP157'};
dirs = le_dirs(task);
for sub=1:length(subjList)
    load(fullfile(dirs.events,[subjList{sub} '_events']));
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