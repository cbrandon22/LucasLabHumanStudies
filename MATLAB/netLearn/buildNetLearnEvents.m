function buildNetLearnEvents
% builds events.mat and eeg.eeglog.up for CCDT session. Stores in session folder.
% events structure with fields:
% subj,session,trial#,eventType,item,msTime,msOffset
subj = 'HUP145';
session = 'Session_0';
dirs = le_dirs('netLearn');

sessDir = fullfile(dirs.data,'eeg',subj,'behavioral',session);
evFile = fullfile(sessDir,'events.mat');
eegLogFile = fullfile(sessDir,'eeg.eeglog.up');

% print opening line
fprintf('  Making netLearn EVENTS for ')
fprintf('%s, %s: \n',subj,session)


%--------------------------------------------
if ~exist(evFile,'file')
    % Read log files and extract events/pulses
    events=struct;
    pulses=[];
    [pulses,events] = parseTrainCSV(subj,sessDir,events,pulses);
    [pulses,events] = parsePracCSV(subj,sessDir,events,pulses);
    [pulses,events] = parseRunCSVs(subj,sessDir,events,pulses);
    if isempty(events)
        fprintf('BAD:  EMPTY EVENTS... exiting.\n\n')
        return
    end
    save(evFile,'events');
    clear events
    fprintf('DONE.\n')
else
  fprintf('SKIPPING (events exist).\n')
end

%--------------------------------------------
fprintf('  Making CCDT EEGlog for ')
fprintf('%s, %s: \n',subj,session)
if ~exist(eegLogFile,'file')
    if ~exist('pulses','var')
        events=struct;
        pulses=[];
        [pulses,events] = parseTrainCSV(subj,sessDir,events,pulses);
        [pulses,events] = parsePracCSV(subj,sessDir,events,pulses);
        [pulses,events] = parseRunCSVs(subj,sessDir,events,pulses);
    end
    buildEEGlog(pulses,eegLogFile);
else
    fprintf('SKIPPING (EEGlog exists).\n')
end
clear
fprintf('DONE.\n')
end

function [pulses,events] = parseTrainCSV(subj,sessDir,events,pulses)
    if isempty(fieldnames(events))
        i=0;
    else
        i=length(events);
    end
    evInd = [i+1 i+1];
    sessNum = strsplit(sessDir,'/');
    sessNum = str2double(regexp(sessNum{length(sessNum)},'\d*','Match'));
    trainFile = dir(fullfile(sessDir,'subj*_log_trainSoc.csv'));
    fid = fopen(fullfile(sessDir,trainFile.name));
    if fid ==-1,error(['log file not found: ' trainFile.name]);end
    thisLine = fgetl(fid);%skip header
    thisLine = fgetl(fid);%first data line
    while ischar(thisLine)
        logLine = strsplit(thisLine,',');
        if strcmp(logLine{1},'extraInfo')
            thisLine = fgetl(fid);% next line
            while ischar(thisLine)
                logLine = strsplit(thisLine,',');
                if strcmp(logLine{1},'START')
                    shiftTimes = str2double(logLine{2})/1000;
                    thisLine = fgetl(fid);% next line
                    thisLine = fgetl(fid);% next line
                    thisLine = fgetl(fid);% next line
                end
                thisLine = fgetl(fid);% next line
            end
            continue
        end
        [~,rspInd]=ismember({'[''f'']','[''j'']'},logLine);
        numAttempts = str2double(logLine{min(rspInd(rspInd>0))-1});
        pulseStart = 8+numAttempts*2;
        pulseEnd = min(rspInd(rspInd>0))-2;
        for ii=pulseStart:pulseEnd
            pulses = [pulses;str2double(logLine{ii})];
        end
        if strcmp(logLine{5},'TRUE')
            target = 'j';
        else
            target = 'f';
        end
        correct=0;
        for j=1:numAttempts
            if j==numAttempts,correct=1;end
            % Onset 1
            i=i+1;
            events(i).subj = subj;
            events(i).pid = str2double(logLine{4});
            events(i).session = sessNum;
            events(i).section = 'train';
            events(i).trial = str2double(logLine{3})+1;
            events(i).image = str2double(regexp(logLine{2},'\d*','Match'));
            events(i).type = 'onset1';
            events(i).target = target;
            events(i).response = char(regexp(logLine{min(rspInd(rspInd>0))+j-1},'\w','Match'));
            events(i).correct = correct;
            events(i).sound = 0;
            events(i).walk = str2double(logLine{1});
            events(i).mstime = str2double(logLine{8+j-1});
            % Response 1
            i=i+1;
            events(i).subj = subj;
            events(i).pid = str2double(logLine{4});
            events(i).session = sessNum;
            events(i).section = 'train';
            events(i).trial = str2double(logLine{3})+1;
            events(i).image = str2double(regexp(logLine{2},'\d*','Match'));
            events(i).type = 'response1';
            events(i).target = target;
            events(i).response = char(regexp(logLine{min(rspInd(rspInd>0))+j-1},'\w','Match'));
            events(i).correct = correct;
            events(i).sound = 0;
            events(i).walk = str2double(logLine{1});
            events(i).mstime = str2double(logLine{8+j-1})+str2double(logLine{min(rspInd(rspInd>0))+numAttempts+j-1});
            % Onset 2
            i=i+1;
            events(i).subj = subj;
            events(i).pid = str2double(logLine{4});
            events(i).session = sessNum;
            events(i).section = 'train';
            events(i).trial = str2double(logLine{3})+1;
            events(i).image = str2double(regexp(logLine{2},'\d*','Match'));
            events(i).type = 'onset2';
            events(i).target = target;
            events(i).response = char(regexp(logLine{min(rspInd(rspInd>0))+j-1},'\w','Match'));
            events(i).correct = correct;
            events(i).sound = 0;
            events(i).walk = str2double(logLine{1});
            events(i).mstime = str2double(logLine{8+numAttempts+j-1});
            % Response 2
            i=i+1;
            events(i).subj = subj;
            events(i).pid = str2double(logLine{4});
            events(i).session = sessNum;
            events(i).section = 'train';
            events(i).trial = str2double(logLine{3})+1;
            events(i).image = str2double(regexp(logLine{2},'\d*','Match'));
            events(i).type = 'response2';
            events(i).target = target;
            events(i).response = char(regexp(logLine{min(rspInd(rspInd>0))+j-1},'\w','Match'));
            events(i).correct = correct;
            events(i).sound = 0;
            events(i).walk = str2double(logLine{1});
            events(i).mstime = str2double(logLine{8+numAttempts+j-1})+str2double(logLine{min(rspInd(rspInd>0))+2*numAttempts+j-1});
        end
        thisLine = fgetl(fid);% next line
    end
    evInd(2) = i;
    for i=evInd(1):evInd(2)
        events(i).mstime = events(i).mstime+shiftTimes;
    end
    pulses = pulses+shiftTimes;
    pulses(1)
end

function [pulses,events] = parsePracCSV(subj,sessDir,events,pulses)
    if isempty(fieldnames(events))
        i=0;
    else
        i=length(events);
    end
    evInd = [i+1 i+1];
    pInd = length(pulses)+1;
    sessNum = strsplit(sessDir,'/');
    sessNum = str2double(regexp(sessNum{length(sessNum)},'\d*','Match'));
    pracFile = dir(fullfile(sessDir,'subj*_log_pracSoc.csv'));
    fid = fopen(fullfile(sessDir,pracFile.name));
    if fid ==-1,error(['log file not found: ' pracFile.name]);end
    thisLine = fgetl(fid);%skip header
    thisLine = fgetl(fid);%first data line
    while ischar(thisLine)
        logLine = strsplit(thisLine,',');
        if strcmp(logLine{1},'extraInfo')
            thisLine = fgetl(fid);% next line
            while ischar(thisLine)
                logLine = strsplit(thisLine,',');
                if strcmp(logLine{1},'START')
                    shiftTimes = str2double(logLine{2})/1000;
                    thisLine = fgetl(fid);% next line
                    thisLine = fgetl(fid);% next line
                    thisLine = fgetl(fid);% next line
                end
                thisLine = fgetl(fid);% next line
            end
            continue
        end
        [~,rspInd]=ismember({'''f''','''j'''},logLine);
        pulseStart = 8;
        pulseEnd = min(rspInd(rspInd>0))-1;
        for ii=pulseStart:pulseEnd
            pulses = [pulses;str2double(logLine{ii})];
        end
        if str2double(logLine{5})
            target = 'j';
        else
            target = 'f';
        end
        if str2double(logLine{6})
            sound = 0;
        else
            sound = 1;
        end
        % Onset 1
        i=i+1;
        events(i).subj = subj;
        events(i).pid = str2double(logLine{2});
        events(i).session = sessNum;
        events(i).section = 'practice';
        events(i).trial = str2double(logLine{4})+1;
        events(i).image = str2double(regexp(logLine{1},'\d*','Match'));
        events(i).type = 'onset';
        events(i).target = target;
        events(i).response = char(regexp(logLine{min(rspInd(rspInd>0))},'\w','Match'));
        events(i).correct = str2double(logLine{6});
        events(i).sound = sound;
        events(i).walk = str2double(logLine{3});
        events(i).mstime = str2double(logLine{7});
        % Response 1
        i=i+1;
        events(i).subj = subj;
        events(i).pid = str2double(logLine{2});
        events(i).session = sessNum;
        events(i).section = 'practice';
        events(i).trial = str2double(logLine{4})+1;
        events(i).image = str2double(regexp(logLine{1},'\d*','Match'));
        events(i).type = 'response';
        events(i).target = target;
        events(i).response = char(regexp(logLine{min(rspInd(rspInd>0))},'\w','Match'));
        events(i).correct = str2double(logLine{6});
        events(i).sound = sound;
        events(i).walk = str2double(logLine{3});
        events(i).mstime = str2double(logLine{7})+str2double(logLine{min(rspInd(rspInd>0))+1});
        
        thisLine = fgetl(fid);% next line
    end
    evInd(2) = i;
    for i=evInd(1):evInd(2)
        events(i).mstime = events(i).mstime+shiftTimes;
    end
    pulses(pInd:end) = pulses(pInd:end)+shiftTimes;
    pulses(pInd)
end

function [pulses,events] = parseRunCSVs(subj,sessDir,events,pulses)
    if isempty(fieldnames(events))
        i=0;
    else
        i=length(events);
    end
    sessNum = strsplit(sessDir,'/');
    sessNum = str2double(regexp(sessNum{length(sessNum)},'\d*','Match'));
    runFiles = dir(fullfile(sessDir,'subj*_logSoc_run*.csv'));
    for j=1:length(runFiles)
        fid = fopen(fullfile(sessDir,runFiles(j).name));
        if fid ==-1,error(['log file not found: ' runFiles(j).name]);end
        evInd = [i+1 i+1];
        pInd = length(pulses)+1;
        thisLine = fgetl(fid);%skip header
        thisLine = fgetl(fid);%first data line
        while ischar(thisLine)
            logLine = strsplit(thisLine,',');
            if ismember(logLine{1},{'extraInfo',''})
                thisLine = fgetl(fid);% next line
                while ischar(thisLine)
                    logLine = strsplit(thisLine,',');
                    if strcmp(logLine{1},'START')
                        shiftTimes = str2double(logLine{2})/1000;
                        thisLine = fgetl(fid);% next line
                        thisLine = fgetl(fid);% next line
                        thisLine = fgetl(fid);% next line
                    end
                    thisLine = fgetl(fid);% next line
                end
                continue
            end
            pulses = [pulses;[str2double(logLine{8});str2double(logLine{9})]];
            if isnan(str2double(logLine{11}))&& ~strcmp(logLine{11}(2:end-1),'NA')% For some reason, RT is sometimes saved as a string (within string), still appears accurate
                logLine{11} = logLine{11}(2:end-1);
            end
            if strcmp(logLine{5},'TRUE')
                target = 'j';
            else
                target = 'f';
            end
            if str2double(logLine{6})
                sound = 0;
            elseif strcmp(char(regexp(logLine{10},'\w','Match'))','NA')
                sound = 2;%non response, low freq tone
                logLine{11} = '2';
            else
                sound =1;%incorrect response, high freq tone
            end
            % Onset 1
            i=i+1;
            events(i).subj = subj;
            events(i).pid = str2double(logLine{2});
            events(i).session = sessNum;
            events(i).section = ['run' num2str(j)];
            events(i).trial = str2double(logLine{4})+1;
            events(i).image = str2double(regexp(logLine{1},'\d*','Match'));
            events(i).type = 'onset';
            events(i).target = target;
            events(i).response = char(regexp(logLine{10},'\w','Match'))';
            events(i).correct = str2double(logLine{6});
            events(i).sound = sound;
            events(i).walk = str2double(logLine{3});
            events(i).mstime = str2double(logLine{7});
            % Response 1
            i=i+1;
            events(i).subj = subj;
            events(i).pid = str2double(logLine{2});
            events(i).session = sessNum;
            events(i).section = ['run' num2str(j)];
            events(i).trial = str2double(logLine{4})+1;
            events(i).image = str2double(regexp(logLine{1},'\d*','Match'));
            events(i).type = 'response';
            events(i).target = target;
            events(i).response = char(regexp(logLine{10},'\w','Match'))';
            events(i).correct = str2double(logLine{6});
            events(i).sound = sound;
            events(i).walk = str2double(logLine{3});
            events(i).mstime = str2double(logLine{7})+str2double(logLine{11});
            
            thisLine = fgetl(fid);%next line
        end
        evInd(2) = i;
        for i=evInd(1):evInd(2)
            events(i).mstime = events(i).mstime+shiftTimes;
        end
        pulses(pInd:end) = pulses(pInd:end)+shiftTimes;
        pulses(pInd)
    end
end

function buildEEGlog(syncTimes,eegLogFile)
% Write eeg.eeglog.up file similar to PyEPL processing
fid = fopen(eegLogFile,'w');
for i=1:length(syncTimes)    
    fprintf(fid,'%s\t%d\t%s\n',num2str(syncTimes(i,1)),1,'UP');
end
fclose(fid);
end