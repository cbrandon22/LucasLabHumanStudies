function[events] = an_makeEvents2(subj,sess)
% This function makes events for a particular session
% e.g., an_makeEvents('NR_062415','post')

% go to beh dir
[~,dataDir] = an_subjInfo2(subj);
% Read config
fid=fopen(fullfile(dataDir,subj,'behavioral','config.py'));
tline=fgetl(fid);
while ischar(tline)
    if strfind(tline,'flipTones')>0
        if strfind(tline,'True')>0
            flipTones = 1;
        elseif strfind(tline,'False')>0
            flipTones = 0;
        else
            disp('flipTones not found in config. Setting to true...');
            flipTones = 1;
        end
    end
    tline = fgetl(fid);
end
        
cd(fullfile(dataDir,subj,'behavioral',strcat('session_',num2str(sess))));

% Read in session.log file.
fid=fopen('session.log');

events=[];
count=0;
trialCount = 0; % keep running gount of trial number
passiveCount = 0;

%loop through each line and generate event entries
while true
    tline=fgetl(fid);  
    if ~ischar(tline)
        break
    end
    count=count+1;
    c=textscan(tline,'%f%d%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter','\t');
    events(count).subject=subj;
    events(count).session=sess;
    events(count).mstime=(c{1});
    events(count).type=c{3}{1};

    switch events(count).type
            
        case {'BACKGROUND','TARGET'}
            trialCount = trialCount+1;
            c=textscan(tline,'%f%f%s%d','delimiter','\t');
            events(count).trial=trialCount;
            if rBlock && strcmp(events(count).type,'BACKGROUND')
                events(count).type = 'NO_GO';
            elseif rBlock && strcmp(events(count).type,'TARGET')
                events(count).type = 'GO';
            elseif targTone ==1 && strcmp(events(count).type,'BACKGROUND')
                events(count).type = 'BACKGROUNDLF';
            elseif targTone ==1 && strcmp(events(count).type,'TARGET')
                events(count).type = 'TARGETHF';
            elseif targTone ==2 && strcmp(events(count).type,'BACKGROUND')
                events(count).type = 'BACKGROUNDHF';
            elseif targTone ==2 && strcmp(events(count).type,'TARGET')
                events(count).type = 'TARGETLF';
            end
            
        case {'END_SOUND','PASSIVE_START','RESPONSE_START','PASSIVE_END','RESPONSE_END'}
            %Key line of code: Read events.
            c=textscan(tline,'%f%f%s%d','delimiter','\t');
            if strcmp(events(count).type,'PASSIVE_START')
                passiveCount = passiveCount+1;
                rBlock = 0;
                targTone = 1;
                if flipTones && rem(passiveCount,2)==0
                    targTone = 2;
                end
            elseif strcmp(events(count).type,'RESPONSE_START')
                rBlock = 1;
                targTone = 0;
            end
           
            %Make event
            events(count).trial=trialCount;
            
        case 'RESPONSE'
            c=textscan(tline,'%f%f%s%d%s%s','delimiter','\t');
            respBool = regexp(c{6},'=','split');
            respBool = respBool{1}{2};
            if strcmp(respBool,'False')
                respBool = 0;
            elseif strcmp(respBool,'True')
                respBool = 1;
            end
            events(count).trial=trialCount;
            events(count).response=respBool;
            events(count-1).response=respBool;
            events(count-2).response=respBool;
            
    otherwise
        %fill in empty quotations
        events(count).trial='';
        events(count).trial='';
        events(count).response='';
    end
end     

% save events 
%within session folder
save('events','events')