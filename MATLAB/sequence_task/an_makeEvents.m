function[events] = an_makeEvents(subj,sess)
% This function makes events for a particular session
% e.g., an_makeEvents('NR_062415','post')

% go to beh dir
[~,dataDir] = an_subjInfo(subj);
cd(fullfile(dataDir,subj,'behavioral',strcat('session_',num2str(sess))));

% Read in session.log file.
fid=fopen('session.log');

events=[];
count=0;
trialCount = 0; % keep running gount of trial number
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
        case {'TRIAL_START'}
            trialCount = trialCount+1;
            
        case {'SOUND','RESP','RECALL','RECALL_END'}
         %Key line of code: Read events.
         c=textscan(tline,' %f%d%s%s%s%s','delimiter','\t','returnOnError',0); 
           
            %Make event
            events(count).trial=trialCount;
            events(count).seq_pos=str2double(c{4}{1});
            
            % process target (sound/command/target)
            if ~isempty(c{5})
                if strfind(c{5}{1},'right')
                    events(count).target= 'right';
                elseif strfind(c{5}{1},'left')
                    events(count).target= 'left';
                else
                    events(count).target= '';
                end
            else
                events(count).target = '';
            end
            
            % process response (response/recall)
            if ~isempty(c{6})
                if strfind(c{6}{1},'right')
                    events(count).response= 'right';
                elseif strfind(c{6}{1},'left')
                    events(count).response= 'left';
                elseif strfind(c{6}{1},'pass')
                    events(count).response= 'pass';
                else
                    events(count).response= '';
                end
            else
                events(count).response = '';
            end
            
            % process feedback (special case)
            if strcmp(events(count).type,'RECALL_END')
                events(count).type = 'feedback'; 
                events(count).seq_pos = nan;
                if strfind(c{4}{1},'isCorrect=1')
                    events(count).feedback = 1;
                elseif strfind(c{4}{1},'isCorrect=0')
                    events(count).feedback = 0;
                else
                    events(count).feedback = '';
                end
            else
                events(count).feedback = '';
            end
            
    otherwise
        %fill in empty quotations
        events(count).type=''; 
        events(count).trial=''; 
        events(count).seq_pos='';
        events(count).target='';
        events(count).response='';
        events(count).feedback='';
    end
end     

% save events 
%within session folder
save('events','events')