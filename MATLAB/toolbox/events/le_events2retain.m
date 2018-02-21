
function[retain1,retain2,retain3] = le_events2retain(events,task,alignment,comparison,thresh)
%This function maintains a master list of comparisons that can be performed
%when making reports
%Inputs
%events ... a set of events)
%task ... the task (e.g., 'motormap'
%comparison ... the type of comparison format is task_compar (eg, 'pyFR_SME'
            %(see below))
%Outputs
%retain1 ... the indices with the first subset of events
%retain2 ... the indices with the second subset of events (can be empty)

switch task
    case 'motormap'
        switch comparison
            %% motormap
            case 'moveWait'
                retain1 = strcmp({events.type},'MOVE');
                retain2 = strcmp({events.type},'WAIT');
            case 'moveInstruct'
                retain1 = strcmp({events.type},'MOVE');
                retain2 = strcmp({events.type},'INSTRUCT');
            case 'leftWait'
                retain1 = strcmp({events.type},'MOVE')&strcmp({events.item},'Left Hand');
                retain2 = strcmp({events.type},'WAIT')&strcmp({events.item},'Left Hand');
            case 'rightWait'
                retain1 = strcmp({events.type},'MOVE')&strcmp({events.item},'Right Hand');
                retain2 = strcmp({events.type},'WAIT')&strcmp({events.item},'Right Hand');
            case 'mouthWait'
                retain1 = strcmp({events.type},'MOVE')&strcmp({events.item},'Mouth and Tongue');
                retain2 = strcmp({events.type},'WAIT')&strcmp({events.item},'Mouth and Tongue');
            case 'handWait'
                retain1 = (strcmp({events.type},'MOVE')&strcmp({events.item},'Left Hand'))|(strcmp({events.type},'MOVE')&strcmp({events.item},'Right Hand'));
                retain2 = (strcmp({events.type},'WAIT')&strcmp({events.item},'Left Hand'))|(strcmp({events.type},'WAIT')&strcmp({events.item},'Right Hand'));
            case 'instructWait'
                retain1 = strcmp({events.type},'INSTRUCT');
                retain2 = strcmp({events.type},'WAIT');
        end
    case 'CCDT'
        switch alignment
            case 'cue'
                switch comparison
                    case 'delay'
                        retain1 = strcmp({events.type},'FIX_START')&cell2mat({events.delay})==500;
                        %retain2 = strcmp({events.type},'FIX_START')&cell2mat({events.delay})==1000;
                        retain2 = strcmp({events.type},'FIX_START')&cell2mat({events.delay})==1500;
                    case 'RTshort'
                        retain1 = strcmp({events.type},'FIX_START')&cell2mat({events.rt})<=thresh(1) & cell2mat({events.delay})==500;
                        retain2 = strcmp({events.type},'FIX_START')&cell2mat({events.rt})>thresh(2) & cell2mat({events.delay})==500;
                    case 'RTlong'
                        retain1 = strcmp({events.type},'FIX_START')&cell2mat({events.rt})<=thresh(1) & cell2mat({events.delay})==1500;
                        retain2 = strcmp({events.type},'FIX_START')&cell2mat({events.rt})>thresh(2) & cell2mat({events.delay})==1500;
                    case 'timeShort'
                        earlyThresh = round((length(events)/3)/3);
                        lateThresh = (length(events)/3)-earlyThresh;
                        retain1 = strcmp({events.type},'FIX_START')&cell2mat({events.trial})<=earlyThresh&cell2mat({events.delay})==500;
                        retain2 = strcmp({events.type},'FIX_START')&cell2mat({events.trial})>=lateThresh&cell2mat({events.delay})==500;
                    case 'timeLong'
                        earlyThresh = round((length(events)/3)/3);
                        lateThresh = (length(events)/3)-earlyThresh;
                        retain1 = strcmp({events.type},'FIX_START')&cell2mat({events.trial})<=earlyThresh&cell2mat({events.delay})==1500;
                        retain2 = strcmp({events.type},'FIX_START')&cell2mat({events.trial})>=lateThresh&cell2mat({events.delay})==1500;
                end
            case 'CC'
                switch comparison
                    case 'delay'
                        retain1 = strcmp({events.type},'CC');
                    case 'RT'
                        retain1 = strcmp({events.type},'CC');
                    case 'time'
                        retain1 = strcmp({events.type},'CC');
                end
            case 'response'
                switch comparison
                    case 'delay'
                        retain1 = strcmp({events.type},'RESPONSE')&cell2mat({events.delay})==500;
                        retain2 = strcmp({events.type},'RESPONSE')&cell2mat({events.delay})==1500;
                    case 'RT'
                        retain1 = strcmp({events.type},'RESPONSE');
                    case 'timeLong'
                        earlyThresh = round((length(events)/3)/3); %first and last 3rd of trails (only 1500ms delay)
                        lateThresh = earlyThresh*2;
                        retain1 = strcmp({events.type},'RESPONSE')&cell2mat({events.trial})<=earlyThresh&cell2mat({events.delay})==1500;
                        retain2 = strcmp({events.type},'RESPONSE')&cell2mat({events.trial})<=lateThresh&cell2mat({events.delay})==1500;
                    case 'timeShort'
                        earlyThresh = round((length(events)/3)/3); %first and last 3rd of trails (only 500ms delay)
                        lateThresh = earlyThresh*2;
                        retain1 = strcmp({events.type},'RESPONSE')&cell2mat({events.trial})<=earlyThresh&cell2mat({events.delay})==500;
                        retain2 = strcmp({events.type},'RESPONSE')&cell2mat({events.trial})<=lateThresh&cell2mat({events.delay})==500;
                end
        end
    case('oddball')
        switch(alignment)
            case 'startSound'
                switch(comparison)
                    case 'TARGETHFBACKGROUNDHF'
                        retain1 = strcmp({events.type},'TARGETHF');
                        retain2 = strcmp({events.type},'BACKGROUNDHF');
                    case 'TARGETHFBACKGROUNDLF'
                        retain1 = strcmp({events.type},'TARGETHF');
                        retain2 = strcmp({events.type},'BACKGROUNDLF');
                    case 'TARGETLFBACKGROUNDLF'
                        retain1 = strcmp({events.type},'TARGETLF');
                        retain2 = strcmp({events.type},'BACKGROUNDLF');
                    case 'TARGETLFBACKGROUNDHF'
                        retain1 = strcmp({events.type},'TARGETLF');
                        retain2 = strcmp({events.type},'BACKGROUNDHF');
                    case 'GONOGO'
                        retain1 = strcmp({events.type},'GO');
                        retain2 = strcmp({events.type},'NOGO');
                end
        end
end    
end