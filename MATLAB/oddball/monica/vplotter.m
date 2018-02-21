function [retain1, retain2] = vplotter(events,comparison)

switch comparison 
    case 'backgroundTARGET'
        retain1 = strcmp({events.type}, 'BACKGROUND');
        retain2 = strcmp({events.type}, 'TARGET');  
    case 'passiveEARLY'
%        [passEvents] = splitevents(events);
        earlyThresh = round((length(events)/3));
        lateThresh = earlyThresh*2;
%         retain1 = strcmp({passEvents(1:earlyThresh).type},'TARGET');
%         retain2 = strcmp({passEvents(1:earlyThresh).type},'BACKGROUND');
        retain1 = strcmp({events(1:earlyThresh).type},'TARGET');
        retain2 = strcmp({events(1:earlyThresh).type},'BACKGROUND');
    case 'passiveLATE'
 %       [passEvents] = splitevents(events);
        earlyThresh = round((length(events)/3));
       % etot = length(passEvents);
        etot = length(events);
        lateThresh = earlyThresh*2;
%         retain1 = strcmp({passEvents(lateThresh:etot).type},'TARGET');
%         retain2 = strcmp({passEvents(lateThresh:etot).type},'BACKGROUND');     
        retain1 = strcmp({events(lateThresh:etot).type},'TARGET');
        retain2 = strcmp({events(lateThresh:etot).type},'BACKGROUND');      
    case 'goNOGO'
        retain1 = strcmp({events.type},'GO')&cell2mat({events.response})==1;
        retain2 = strcmp({events.type},'NO_GO')&cell2mat({events.response})==0;
    case 'goEARLY'
        earlyThresh = round((length(events)/3));
        lateThresh = earlyThresh*2;
        retain1 = strcmp({events(1:earlyThresh).type},'GO');
        retain2 = strcmp({events(1:earlyThresh).type},'NO_GO');
    case 'goLATE'
        earlyThresh = round((length(events)/3));
        etot = length(events);
        lateThresh = earlyThresh*2;  
        retain1 = strcmp({events(lateThresh:etot).type},'GO');
        retain2 = strcmp({events(lateThresh:etot).type},'NO_GO');   
    case 'goRT'
        retain1 = strcmp({events.type},'GO')&cell2mat({events.response})==1;
        retain2 = strcmp({events.type},'RESPONSE')&cell2mat({events.response})==1;
    case 'nogoERRORS'
        retain1 = strcmp({events.type}, 'NO_GO')&cell2mat({events.response})==0;
        retain2 = strcmp({events.type}, 'NO_GO')&cell2mat({events.response})==1;
    case 'goERRORS'
        retain1 = strcmp({events.type}, 'GO')&cell2mat({events.response})==0;
        retain2 = strcmp({events.type}, 'GO')&cell2mat({events.response})==1; 
    case 'targetERRORS'
        retain1 = strcmp({events.type}, 'TARGET')&cell2mat({events.response})==0;
        retain2 = strcmp({events.type}, 'TARGET')&cell2mat({events.response})==1;
    case 'targRESP'
        retain1 = strcmp({events.type}, 'TARGET')&cell2mat({events.response})==0;
        retain2 = strcmp({events.type}, 'RESPONSE')&cell2mat({events.response})==1;
    case 'passiveGONOGO'
        retain1 = strcmp({events.type}, 'END_SOUND')&cell2mat({events})<=100;
        retain2 = strcmp({events.type}, 'RESPONSE_START');
end


