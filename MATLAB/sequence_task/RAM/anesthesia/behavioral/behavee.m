function events = behavee(ddir,subj)
dobehave = 1;
cd([ddir filesep subj filesep 'lfp.noreref']);
fname = ls('*.001');
[~,subjF,~] = fileparts(fname); % automatically obtain the correct filename for electrodes


sr = eegparams('samplerate',[ddir  filesep subj filesep 'lfp.noreref' filesep 'params.txt']); %obtain sample rate from params text file
timefix = sr*60; % convert sample rate to samples/min
% fs = sr;

events = [];
load([ddir filesep subj filesep 'behavioral' filesep 'session_0' filesep 'events.mat']);
rootF = [ddir filesep subj filesep 'lfp.noreref' filesep];
etime = events(length(events)-1).lfpoffset/timefix; %obtain length of behavioral event trials (last trial is usually incomplete, hence length(events)-1
% etime = (events(length(events)-1).lfpoffset-events(1).lfpoffset)/timefix;
% etime = 11;
% subj = 'HUP119_i_24Jun16_0959'

for n=1:length(events)
    events(n).lfpfile = [rootF subjF];
end


for n=1:length(events)
    events(n).lfpoffset2 = events(n).lfpoffset/timefix;
end
if dobehave == 1
    
% [correctR incorrectR] = eventStrategize(events); % find correct, incorrect, and pass (LOC) recall trials
% cor = [correctR(1).lfpoffset/timefix correctR(size(correctR,2)).lfpoffset/timefix];
% LWM = incorrectR(1).lfpoffset/timefix;
% LOC = correctR(size(correctR,2)).lfpoffset/timefix; % estmated based on final correct trial

ci = 1;
ii = 1;
rti = 1;
rti1 = 1;
resorte = events(1).lfpoffset2; %0; % 
% sr = 2.047999e+03 %sample rate from 
yRT1 = []; xpc = []; ypc = []; xpcR = []; ypcR = []; yRT = []; yRT2 = []; xRT1 = []; yRT3 = [];
for i = 1:length(events)
    if i+2<length(events) && strcmp(events(i).type,'RECALL') && strcmp(events(i+2).type,'RECALL');
        events(i).correct = strcmp(events(i).response,events(i).target);
        events(i+1).correct = strcmp(events(i+1).response,events(i+1).target);
        events(i+2).correct = strcmp(events(i+2).response,events(i+2).target);
        events(i).percentcorrect = (events(i).correct + events(i+2).correct + events(i+1).correct)/3*100;
        events(i+1).percentcorrect = events(i).percentcorrect;
        events(i+2).percentcorrect = events(i).percentcorrect;
%         pc(ci) = events(i);
%         pc(ci+1) = events(i+1);
%         pc(ci+2) = events(i+2);
        xpc(ci) = events(i).lfpoffset2-resorte;
        xpc(ci+1) = events(i+1).lfpoffset2-resorte;
        xpc(ci+2) = events(i+2).lfpoffset2-resorte;
        ypc(ci) = events(i).percentcorrect;
        ypc(ci+1) = events(i+1).percentcorrect;
        ypc(ci+2) = events(i+2).percentcorrect;
        events(i).RT1 = events(i+1).lfpoffset-events(i).lfpoffset;
        events(i).RT2 = events(i+2).lfpoffset-events(i+1).lfpoffset;
        xRT1(rti1) = events(i).lfpoffset2-resorte;
        yRT1(rti1) = mean([events(i).RT1 events(i).RT2])/sr;
        rti1 = rti1+1;
        ci = ci + 3;   

        end
        if i+4<length(events) && strcmp(events(i).type,'RESP') && strcmp(events(i+2).type,'RESP') && strcmp(events(i+4).type,'RESP');
        events(i).correctRESP = strcmp(events(i).response,events(i).target);
        events(i+2).correctRESP = strcmp(events(i+2).response,events(i+2).target);
        events(i+4).correctRESP = strcmp(events(i+4).response,events(i+4).target);
        events(i).percentcorrectRESP = (events(i).correctRESP + events(i+2).correctRESP + events(i+4).correctRESP)/3*100;
        events(i+2).percentcorrectRESP = events(i).percentcorrectRESP;
        events(i+4).percentcorrectRESP = events(i).percentcorrectRESP;
%         pcR(ii) = events(i);
%         pcR(ii+1) = events(i+2);
%         pcR(ii+2) = events(i+4);
        xpcR(ii) = events(i).lfpoffset2-resorte;
        xpcR(ii+1) = events(i+2).lfpoffset2-resorte;
        xpcR(ii+2) = events(i+4).lfpoffset2-resorte;
        ypcR(ii) = events(i).percentcorrectRESP;
        ypcR(ii+1) = events(i+2).percentcorrectRESP;
        ypcR(ii+2) = events(i+4).percentcorrectRESP;
        ii = ii + 3;   
        end
        if i+1<length(events) && strcmp(events(i).type,'SOUND') && strcmp(events(i+1).type,'RESP')
            events(i).RT = events(i+1).lfpoffset-events(i).lfpoffset;
            xRT(rti) = events(i).lfpoffset2-resorte;
            yRT(rti) = events(i).RT/sr;
            rti = rti + 1;
        end
        
end

    yRT2 = yRT/max(yRT)*100; % percent maximal RT
    if yRT1
    yRT3 = yRT1/max(yRT1)*100; % percent maximal RT
    end
    RTmax = find(yRT3 > 90);
    
%     seq = RTmax(diff(RTmax) == 1);
%     ret = xRT1(seq(1)); % time of max RT for resorting

    seq = RTmax(diff(RTmax) == 1);
    if seq(1) == 1
    ret = xRT1(seq(2)); % time of max RT for resorting
    else
    ret = xRT1(seq(1)); % time of max RT for resorting 
    end
    
    retsr = ret*timefix;
end
events = events(1); % take only the first event to look at entire time sequence.
events.lfpoffset = events.lfpoffset+round(retsr);
end