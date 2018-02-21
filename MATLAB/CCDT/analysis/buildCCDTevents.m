function buildCCDTevents
% builds events.mat and eeg.eeglog.up for CCDT session. Stores in session folder.
% events structure with fields:
% subj,session,trial#,eventType,item,msTime,msOffset
subj = 'HUP157';
session = 'Session_0';

dirs = le_dirs('CCDT');
sessDir = fullfile(dirs.data,'eeg',subj,'behavioral',session);
evFile = fullfile(sessDir,'events.mat');
eegLogFile = fullfile(sessDir,'eeg.eeglog.up');

% print opening line
fprintf('  Making CCDT EVENTS for ')
fprintf('%s, %s: \n',subj,session)

% Get RTs and pulses
[allRTs,allSyncs] = CCDTcatBehavioral(subj,{session});

%--------------------------------------------
if ~exist(evFile,'file')
  events=buildCCDTevents_local(subj,session,allRTs);
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
    buildEEGlog(allSyncs,eegLogFile);
else
    fprintf('SKIPPING (EEGlog exists).\n')
end
clear
fprintf('DONE.\n')
end

function events=buildCCDTevents_local(subj,session,sessRTs)
for trial=1:length(sessRTs)
    fixMS = sessRTs(trial,2) - sessRTs(trial,1) - sessRTs(trial,3);
    ccMS = sessRTs(trial,2) - sessRTs(trial,1);
    respMS = sessRTs(trial,2);
    setNum = sessRTs(trial,4);
    events(trial*3-2).subject = subj;
    events(trial*3-1).subject = subj;
    events(trial*3).subject = subj;
    events(trial*3-2).session = session;
    events(trial*3-1).session = session;
    events(trial*3).session = session;
    events(trial*3-2).set = setNum;
    events(trial*3-1).set = setNum;
    events(trial*3).set = setNum;
    events(trial*3-2).trial = trial;
    events(trial*3-1).trial = trial;
    events(trial*3).trial = trial;
    events(trial*3-2).type = 'FIX_START';
    events(trial*3-1).type = 'CC';
    events(trial*3).type = 'RESPONSE';
    events(trial*3-2).rt = sessRTs(trial,1);
    events(trial*3-1).rt = sessRTs(trial,1);
    events(trial*3).rt = sessRTs(trial,1);
    events(trial*3-2).delay = sessRTs(trial,3);
    events(trial*3-1).delay = sessRTs(trial,3);
    events(trial*3).delay = sessRTs(trial,3);
    events(trial*3-2).mstime = fixMS;
    events(trial*3-1).mstime = ccMS;
    events(trial*3).mstime = respMS;
    events(trial*3-2).msoffset = 0;
    events(trial*3-1).msoffset = 0;
    events(trial*3).msoffset = 0;
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