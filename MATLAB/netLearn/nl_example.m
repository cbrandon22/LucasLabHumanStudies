function nl_example
%Inputs:
subj = 'HUP143';
evFile = '/Volumes/LUCAS_DRIVE/HumanStudies/netLearn/events/HUP143_events.mat';
eegDir = '/Volumes/LUCAS_DRIVE/HumanStudies/netLearn/eeg/HUP143/eeg.noreref';
elec = 2;
signalStartMS = -200;
signalDuration = 800;

%Load events.mat
load(evFile,'events');

%Select events to average
trainEv = strcmp({events.section},'train') & cell2mat({events.correct})==1;
pracEv = strcmp({events.section},'practice');
run1Ev = strcmp({events.section},'run1');
run2Ev = strcmp({events.section},'run1');
allRunEv = logical(run1Ev+run2Ev);

runOnsets = allRunEv & strcmp({events.type},'onset');
%behavioral examples

%Use boolean array to index into eeg
eeg = gete_ms_wrapper(elec,events(runOnsets),signalDuration,signalStartMS,[],[],[],[]);
plot(nanmean(eeg,1))
