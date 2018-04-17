function nl_example
%Inputs:
subj = 'HUP157';
evFile = '/Volumes/HumanStudies/HumanStudies/netLearn/eeg/HUP157/behavioral/Session_0/events.mat';
eegDir = '/Volumes/HumanStudies/HumanStudies/netLearn/eeg/HUP157/eeg.noreref';
elec = 18;
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
figure;
plot(nanmean(eeg,1))
