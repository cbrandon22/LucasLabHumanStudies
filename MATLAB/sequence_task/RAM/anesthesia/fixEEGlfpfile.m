% This script simply changes the lfpfile field of events
% su = 'HUP119_i';
rootF = ['E:\data\C\data\eeg' filesep subj filesep 'lfp.noreref' filesep];
% subj = 'HUP119_i_24Jun16_0959'

for n=1:length(events)
    events(n).lfpfile = [rootF subjF];
end
