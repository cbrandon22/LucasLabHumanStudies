clear
% Inputs:
subj = 'HUP155_i';

p_dir = fullfile('D:\TNL_Data\oddball\eeg',subj,'processed');
load(fullfile(p_dir,'sessInfo.mat'));

% nlxEventsToPlot = [16]; %set per subject

go_ev = ismember({events.type},'GO');
nogo_ev = ismember({events.type},'NOGO');
resp_ev = cell2mat({events.response});
resp_ev(isnan(resp_ev)) = 0;
correct_go = go_ev & cell2mat({events.response})==1;
incorrect_go = go_ev & ~(cell2mat({events.response})==1);
rts = zeros(length(events),1);
for i=1:length(events)
    if strcmp(events(i).type,'RESPONSE')
        rts(i-2:i) = (events(i).lfpoffset-events(i-2).lfpoffset)/(srate/1000);
    end
end
% tstart=events.mstime;
% time=events.mstime-tstart
%% make plot
figure;
go_trial = cell2mat({events(correct_go).trial});
nogo_trial = cell2mat({events(incorrect_go).trial});
scatter(go_trial,rts(correct_go),'ob');
hold on
scatter(nogo_trial,rts(incorrect_go),'xr');
ylabel('Reaction time (ms)', 'fontsize', 15);
% title('''Go'' Responses', 'fontsize', 20);
[lh,~]=legend('Correct','Incorrect','Location','southeast');
xlabel('Trials', 'fontsize', 15);