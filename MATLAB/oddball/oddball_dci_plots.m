%% Make plots from oddball_dci.m saved output
% Inputs:
subj = 'HUP155_i';
ddir = '/Volumes/HumanStudies/HumanStudies/oddball';
tw = '5462';
order = '1';

%% load in all related .mat files
p_dir = fullfile(ddir,'eeg',subj,'processed');
dci_dir = fullfile(ddir,'scratch/dci',subj,['tw' tw],'Analysis');
load(fullfile(p_dir,'sessInfo.mat'));
load(fullfile(p_dir,'parameters.mat'));
load(fullfile(dci_dir,'grid_win.mat'));
load(fullfile(dci_dir,['order_' order],'AR_mod.mat'));
load(fullfile(dci_dir,['order_' order],'eig_modes.mat'));
load(fullfile(dci_dir,['order_' order],'medians.mat'));
load(fullfile(dci_dir,['order_' order],'time_s.mat'));

%% Behavioral analysis
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

% scatter(cell2mat({events(correct_go).lfpoffset}),rts(correct_go),'ob')
% hold on
% scatter(cell2mat({events(incorrect_go).lfpoffset}),rts(incorrect_go),'xr')
%% Median Eigenmode plot
figure;
subplot(2,1,1)
scatter(cell2mat({events(correct_go).lfpoffset}),rts(correct_go),'ob');
hold on
scatter(cell2mat({events(incorrect_go).lfpoffset}),rts(incorrect_go),'xr');
xlabel('Time (samples)', 'fontsize', 15);
ylabel('Reaction time (ms)', 'fontsize', 15);
title('''Go'' Responses', 'fontsize', 20);
legend('Correct','Incorrect');
subplot(2,1,2)
nlxEventsToPlot = [3,4,5,6,8,9,11,12,19]; % set manually for each subject!
plot(time_s, med);
hold on
nlxEvHandles = [];
for i=1:length(nlxEventsToPlot)
    ts = time_s(floor(nlxEvents(nlxEventsToPlot(i)).lfpoffset/win));
    h=plot([ts, ts], [min(med), max(med)], 'linewidth', 2);
    nlxEvHandles = [nlxEvHandles,h];
end
nlxEvLabels = {nlxEvents(nlxEventsToPlot).type};
xlim([time_s(1), time_s(end)])
ylim([min(med), max(med)])
xlabel('Time (windows)', 'fontsize', 15)
ylabel('Median', 'fontsize', 15);
title('Eigenmode Median', 'fontsize', 20)
legend(nlxEvHandles,nlxEvLabels)
%% Next Plot