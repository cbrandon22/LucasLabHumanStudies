%% Make plots from oddball_dci.m saved output
clear
% Inputs:
subj = 'HUP155_i';
ddir = '/Volumes/HumanStudies/HumanStudies/oddball';
%tw = '1024';
%tw = '5462';
tw = '547';
decimate_factor = 10;
order = '1';
saveFigs = 0;
 
%% load in all related .mat files
p_dir = fullfile(ddir,'eeg',subj,'processed');
dci_dir = fullfile(ddir,'scratch/dci',subj,'lead_average',num2str(decimate_factor),['tw' tw],'Analysis');
load(fullfile(p_dir,'sessInfo.mat'));
load(fullfile(p_dir,'parameters.mat'));
load(fullfile(dci_dir,'grid_win.mat'));
load(fullfile(dci_dir,['order_' order],'AR_mod.mat'));
load(fullfile(dci_dir,['order_' order],'eig_modes.mat'));
load(fullfile(dci_dir,['order_' order],'medians.mat'));
load(fullfile(dci_dir,['order_' order],'time_s.mat'));
load(fullfile(dci_dir,['order_' order],'filtTime.mat'));
time_s = ((1:size(AR_mod,2)).*win - (win/2))./ceil(srate/decimate_factor);
load(fullfile(ddir,'scratch/dci/plotEventsStruct.mat'));
if ismember(subj,{plotEventsStruct.subj})
    [~, sInd] = ismember(subj,{plotEventsStruct.subj});
    nlxEventsToPlot = plotEventsStruct(sInd).events;
else
    disp('Select neuralynx events to plot for this patient')
    keyboard; % Edit the events field and run the next 3 lines:
    plotEventsStruct(length(plotEventsStruct)+1).subj = subj;
    plotEventsStruct(length(plotEventsStruct)).events = [7,8,14,16,18,19,20,21,24,26,27];
    save(fullfile(ddir,'scratch/dci/plotEventsStruct.mat'),'plotEventsStruct');
end
 
%% Behavioral analysis
if sum(ismember({events.type},'GO'))>1
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
end
 
% scatter(cell2mat({events(correct_go).lfpoffset}),rts(correct_go),'ob')
% hold on
% scatter(cell2mat({events(incorrect_go).lfpoffset}),rts(incorrect_go),'xr')
%% Median Eigenmode plot
if sum(ismember({events.type},'GO'))>1
    %behavioral subplot
    %t_min = (time_s*win)/srate/60;
    figure;
    ax1=subplot(3,1,1);
    if decimate_factor > 0
        go_ts = time_s(floor(ceil(cell2mat({events(correct_go).lfpoffset})/decimate_factor)/win))/60;
        nogo_ts = time_s(floor(ceil(cell2mat({events(incorrect_go).lfpoffset})/decimate_factor)/win))/60;
    else
        go_ts = time_s(floor(cell2mat({events(correct_go).lfpoffset})/win))/60;
        nogo_ts = time_s(floor(cell2mat({events(incorrect_go).lfpoffset})/win))/60;
    end
    scatter(go_ts,rts(correct_go),'ob');
    hold on
    scatter(nogo_ts,rts(incorrect_go),'xr');
    ylabel('Reaction time (ms)', 'fontsize', 15);
    title('''Go'' Responses', 'fontsize', 20);
    [lh,~]=legend('Correct','Incorrect','Location','southeast');
    
    % eigenmodes
    ax2=subplot(3,1,[2 3]);
else
    figure;
end
plot(time_s/60, med);
hold on
nlxEvHandles = [];
for i=1:length(nlxEventsToPlot)
    if decimate_factor > 0
        ts = time_s(floor(ceil(nlxEvents(nlxEventsToPlot(i)).lfpoffset/decimate_factor)/win))/60;
    else
        ts = time_s(floor(nlxEvents(nlxEventsToPlot(i)).lfpoffset/win))/60;
    end
    h=plot([ts, ts], [min(med), max(med)], 'linewidth', 2);
    nlxEvHandles = [nlxEvHandles,h];
end
nlxEvLabels = {nlxEvents(nlxEventsToPlot).type};
xlim([time_s(1)/60, time_s(end)/60])
ylim([min(med), max(med)])
xlabel('Time (minutes)', 'fontsize', 15)
ylabel('Median', 'fontsize', 15);
win_ms = num2str(round(win/srate*1000));
title(['Eigenmode Median (window= ' win_ms 'ms, order= ' order ')'], 'fontsize', 20)
legend(nlxEvHandles,nlxEvLabels,'Location','southwest')
if sum(ismember({events.type},'GO'))>1
    linkaxes([ax2,ax1],'x');
end
if saveFigs
    saveas(gcf,fullfile(dci_dir,['order_' order],'images/medians.png'));
end
%% frequency dampening plot
dt = 1/srate;
freq = zeros(length(AR_mod),size(AR_mod(1).imag,2));
dectime = zeros(length(AR_mod),size(AR_mod(1).imag,2));
for i=1:length(AR_mod)
    theta = angle(AR_mod(i).imag);
    r = abs(AR_mod(i).imag); 
    freq(i,:) = abs(theta)./(2*pi*dt);
    dectime(i,:) = log(r)./dt;
end
dectime2 = [];
freq2 = [];
for i=1:size(dectime,2)
    dectime2 = vertcat(dectime2,dectime(:,i));
    freq2 = vertcat(freq2,freq(:,i));
end
good = find(dectime2>-300);
tempFreq = freq2(good);
figure
histogram(tempFreq, logspace(0, 3, 100));
xlabel('Frequency (Hz)');
ylabel('Number of Eigenmodes');
set(gca,'XScale','log');
figure
histogram(freq2, logspace(0, 4, 800));
xlabel('Frequency (Hz)');
ylabel('Number of Eigenmodes');
set(gca,'XScale','log');

%% p-value plots (t-test from baseline pre-drug window)
baseInd = floor(nlxEvents(3).lfpoffset/win);
for i=1:length(AR_mod)
    [p, h, stats]=ranksum(med(1:baseInd),med(baseInd:end));
end

%% Plot EM distribution
 
% number of AR models calculated
n_tstamps = size(AR_mod,2);
 
% get initial edges for histogram, so that they are all the same
ev = AR_mod(1).ev;
% automatically bins data, and gets count for data points in each bin, with
% specified number of bin
[N, edges] = histcounts(ev,1000);
 
% number of bins by number of timestamps
binned_em = zeros(numel(edges)-1,n_tstamps);
% just matrix of eigen modes
eig_modes = zeros(numel(ev), n_tstamps);
 
for i = 1:n_tstamps-1
    % add to eig_modes
    ev = AR_mod(i).ev;
    eig_modes(1:numel(ev),i) = ev;
    % add to binned em
    tmp = histcounts(ev, edges);
    binned_em(:,i) = tmp;
end
 
% since we're only interested in values around 1, only plot .8 and up
idx = (edges >= .85);
edges = edges(idx);
plot_eig_modes = binned_em(idx(1:end-1),:);
 
% plot
figure
clf
hold on
imagesc(time_s, edges, plot_eig_modes)
colormap('hot')
%caxis([0,5])
colorbar
ylim([min(edges),max(edges)])
xlim([time_s(1), time_s(end)])
 
% plot behavioral COC nd anesthesia
%plot([coc, coc], [1,numel(edges)], 'w', 'linewidth', 3)
%plot([ana_start, ana_start], [1,numel(edges)], 'w')
 
% edges as a string
xlabel('Time Windows (s)', 'fontsize', 15)
ylabel('Criticality Index','fontsize', 15)
title('Histogram of Criticality Over Time','fontsize', 20)
hold off
