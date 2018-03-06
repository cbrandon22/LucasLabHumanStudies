%% Criticality Index
% This script calculates the criticality at every time point for a given
% data set. Generel workflow is as follows:
%
% Fit an auto regressive model to the data using th ARfit toolbox.
% Find the eigenvalues of this model
% Plot
%
% @author JStiso 02/2017

%% Add paths of toolboxes
addpath(genpath('/Users/tnl/Desktop/MATLAB/arfit/'))
addpath('/Users/tnl/Desktop/MATLAB/fieldtrip-master/')
addpath(genpath('/Users/tnl/Desktop/MATLAB/eeglab14_0_0b/'))

%% Load and define global variables

% define global variables
top_dir = '/Users/tnl/Desktop/C/data/eeg/';
subj = 'HUP119_i';
win = 500; % in ms
cond = subj(end); % inductance or emergence

data_dir = [top_dir, subj, '/processed_data/dci/']; % where data is
%data_dir = [top_dir, subj, '/lfp.noreref'];
srate = 1000; % what the data is sampled at
ev_srate = eegparams('samplerate',[top_dir, subj, '/lfp.noreref', '/params.txt']); %obtain sample rate from params text file
ev_dir = [top_dir, subj, '/behavioral/session_0/']; % where behavioral data is
load([ev_dir, 'events'])
load([data_dir, 'elecs']);
elecs = elec_info.good;

% for plotting
ref_num = 16; % number of time bins to use as a conscious reference distribution

%% Behavioral Data anlysis
% This will calculate the time window we wish to use, as well as the
% rection time percentage to data, which will be compared to the
% criticality index result

% find correct, incorrect, and pass (LOC) recall trials. Result is a subset
% of events
[correctR, incorrectR] = eventStrategize(events);


% change of consciousness, as indexed by the last correct trial
coc = (correctR(size(correctR,2)).lfpoffset/ev_srate) * 1000;% convert from s to ms
% convert to plot points, within which window did this occur
coc = floor(coc/win);
% start of anasthesia (for HUP119)
load(['/Users/tnl/Dropbox/NEV/', subj, 'NEV.mat']);
% start time - first TTL pulse
%ana_start = 810058533-215667820; % determined by looking through NEV file
%ana_start = 158612032 - 125250496; % HUP121i
%ana_start = 2.085477936000000e9 - 1.927974984000000e9; % HUP117, using propofol increase
%ana_start = 903788448 - 478944116; %HUP 108
ana_start = 302486416 - 45939498; % HUP 60
ana_start = (ana_start/ev_srate); % to ms
ana_start = ceil(ana_start/win); % convert to time windows

% load matrix
% data will be an E x T matrix, where E is the number of electrodes and T
% is the number of timepoints
data = load([data_dir, 'data']);
% for some reson this loads as a struct?
data = data.data_good;

% save parameters for later
parameters.ana = ana_start;
parameters.coc = coc;
save([data_dir, 'Analysis/parameters'], 'parameters')

%% Fit AR model
% arfit: optimal model order; size(Aest,2)/m gives optimal parameters

% calculate order range
% pmax is the time window divided by the electrodes
% pmax could be floor(win/numel(elecs)) if you want to optimize the order;
pmax = 1; %floor(win/numel(elecs));
pmin = 1;

% half window size, to be more concise in code
hw = win/2;

% initialize data structure
AR_mod = struct();

cnt = 0; % for indexing purposes
fprintf('\nCreating AR model for time chunk starting at ')
% for every 300ms time window, get coefficient matrix
for i = hw+1:win:( size(data,2) - hw);
    fprintf('\n...%d out of %d', i, size(data,2) - hw)
    cnt=  cnt + 1;
    % use: arfit(vector, pmin, pmax, selector)
    
    % detrend and subtract mean from data
    n_data = detrend(data(:,i-hw:i+hw)');
    n_data = n_data - mean(mean(n_data))';
    
    
    % sbc stands for schwartz criterion
    % w = intercepts (should be 0), A is coefficient matrix, C is noise
    % covariance, th is needed for confidence interval calculation
    % w is returning non zero val, not sure why
    [w,A,C,SBC,FPE,th] = arfit(n_data,pmin,pmax,'sbc'); % data_n should be a time chunk;
    AR_mod(cnt).w = w; AR_mod(cnt).A = A; AR_mod(cnt).C = C; AR_mod(cnt).th = th;
    
    % get optimal order
    popt = size(A,2)/numel(elecs);
    
    % test residuls of model
    % arres: tests significance of residuals
    % acf: plots autocorrelations, shouls lie within dashed confidence interval
    
    % there is always autocorrelation, commenting this out to save time
    % siglev is significance level, res is a time series of residuals
    %     [siglev,res] = arres(w,A,n_data);
    %     % check for autocorrelation
    %     AR_mod(cnt).siglev = siglev;
    %     if siglev > .05
    %         AR_mod(cnt).ac = 0;
    %     else
    %         AR_mod(cnt).ac = 1;
    %     end
    %     % also save confidence intervals
    %     [Aerr, werr] = arconf(A, C, w, th);
    %     AR_mod(cnt).Aci = Aerr;
    
    % calculate and store eigenvalues
    % armode: computes eigen decomposition, tau gives damping times
    % max A square if mode is greater than 1
    if pmax > 1
        r = size(A,1);
        c = size(A,2);
        [s, ev] = eig([A; eye(c-r) zeros(c-r,r)]);
        % get true eigen values
        ev = diag(ev)';
        % adjust phase of eigenmodes
        s = adjph(s);
    else
        [s, ev]  = eig(A);
        % get true eigen values
        ev = diag(ev)';
        % adjust phase of eigenmodes
        s = adjph(s);
    end
    % take the absolute value to get criticality index
    AR_mod(cnt).ev = abs(ev);
    AR_mod(cnt).evect = s;
end

save([data_dir, 'Analysis/order_', num2str(popt), '/AR_mod_tw', num2str(win)], 'AR_mod', '-v7.3');

%% Plot median over time

med = zeros(1,numel(AR_mod));
for i = 1:numel(med)
    med(i) = median(AR_mod(i).ev);
end

plot(med)
xlabel('Time (bins)')
ylabel('Median');
title('Eigenmode Median')

save([data_dir, 'Analysis/order_', num2str(popt), '/medians'], 'med', '-v7.3');
saveas(gca, [data_dir, 'Analysis/order_', num2str(popt), '/images/medians.jpg'], 'jpg');

%% Plot EM distribution
% not sure what the best way to visualize this is, going to start with 2d
% histogram

% number of AR models calculated
n_tstamps = size(AR_mod,2);

% get initial edges for histogram, so that they are all the same
ev = AR_mod(1).ev;
% automatically bins data, and gets count for data points in each bin, with
% specified number of bin
[N, edges] = histcounts(ev,300);

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
colormap('bone')

% plot
figure(1)
hold on
imagesc(plot_eig_modes)
%caxis([0,5])
colorbar
ylim([1,numel(edges)])
xlim([1,size(plot_eig_modes,2)])

% plot behavioral COC nd anesthesia
plot([coc, coc], [1,numel(edges)], 'w')
plot([ana_start, ana_start], [1,numel(edges)], 'w')

set(gca, 'ytick', 1:numel(edges))
% edges as a string
set(gca, 'yticklabels',(strtrim(cellstr(num2str(edges'))')))
xlabel(['Time Windows (', num2str(win), ' ms each)'])
ylabel('Criticality Index')
hold off

save([data_dir, 'Analysis/order_', num2str(popt), '/eig_modes'], 'eig_modes', '-v7.3');
saveas(gca, [data_dir, 'Analysis/order_', num2str(popt), '/images/ci.jpg'], 'jpg');

%% Plot Reference Distribution against Individual Time Points

% get conscious vs unconscious distribution
if strcmp(cond,'i')
    %vectorize
    ref_em1 = reshape(eig_modes(:,1:ref_num), [], 1);
    [n1, x1] = hist(ref_em1,100);
    ref_em2 = reshape(eig_modes(:,end-ref_num:end-1), [], 1);
    [n2, x2] = hist(ref_em2,100);
    
else
    ref_em1 = reshape(eig_modes(:,end-ref_num:end-1), [], 1);
    [n1, x1] = hist(ref_em1,100);
    ref_em2 = reshape(eig_modes(:,1:ref_num), [], 1);
    [n2, x2] = hist(ref_em2,100);
end
% normalize
n1 = n1/max(n1);
n2 = n2/max(n2);
% back to plotting
bar(x1, n1,'FaceColor','b','EdgeColor','none', 'FaceAlpha', .5);
hold on
bar(x2,n2, 'Facecolor', 'r', 'EdgeColor','none', 'FaceAlpha', .5);
title('Eigenmode Distribution')
xlabel('Eigenvalue')
ylabel('Number of Modes')
legend([{'Conscious'}, {'Unconscious'}])
hold off
saveas(gca, [data_dir, 'Analysis/order_', num2str(popt), '/stats/images/discrete/group_hist.jpg'], 'jpg');
close

% get individual time points
for i = 1:(size(AR_mod,2)-ref_num);
    figure
    % change if the reference is taken from the beginning or the end, depending
    % on if this is inductance or emergence
    if strcmp(cond,'i')
        %vectorize
        ref_em = reshape(eig_modes(:,1:ref_num), [], 1);
        [n1, x1] = hist(ref_em,100);
        [n2, x2] = hist(eig_modes(:,i),100);       
    else
        ref_em = reshape(eig_modes(:,end-ref_num:end-1), [], 1);
        [n1, x1] = hist(ref_em,100);
        [n2, x2] = hist(eig_modes(:,end-i),100);
    end
    % normalize
    n1 = n1/max(n1);
    n2 = n2/max(n2);
    % back to plotting
    bar(x1, n1,'FaceColor','b','EdgeColor','none', 'FaceAlpha', .5);
    hold on
    bar(x2,n2, 'Facecolor', 'r', 'EdgeColor','none', 'FaceAlpha', .5);
    title(['Eigenmode Distribution ', num2str(i)])
    xlabel('Eigenvalue')
    ylabel('Number of Modes')
    legend([{'Conscious'}, {'Unconscious'}])
    hold off
    saveas(gca, [data_dir, 'Analysis/order_', num2str(popt), '/stats/images/discrete/', num2str(i)], 'jpg');
    close
end

