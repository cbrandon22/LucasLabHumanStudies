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
%addpath('/Users/tnl/Desktop/MATLAB/fieldtrip-master/')

%% Load and define global variables

subjects = [{'HUP119_i'}];%[{'HUP133_e'}, {'HUP119_i'}, {'HUP121_i'}, {'HUP117_i'}, {'HUP132_i'}, {'HUP138_i'}, {'HUP121_e'}
% loop through subjects
for s = 1:numel(subjects)
    % define global variables
    top_dir = '/Users/tnl/Desktop/C/data/eeg/';
    subj = subjects{s};
    cond = subj(end); % inductance or emergence
    order = 1; %maximum order to try
    
    data_dir = [top_dir, subj, '/processed_data/dci/']; % where data is
    %data_dir = [top_dir, subj, '/lfp.noreref'];
    srate = load([data_dir, 'srate']); % what the data is sampled at
    srate = srate.srate;
    load([data_dir, 'parameters']) % behavioral info
    load([data_dir, 'elecs']);
    elecs = elec_info.good;
    % window size: N times the number of electrodes
    win = 500*(srate/1000);
    
    % for plotting
    ref_num = floor(5*srate/win); % number of time bins to use as a conscious reference distribution: should equal 5 seconds
    if strcmp(cond,'i')
        coc = floor(parameters.coc/win); % time bin where change in consciousness occurs occurs (arbitrary value right now) 
    else
        coc = 1;
    end
    %ana_start = round(parameters.ana/win);
    
    % load matrix
    % data will be an E x T matrix, where E is the number of electrodes and T
    % is the number of timepoints
    data = load([data_dir, 'data']);
    % for some reson this loads as a struct?
    data = data.data_good;
    
    % save window
    save([data_dir, 'grid_win'], 'win')
    
    %% make directories
    
    if ~exist([data_dir, 'tw', num2str(win), '/Analysis/order_', num2str(order), '/images'], 'dir')
        mkdir([data_dir, 'tw', num2str(win), '/Analysis/order_', num2str(order), '/images']);
    end
    if ~exist([data_dir, 'tw', num2str(win), '/Analysis/order_', num2str(order), '/stats/images/discrete'], 'dir')
        mkdir([data_dir, 'tw', num2str(win), '/Analysis/order_', num2str(order), '/stats/images/discrete']);
    end
    
    save_dir = [data_dir, 'tw', num2str(win), '/Analysis'];
    
    
    %% Fit AR model
    % arfit: optimal model order; size(Aest,2)/m gives optimal parameters
    
    % calculate order range
    % pmax is the time window divided by the electrodes
    % pmax could be floor(win/numel(elecs)) if you want to optimize the order;
    pmax = order; %floor(win/numel(elecs));
    pmin = 1;
    
    % half window size, to be more concise in code
    hw = win/2;
    
    % initialize data structure
    AR_mod = struct();
    
    cnt = 0; % for indexing purposes
    fprintf('\nCreating AR model for time chunk starting at ')
    % for every 300ms time window, get coefficient matrix
    for i = hw+1:win:(size(data,2) - hw);
        fprintf('\n...%d out of %d', i, size(data,2) - hw)
        cnt =  cnt + 1;
        % use: arfit(vector, pmin, pmax, selector)
        
        % detrend and subtract mean from data
        n_data = detrend(data(:,i-hw:i+hw)');
        n_data = n_data - mean(mean(n_data))';
        
        
        % sbc stands for schwartz criterion
        % w = intercepts (should be 0), A is coefficient matrix, C is noise
        % covariance, th is needed for confidence interval calculation
        % w is returning non zero val, not sure why
        [w,A,C,SBC,FPE,th] = arfit(n_data,pmin,pmax,'sbc'); % data_n should be a time chunk;
        %AR_mod(cnt).w = w; AR_mod(cnt).A = A; AR_mod(cnt).C = C; AR_mod(cnt).th = th;
        
        % get optimal order
        popt = size(A,2)/numel(elecs);
        
        % test residuals of model
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
        AR_mod(cnt).imag = ev;
        AR_mod(cnt).evect = s;
        
        % get frequency and samping time
        % lambda = p*exp(i*phi); lambda = p(cos(phi) + isin(phi))
        % freq = phi/2pi*delta_t
        % damping rate = log(p)/delta_t
        %     a = real(ev); % real part
        %     b = imag(ev); % imaginary part
        %     p = sqrt(a.^2 + b.^2); % magnitude
        %     phi = atan2(b,a);
        %     freq = phi./(2*pi*srate/1000);
        %     dampr = log(p)./(srate/1000);
        %     % add to struct
        %     AR_mod(cnt).freq = freq;
        %     AR_mod(cnt).dampr = dampr;
        
        %     [S, Serr, per, tau, exctn, lambda] = armode(A, C, th);
        %     AR_mod(cnt).freq = 1./per;
        %     AR_mod(cnt).ev = lambda;
        %     AR_mod(cnt).evect = S;
        %     AR_mod(cnt).dampr = 1./tau;
    end
    
    if strcmp(subj, 'HUP119_i')
        idx1 = (((1:size(AR_mod,2)).*win - hw) >= parameters.artifact(1));
        idx2 = (((1:size(AR_mod,2)).*win - hw) <= parameters.artifact(end));
        art_idx = idx1 & idx2;
        
        % for plotting
        time_s = ([1:size(AR_mod,2)].*win - hw)./srate;
        time_s = time_s(~art_idx);
        
        % cut out artifact
        AR_mod = AR_mod(~art_idx);
    else
        % for plotting
        time_s = ((1:size(AR_mod,2)).*win - hw)./srate;
    end
    % get change of consciousness in seconds
    coc = time_s(coc);
    
    save([save_dir, '/order_', num2str(popt), '/AR_mod'], 'AR_mod', '-v7.3');
    save([save_dir, '/order_', num2str(popt), '/time_s'], 'time_s', '-v7.3');
    
    %% Plot median over time
    
    med = zeros(1,numel(AR_mod));
    for i = 1:numel(med)
        med(i) = median(AR_mod(i).ev);
    end
    
    clf
    plot(time_s, med)
    hold on
    plot([coc, coc], [min(med), max(med)], 'r', 'linewidth', 3)
    %plot([ana_start, ana_start], [min(med), max(med)], 'k')
    xlim([time_s(1), time_s(end)])
    ylim([min(med), max(med)])
    xlabel('Time (windows)', 'fontsize', 15)
    ylabel('Median', 'fontsize', 15);
    title('Eigenmode Median', 'fontsize', 20)
    
    save([save_dir, '/order_', num2str(popt), '/medians'], 'med', '-v7.3');
    saveas(gca, [save_dir, '/order_', num2str(popt), '/images/medians.jpg'], 'jpg');
    
    %% Plot EM distribution
    
    % not sure what the best way to visualize this is, going to start with 2d
    % histogram
    
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
    figure(1)
    clf
    hold on
    imagesc(time_s, edges, plot_eig_modes)
    colormap('hot')
    %caxis([0,5])
    colorbar
    ylim([min(edges),max(edges)])
    xlim([time_s(1), time_s(end)])
    
    % plot behavioral COC nd anesthesia
    plot([coc, coc], [1,numel(edges)], 'w', 'linewidth', 3)
    %plot([ana_start, ana_start], [1,numel(edges)], 'w')
    
    % edges as a string
    xlabel('Time Windows (s)', 'fontsize', 15)
    ylabel('Criticality Index','fontsize', 15)
    title('Histogram of Criticality Over Time','fontsize', 20)
    hold off
    
    save([save_dir, '/order_', num2str(popt), '/eig_modes'], 'eig_modes', '-v7.3');
    saveas(gca, [save_dir, '/order_', num2str(popt), '/images/ci.jpg'], 'jpg');
    
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
    clf
    bar(x1, n1,'FaceColor','b','EdgeColor','none', 'FaceAlpha', .5);
    hold on
    bar(x2,n2, 'Facecolor', 'r', 'EdgeColor','none', 'FaceAlpha', .5);
    title('Eigenmode Distribution', 'fontsize', 20)
    xlabel('Eigenvalue', 'fontsize', 15)
    ylabel('Number of Modes', 'fontsize', 15)
    legend([{'Conscious'}, {'Unconscious'}], 'fontsize', 10)
    hold off
    saveas(gca, [save_dir, '/order_', num2str(popt), '/stats/images/discrete/group_hist.jpg'], 'jpg');
    
    
    % % get individual time points
    % for i = 1:(size(AR_mod,2)-ref_num);
    %     figure
    %     % change if the reference is taken from the beginning or the end, depending
    %     % on if this is inductance or emergence
    %     if strcmp(cond,'i')
    %         %vectorize
    %         ref_em = reshape(eig_modes(:,1:ref_num), [], 1);
    %         [n1, x1] = hist(ref_em,100);
    %         [n2, x2] = hist(eig_modes(:,i),100);
    %     else
    %         ref_em = reshape(eig_modes(:,end-ref_num:end-1), [], 1);
    %         [n1, x1] = hist(ref_em,100);
    %         [n2, x2] = hist(eig_modes(:,end-i),100);
    %     end
    %     % normalize
    %     n1 = n1/max(n1);
    %     n2 = n2/max(n2);
    %     % back to plotting
    %     bar(x1, n1,'FaceColor','b','EdgeColor','none', 'FaceAlpha', .5);
    %     hold on
    %     bar(x2,n2, 'Facecolor', 'r', 'EdgeColor','none', 'FaceAlpha', .5);
    %     title(['Eigenmode Distribution ', num2str(i)])
    %     xlabel('Eigenvalue')
    %     ylabel('Number of Modes')
    %     legend([{'Conscious'}, {'Unconscious'}])
    %     hold off
    %     saveas(gca, [save_dir, '/order_', num2str(popt), '/stats/images/discrete/', num2str(i)], 'jpg');
    %     close
    % end
    close all
end
