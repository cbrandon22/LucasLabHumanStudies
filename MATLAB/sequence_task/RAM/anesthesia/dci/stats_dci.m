%% DCI stats Test
%
% Qualitative: This script plots medians between 'early' and 'late' groups,
% and also looks at the number of modes above and below .98
% quantitative: This script calculates the KS statistic for every
% timepoint, and also calculates a permutation test and nonparametric
% Mann-Whitney U test between the the first time point, and every time
% window following it
%
% All tests save data structures and make plots. Run SUBJ_dci before this
%
% @ author JStiso
%
%% Loop through subjects

subjects = [{'HUP125_i'}]; %[{'HUP133_e'}, {'HUP138_i'}, {'HUP121_i'}, {'HUP119_i'}, {'HUP117_i'}, {'HUP132_i'}, {'HUP121_e'}];

for s = 1:numel(subjects)
    %% Global Varibles
    top_dir = '/Users/tnl/Desktop/C/data/eeg/';
    subj = subjects{s};
    order = 1; % order od model to use
    cond = subj(end); % emergence or inductance
    
    data_dir = [top_dir, subj, '/processed_data/dci/']; % where data is
    % loada window size
    load([data_dir, 'grid_win']);
    load([data_dir, 'srate']); % what the data is sampled at
    ev_dir = [top_dir, subj, '/behavioral/session_0/']; % where behavioral data is
    load([data_dir, 'elecs']);
    elecs = elec_info.good;
    load([data_dir, 'tw', num2str(win), '/Analysis/order_', num2str(order), '/AR_mod']);
    save_dir = [data_dir, 'tw', num2str(win), '/Analysis/'];
    load([data_dir, 'tw', num2str(win), '/Analysis/order_', num2str(order), '/time_s']);
    
    
    % load anethesia start and change of consciousness
    if strcmp(cond, 'i')
        load([data_dir, 'parameters']);
        coc_time = floor(parameters.coc/win); % time bin where change in consciousness occurs occurs (arbitrary value right now)
        coc_time = time_s(coc_time); % time in s where this occurs
    else
        coc_time = time_s(1);
    end
    %ana = parameters.ana/win; % approximate time bin when anasthesia ws delivered/ stopped (arbitrary right now)
    %ana_test = ana; % if you want to pad with some time points, do so here
    
    % stats variables
    n_perm = 1000; % number of permutations
    % if behavioral data isnt working
    ref_num = floor(20*srate/win); % in seconds, and alternative to using change of consciousness
    ana_test = ref_num;
    coc = numel(AR_mod) - ref_num;
    thresh = .98;
    
    %% make directories
    
    if ~exist([data_dir, 'tw', num2str(win), '/Analysis/order_', num2str(order), '/stats/images/continuous'], 'dir')
        mkdir([data_dir, 'tw', num2str(win), '/Analysis/order_', num2str(order), '/stats/images/continuous']);
    end
    
    %% Qualitative: look at medians, and number of modes below thresh
    
    % get medians of conscious and unconcscious state
    % use only data that is definitely conscious, or definitely unconscious
    min_param = min([numel(AR_mod) - coc, ana_test]);
    if min_param == ana_test
        unc = zeros(1, min_param + 1);
        con = zeros(1, min_param + 1);
    else
        unc = zeros(1,numel(AR_mod)- min_param + 1);
        con = zeros(1,numel(AR_mod)- min_param + 1);
    end
    
    % get timeseries of number of data points above the threshold
    thresh_ts = zeros(1,numel(AR_mod));
    for i = 1:numel(AR_mod)
        thresh_ts(i) = sum(AR_mod(i).ev > thresh);
    end
    
    % put into conscious and unconscious bins
    if strcmp(cond, 'i')
        % for getting median
        cnt = 1;
        % unconscious is at the end for induction
        unc_high = thresh_ts((numel(AR_mod) - min_param):numel(AR_mod));
        unc_low = numel(elecs)*order - unc_high;
        % conscious at the beginning
        con_high = thresh_ts(1:numel(con));
        con_low = numel(elecs)*order - con_high;
        for i = (numel(AR_mod) - min_param):numel(AR_mod)
            unc(cnt) = median(AR_mod(i).ev);
            cnt = cnt + 1;
        end
        for j = 1:numel(con)
            con(j) = median(AR_mod(j).ev);
        end
    elseif strcmp(cond, 'e')
        % unconscious is at the beginning for emergence
        unc_high = thresh_ts(1:numel(unc));
        unc_low = numel(elecs)*order - unc_high;
        % conscious at the end
        con_high = thresh_ts((numel(AR_mod) - min_param):numel(AR_mod));
        con_low = numel(elecs)*order - con_high;
        % get medians
        cnt = 1;
        for i = 1:numel(unc)
            unc(i) = median(AR_mod(i).ev);
        end
        for j = (numel(AR_mod) - min_param):numel(AR_mod)
            con(cnt) = median(AR_mod(j).ev);
            cnt = cnt + 1;
        end
    end
    % save threshold data, to plot in excel later (excel is better at bar
    % plots)
    thresh_data = [con_high; con_low; unc_high; unc_low]';
    save([save_dir, 'order_', num2str(order), '/stats/threshold_data'], 'thresh_data', '-v7.3')
    
    figure
    clf
    scatter(1:numel(unc),unc, 'filled', 'MarkerFaceColor', 'red', 'markerfacealpha', .6)
    hold on
    scatter(1:numel(con),con, 'filled', 'MarkerFaceColor', 'blue', 'markerfacealpha', .6)
    legend(gca, 'Unconscious', 'Conscious', 'fontsize', 10)
    title('Distribution of Medians', 'fontsize', 20)
    ylabel('Median Eigenvalue Magnitude', 'fontsize', 15)
    %xlabel('Time, relative to behavioral marker' 'fontsize', 15)
    hold off
    saveas(gca, [save_dir, 'order_', num2str(order), '/images/median_dist.jpg'], 'jpg');
    close
    
    
    % plot  timeseries
    figure('position', [1, 1, 6000 ,400])
    clf
    plot(time_s, thresh_ts)
    title(['Threshold of ', num2str(thresh)], 'fontsize', 20)
    xlabel('Time (s)', 'fontsize', 15)
    ylabel('Number of Eigenvalues Above Threshold', 'fontsize', 15)
    xlim([min(time_s), max(time_s)])
    saveas(gca, [save_dir, 'order_', num2str(order), '/stats/images/continuous/thresh_ts.jpg'], 'jpg');
    
    %% Kolmogorov-Smirnov
    
    % compare ll others to the first eigenmode
    ref_AR_mod = AR_mod(1:ana_test);
    test_dist = zeros(numel(AR_mod(1).ev), numel(ref_AR_mod));
    
    for i = 1:numel(ref_AR_mod)
        test_dist(:,i) = ref_AR_mod(i).ev;
    end
    test_dist_v = reshape(test_dist, 1, []);
    
    % campare all others, save p
    fprintf('\nStarting Kolmogorov-Smirnov...')
    p_ks = zeros(1,numel(AR_mod) - ref_num);
    cnt = 0;
    for i = ref_num:numel(AR_mod);
        cnt = cnt + 1;
        fprintf('\n...%d',i);
        [h,p,ks2stat] = kstest2(test_dist_v,AR_mod(i).ev);
        p_ks(cnt) = p;
    end
    
    % plot
    figure('position', [1, 1, 6000 ,400])
    clf
    plot(time_s(ref_num:end), p_ks)
    hold on
    plot([coc_time coc_time], [max(max(p_ks)) min(min(p_ks))], 'k', 'linewidth', 2);
    %plot([ana-ref_num ana-ref_num], [max(max(p_ks)) min(min(p_ks))], 'k', 'linewidth', 3);
    title('Kolmogorov-Smirnov', 'fontsize', 20)
    xlabel('Time (s)', 'fontsize', 15)
    ylabel('Pval', 'fontsize', 15)
    plot([1,numel(p_ks)], [.05, .05], 'r'); %red line at thresh
    xlim([min(time_s(ref_num:end)), max(time_s(ref_num:end))])
    hold off
    saveas(gca, [save_dir, 'order_', num2str(order), '/stats/images/continuous/ks.jpg'], 'jpg');
    
    
    % save
    save([save_dir, 'order_', num2str(order), '/stats/p_ks'], 'p_ks', '-v7.3')
    
    %% Permutation and Mann-Whitney U tests
    
    % get average reference model, to test all others against
    ref_eig = zeros(numel(AR_mod(1).ev), ref_num);
    %AR_mod = AR_mod(ref_num+1:end);
    
    for i = 1:ref_num
        ref_eig(:,i) = ref_AR_mod(i).ev;
    end
    % vectorize
    ref_eig_v = reshape(ref_eig, 1, []);
    %ref_eig = mean(ref_eig,2);
    
    
    % initialize data structures
    p = zeros(1,numel(AR_mod) - ref_num);
    z = zeros(1,numel(AR_mod) - ref_num);
    
    % for every time point, get P val from permutation between first and
    % current, and mann whitney U test between first and current
    fprintf('\nStrating Mann-Whitney U test and Permutation test for time chunk')
    for i = ref_num:numel(AR_mod)
        curr_mod = AR_mod(i);
        fprintf('\n...%d of %d', i, numel(AR_mod));
        
%         % permutation - takes too long
%         [p_curr, ~, ~] = permtest(ref_eig_v, curr_mod.ev, n_perm);
%         p(i) = p_curr;
        
        % mann-whitney
        [p_mw,h,stats] = ranksum(ref_eig_v, curr_mod.ev);
        z(i) = stats.zval;
    end
    
    % save
    save([save_dir, 'order_', num2str(order), '/stats/p_permtest'], 'p', '-v7.3')
    save([save_dir, 'order_', num2str(order), '/stats/z_mannwhitneyu'], 'z', '-v7.3')
    
    %% plot
    % pval
    figure('position', [1, 1, 6000 ,400])
    clf
    plot(time_s, z)
    hold on
    plot([coc_time coc_time], [max(max(z)) min(min(z))], 'k', 'linewidth', 2);
    %plot([ana-ref_num ana-ref_num], [max(max(z)) min(min(z))], 'k');
    if strcmp(cond, 'i')
        plot([1 numel(z)], [2, 2], 'r')
    else
        plot([1 numel(z)], [-2, -2], 'r')
    end
    title('Mann-Whitney U Test', 'fontsize', 20)
    xlabel('Time (s)', 'fontsize', 15)
    ylabel('zscore', 'fontsize', 15)
    xlim([min(time_s), max(time_s)])
    ylim([min(z) max(z)])
    hold off
    saveas(gca, [save_dir, 'order_', num2str(order), '/stats/images/continuous/mwu.jpg'], 'jpg');
    
%     % perm test plot
%     figure('position', [1, 1, 6000 ,400])
%     clf
%     plot(time_s, p)
%     hold on
%     plot([coc_time coc_time], [max(max(p)) min(min(p))], 'k', 'linewidth', 2);
%     % plot([ana-ref_num ana-ref_num], [max(max(p)) min(min(p))], 'k');
%     title('Permutation Test', 'fontsize', 20)
%     xlabel('Time Chunk', 'fontsize', 15)
%     ylabel('P-value', 'fontsize', 15)
%     plot([1 numel(z)], [.05, .05], 'r')
%     xlim([min(time_s), max(time_s)])
%     hold off
%     saveas(gca, [save_dir, 'order_', num2str(order), '/stats/images/continuous/permtest.jpg'], 'jpg');
    close all
end