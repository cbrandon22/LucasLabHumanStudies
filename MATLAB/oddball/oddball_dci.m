%% Oddball DCI
% This script calculates the criticality at every time point for a given
% data set. Generel workflow is as follows:
%
% Fit an auto regressive model to the data using th ARfit toolbox.
% Find the eigenvalues of this model
% Plot
%
% dynamical_criticality_index.m (JStiso 02/2017) edited by CBrandon 02/2018 for oddball data
% *** see oddball_processing.m and plot_lead_channels.m for prerequisite
% variables (sessInfo, parameters)***

clear
subject_list = {'HUP165_i'};
ddir = '/Volumes/HumanStudies/HumanStudies/oddball/eeg'; %path to folder containing subjects
maxSamplesToLoad = 8000000; % break session into blocks of maxSamplesToLoad for speed
reref_method = 'lead_average'; %'lead_average'
win_ms = 500;
for subi=1:length(subject_list)
    subj = subject_list{subi};
    cond = subj(end); % inductance or emergence
    order = 1; %maximum order to try
    load(fullfile(ddir,subj,'processed/sessInfo.mat'));
    load(fullfile(ddir,subj,'processed/parameters.mat'));
    
    win = round(win_ms*(srate/1000));
    ref_num = floor(5*srate/win); % number of time bins to use as a conscious reference distribution: should equal 5 seconds
    coc = floor(parameters.coc/win);
    ana_start = round(parameters.ana/win);
    
    %load data into electrodes x samples matrix (going to be way too big
    %for high sample rates.) Probably will move this into AR loop to only
    %make matrix of size electrodes x window samples
    good_elecInfo = elecInfo(~ismember(cell2mat(elecInfo(:,1)),bad_channels),:);
    if (win-order)<(size(good_elecInfo,1)*order)
        error('window too small')
    end
    dat = look(events(1).lfpfile,good_elecInfo{1,1},[],1)';
    numSessSamples = length(dat);
    if numSessSamples<maxSamplesToLoad % load full session now
        data = ones(size(good_elecInfo,1),numSessSamples)*NaN;
        for chan = 1:size(good_elecInfo,1)
            dat = look(events(1).lfpfile,good_elecInfo{chan,1},[],1)';
            data(chan,:) = dat;
        end
        sessBlocks = 1;%number of session blocks
        nBlockWin=ceil(numSessSamples/sessBlocks/win);%number of windows per session block (overshoots)
    else
        sessBlocks = ceil(numSessSamples/(maxSamplesToLoad-win));%number of session blocks
        nBlockWin=ceil(numSessSamples/sessBlocks/win);%number of windows per session block (overshoots)
    end
    leads = regexp({good_elecInfo{:,3}},'\D*','match');
    leads = vertcat(leads{:});
    [lead_prefix,~,lead_ind] = unique(leads);
    
    % save window
    save_dir = fullfile('/Volumes/HumanStudies/HumanStudies/oddball/scratch/dci',subj,reref_method,['tw' num2str(win)],'Analysis');
    if ~exist(save_dir),mkdir(save_dir);end
    save(fullfile(save_dir, 'grid_win.mat'), 'win')
    
    pmax = order; 
    pmin = 1;
    
    % initialize data structure
    AR_mod = struct();
    
    % for every 300ms time window, get coefficient matrix
    for sBlock=1:sessBlocks
        fprintf('loading data block %d of %d (%d time windows per block)\n', sBlock,sessBlocks,nBlockWin);
        if sessBlocks~=1
            bStartInd = (sBlock-1)*nBlockWin*win+1;
            bEndInd = sBlock*nBlockWin*win;%last sample in sess block
            if bEndInd>numSessSamples,bEndInd=numSessSamples;end
            data = ones(size(good_elecInfo,1),bEndInd-(bStartInd-1))*NaN;
            for chan = 1:size(good_elecInfo,1)
                dat = look(events(1).lfpfile,good_elecInfo{chan,1},[],1)';
                data(chan,:) = dat(bStartInd:bEndInd); % trim to only include this window (keeps matrix at manageable size)
            end
        else
            bStartInd = 1;
            bEndInd = numSessSamples;
        end
        switch(reref_method)
            case 'common_average'
                disp('Re-referencing to common average');
                data = bsxfun(@minus,data, mean(data,1));
            case 'lead_average'
                disp('Re-referencing to lead averages');
                for l=1:max(lead_ind)
                    data(lead_ind==l,:) = bsxfun(@minus,data(lead_ind==l,:), mean(data(lead_ind==l,:),1));
                end
        end
        disp('Creating AR model')
        for cnt=1:floor(size(data,2)/win)
            winSamples = cnt*win-(win-1):cnt*win;
            % detrend data (avoids interpretation of lf oscillations as exponential growth/decay)
            n_data = detrend(data(:,winSamples)');
            %n_data = n_data - mean(mean(n_data))';
             
            % sbc stands for schwartz criterion
            % w = intercepts (should be 0), A is coefficient matrix, C is noise
            % covariance, th is needed for confidence interval calculation
            % w is returning non zero val, not sure why
            [w,A,C,SBC,FPE,th] = arfit(n_data,pmin,pmax,'sbc'); % data_n should be a time chunk;
            
            % get optimal order
            popt = size(A,2)/size(good_elecInfo,1); % size(A,2) = optimal order * number of channels
            
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
            AR_mod((sBlock-1)*nBlockWin+cnt).ev = abs(ev);
            AR_mod((sBlock-1)*nBlockWin+cnt).imag = ev;
            AR_mod((sBlock-1)*nBlockWin+cnt).evect = s;
            
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
    end
    hw=win/2;
    if ~isempty(parameters.artifact)
        keyboard
        idx1 = (((1:size(AR_mod,2)).*win) >= parameters.artifact(1));
        idx2 = (((1:size(AR_mod,2)).*win - (win-1)) <= parameters.artifact(end));
        art_idx = idx1 & idx2;
        
        % for plotting
        %time_s = ([1:size(AR_mod,2)].*win - hw)./srate;
        time_s = ((1:size(AR_mod,2)).*win - (win/2))./srate;
        time_s = time_s(~art_idx);
        
        % cut out artifact
        AR_mod = AR_mod(~art_idx);
    else
        % for plotting
        %time_s = ((1:size(AR_mod,2)).*win - hw)./srate;
        time_s = ((1:size(AR_mod,2)).*win - (win/2))./srate;
    end
    % get change of consciousness in seconds
    coc = time_s(coc);
    
    % make directories
    if ~exist(fullfile(save_dir,['order_' num2str(popt)], 'images'), 'dir')
        mkdir(fullfile(save_dir,['order_' num2str(popt)], 'images'));
    end
    if ~exist(fullfile(save_dir,['order_' num2str(popt)], 'stats/images/discrete'), 'dir')
        mkdir(fullfile(save_dir,['order_' num2str(popt)], 'stats/images/discrete'));
    end
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