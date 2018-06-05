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
subject_list = {'HUP155_i'};
ddir = '/Volumes/HumanStudies/HumanStudies/oddball/eeg'; %path to folder containing subjects
maxSamplesToLoad = 1000000; % break session into blocks of maxSamplesToLoad for speed
reref_method = 'lead_average'; %'lead_average'
win_ms = 500;
decimate_factor = 10; % set to 0 to skip filter
plot_detrend = 0;
plot_coef_mat = 0;
calc_pow_spect = 1;
for subi=1:length(subject_list)
    subj = subject_list{subi};
    cond = subj(end); % inductance or emergence
    order = 1; %maximum order to try
    load(fullfile(ddir,subj,'processed/sessInfo.mat'));
    load(fullfile(ddir,subj,'processed/parameters.mat'));
    if exist(fullfile(ddir,subj,'processed/manual_bad_channels.mat'),'file')==2
        load(fullfile(ddir,subj,'processed/manual_bad_channels.mat'));
    end
    
    if decimate_factor > 0
        orig_win = round(win_ms*(srate/1000));
        win = ceil(win_ms*(srate/1000)/decimate_factor);
        ref_num = ceil((5*srate/win)/decimate_factor); % number of time bins to use as a conscious reference distribution: should equal 5 seconds
        parameters.ana = ceil(parameters.ana/decimate_factor);
        parameters.coc = ceil(parameters.coc/decimate_factor);
        parameters.artifact = ceil(parameters.artifact/decimate_factor);
    else
        win = round(win_ms*(srate/1000));
        ref_num = floor(5*srate/win); % number of time bins to use as a conscious reference distribution: should equal 5 seconds
    end
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
    originalTime = linspace(0, size(dat,2)/srate, size(dat,2));
    if calc_pow_spect
        for chan = 1:size(good_elecInfo,1)
            dat = look(events(1).lfpfile,good_elecInfo{chan,1},[],1)';
            freq = [0.1 10000];
            type='bandpass';
            [pxx,f] = pwelch(dat,10*fix(srate),1*fix(srate),logspace(log10(freq(1)),log10(freq(2)),5000),srate);
            figure('Name',subj,'NumberTitle','off','Units','normalized','Position',[1/4 1/4 1/2 1/2],'Color','w');
            plot(f,pxx);
            axis tight; set(gca,'Box','off','XScale','log','YScale','log');
            xlabel('frequency (Hz)'); ylabel('power (\muV^2/Hz)'); title([subj ' ' good_elecInfo{chan,3}],'Interpreter', 'none');
            saveas(gcf,fullfile('/Volumes/HumanStudies/HumanStudies/oddball/scratch/POWER/figs',[subj '_' num2str(chan) '.png']));
            close
        end
        return
    end
    if decimate_factor > 0;
        filtTime = originalTime(1:decimate_factor:end);
        numSessSamples = ceil(length(dat)/decimate_factor);
    else
        numSessSamples = length(dat);
    end
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
    [lead_prefix,ia,lead_ind] = unique(leads);
    lead_lbls = cell(1,length(leads));
    for l=1:length(lead_prefix)
        lead_lbls(ia(l) + floor(sum(lead_ind==l)/2)) = lead_prefix(l);
    end
    
    % save window
    save_dir = fullfile('/Volumes/HumanStudies/HumanStudies/oddball/scratch/dci',subj,reref_method,num2str(decimate_factor),['tw' num2str(win)],'Analysis');
    if ~exist(save_dir),mkdir(save_dir);end
    save(fullfile(save_dir, 'grid_win.mat'), 'win')
    
    pmax = order; 
    pmin = 1;
    
    % initialize data structure
    AR_mod = struct();
    
    % for every time window, get coefficient matrix
    for sBlock=1:sessBlocks
        fprintf('loading data block %d of %d (%d time windows per block)\n', sBlock,sessBlocks,nBlockWin);
        if sessBlocks~=1
            bStartInd = (sBlock-1)*nBlockWin*win+1; %first sample in sess block
            bEndInd = sBlock*nBlockWin*win;%last sample in sess block
            if bEndInd>numSessSamples,bEndInd=numSessSamples;end
            data = ones(size(good_elecInfo,1),bEndInd-(bStartInd-1))*NaN;
            for chan = 1:size(good_elecInfo,1)
                dat = look(events(1).lfpfile,good_elecInfo{chan,1},[],1)';
                if decimate_factor > 0
                    temp = decimate(dat,decimate_factor); % trim to only include this window (keeps matrix at manageable size)
                    data(chan,:) = temp(bStartInd:bEndInd);
                    figure
                    hold on
                    plot(originalTime,dat)
                    plot(filtTime,temp)
                    set(gca,'XLim',[1500 1502]);
                else
                    data(chan,:) = dat(bStartInd:bEndInd);
                end
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
            if plot_detrend
                win_t = linspace(1,max(winSamples)/srate*1000,max(winSamples));
                for l=1:max(lead_ind)
                    cShift = 1.1*max(max(abs(data(lead_ind==l,winSamples))));
                    cpos = [];
                    figure('Position',[0 0 1920 1080])
                    hold on
                    lead_chan = find(lead_ind==l);
                    for c=1:length(lead_chan)
                        plot(win_t,data(lead_chan(c),winSamples)+((c-1)*cShift),'Color',[0, 0.4470, 0.7410]);
                        plot(win_t,n_data(:,lead_chan(c))+((c-1)*cShift),'Color',[0.8500, 0.3250, 0.0980]);
                        plot([win_t(1), win_t(end)],[((c-1)*cShift), ((c-1)*cShift)],'k','LineWidth',0.5)
                        cpos= [cpos, ((c-1)*cShift)];
                    end
                    legend('raw','detrend')
                    clbls = good_elecInfo(lead_ind==l,3);
                    set(gca,'YTick',cpos,'YTickLabels',clbls,'FontSize',20);
                    title([subj '   Lead: ' lead_prefix{l} '   Window: ' num2str(cnt)],'Interpreter', 'none');
                    xlabel('Time (ms)');
                    ylabel('Channel');
                    keyboard
                end
            end
             
            % sbc stands for schwartz criterion
            % w = intercepts (should be 0), A is coefficient matrix, C is noise
            % covariance, th is needed for confidence interval calculation
            % w is returning non zero val, not sure why
            [w,A,C,SBC,FPE,th] = arfit(n_data,pmin,pmax,'sbc'); % data_n should be a time chunk;
            
            if plot_coef_mat
                figure;
                imagesc(A)
                set(gca,'XTick',[1:length(lead_lbls)],'XTickLabels',lead_lbls,'YTick',[1:length(lead_lbls)],'YTickLabels',lead_lbls,'Box','off','TickLength',[0.001 0.025])
                c = colorbar;
                ylabel(c,'coefficient')
                caxis([-0.35 0.5])
                title([subj '   Time Window ' num2str(cnt) ' coefficient matrix'],'interpreter','none');
            end
            
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
        if decimate_factor > 0 
            time_s = ((1:size(AR_mod,2)).*win - (win/2))./ceil(srate/decimate_factor);
        else
            time_s = ((1:size(AR_mod,2)).*win - (win/2))./srate;
        end
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
    save([save_dir, '/order_', num2str(popt), '/filtTime'], 'filtTime');
    
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