%% Oddball DCI
% This script calculates the criticality at every time point for a given
% data set. Generel workflow is as follows:
%
% Fit an auto regressive model to the data using th ARfit toolbox.
% Find the eigenvalues of this model
% Plot
%
% @author JStiso 02/2017
% dynamical_criticality_index.m edited by CBrandon 02/2018 for oddball data
% *** see oddball_processing.m and plot_lead_channels.m for prerequisite
% variables (sessInfo, parameters)***

clear
subject_list = {'HUP159_e'};
ddir = '/Volumes/HumanStudies/HumanStudies/oddball/eeg'; %path to folder containing subjects
for subi=1:length(subject_list)
    subj = subject_list{subi};
    cond = subj(end); % inductance or emergence
    order = 1; %maximum order to try
    load(fullfile(ddir,subj,'processed/sessInfo.mat'));
    load(fullfile(ddir,subj,'processed/parameters.mat'));
    
    win = round(500*(srate/1000));
    ref_num = floor(5*srate/win); % number of time bins to use as a conscious reference distribution: should equal 5 seconds
    coc = floor(parameters.coc/win);
    ana_start = round(parameters.ana/win);
    
    %load data into electrodes x samples matrix (going to be way too big
    %for high sample rates.) Probably will move this into AR loop to only
    %make matrix of size electrodes x window samples
%     files = dir(fullfile(ddir,subj,'processed'));
%     goodChan = 0;
%     data_elecInfo = {};
%     data = [];
%     for i=1:length(files)
%         if ~ismember(files(i).name,{'.','..','.DS_Store','sessInfo.mat'})
%             load(fullfile(ddir,subj,'processed',files(i).name),'lead_elecInfo');
%             includeChan = ~ismember(cell2mat(lead_elecInfo(:,1)),bad_channels);
%             leadChannels = cell2mat(lead_elecInfo(:,1));
%             leadChannels = leadChannels(includeChan);
%             leadEEG = [];
%             for j=1:length(leadChannels)
%                 leadEEG(j,:) = look(events(1).lfpfile,leadChannels(j),[],1)';
%             end
%             data = [data;leadEEG];
%             goodChan = goodChan+sum(includeChan);
%             data_elecInfo = [data_elecInfo;lead_elecInfo(includeChan,:)];
%         end
%     end
    good_elecInfo = elecInfo(~ismember(cell2mat(elecInfo(:,1)),bad_channels),:);
    dat = look(events(1).lfpfile,good_elecInfo{1,1},[],1)';
    numSessSamples = length(dat);
    pmax = order; 
    pmin = 1;
    
    % initialize data structure
    AR_mod = struct();
    
    disp('Creating AR model')
    % for every 300ms time window, get coefficient matrix
    for cnt=1:floor(numSessSamples/win)
        winSamples = cnt*win-(win-1):cnt*win;
        fprintf('time window %d of %d', cnt, floor(numSessSamples/win))
        data = ones(size(good_elecInfo,1),win)*NaN;
        for chan = 1:size(good_elecInfo,1)
            dat = look(events(1).lfpfile,good_elecInfo{chan,1},[],1)';
            data(chan,:) = dat(winSamples); % trim to only include this window (keeps matrix at manageable size)
        end
        % detrend and subtract mean from data
        n_data = detrend(data');
        n_data = n_data - mean(mean(n_data))';
        
        
        % sbc stands for schwartz criterion
        % w = intercepts (should be 0), A is coefficient matrix, C is noise
        % covariance, th is needed for confidence interval calculation
        % w is returning non zero val, not sure why
        [w,A,C,SBC,FPE,th] = arfit(n_data,pmin,pmax,'sbc'); % data_n should be a time chunk;
        
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
end