%% Damping times
%
% This script plots the damping time vs frequency from AR models calculated
% in dynamical_criticality_index.m. In short, these are the equations used.
% lambda is the magnitude of the complex eigen value, nd delta t is the
% sampling rate
%
% lambda =          p*exp(i*phi);
% lambda =          p(cos(phi) + isin(phi))
% freq =            phi/2pi*delta_t
% damping rate =    log(p)/delta_t
%
% For moreinformation on calculation of
% individual damping times and frequencies, see Solovey et al 2012 (monkey
% paper), of the DCI script. This script is most just for plotting.

% @author JStiso March 2017

%% Loop through subjects

subjects = [{'HUP125_i'}]%, {'HUP133_e'}, {'HUP138_i'}, {'HUP121_i'}, {'HUP119_i'}, {'HUP117_i'}, {'HUP132_i'}, {'HUP121_e'}];

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
    load([data_dir, 'elecs']);
    elecs = elec_info.good;
    load([data_dir, 'tw', num2str(win), '/Analysis/order_', num2str(order), '/AR_mod']);
    save_dir = [data_dir, 'tw', num2str(win), '/Analysis/'];
    
    % other variables
    load([data_dir, 'parameters']);
    ref_num = floor(20*srate/win); % number of time bins to use as a conscious reference distribution: should equal 5 seconds
    
    %% Get data into workable structure
    % namely, sorted vector of damping times, and matched sorting for frequency
    
    % initialize
    dampr_con = zeros(ref_num,numel(elecs));
    freq_con = zeros(size(dampr_con));
    dampr_unc = zeros(ref_num,numel(elecs));
    freq_unc = zeros(size(dampr_unc));
    % move to matrix for conscious and anesthetized
    if strcmp(cond, 'i')
        % conscious
        for i = 1:ref_num
            ev = AR_mod(i).imag;
            %a = real(ev); % real part
            %b = imag(ev); % imaginary part
            p = abs(ev); % magnitude
            phi = angle(ev); % can also be atan2 of a b
            freq_curr = (phi./(2*pi)).*srate;
            dampr_curr = log2(p).*(srate);
            dampr_con(i,:) = dampr_curr;%AR_mod(i).dampr;
            freq_con(i,:) = freq_curr;%AR_mod(i).freq;
        end
        % reshape into vectors
        dampr_v_con = reshape(dampr_con, 1, []);
        freq_v_con = reshape(freq_con, 1, []);
        
        %anesthetized
        cnt = 0;
        for i = numel(AR_mod) - ref_num:numel(AR_mod)
            cnt = cnt + 1;
            ev = AR_mod(i).imag;
            %a = real(ev); % real part
            %b = imag(ev); % imaginary part
            p = abs(ev); % magnitude
            phi = angle(ev); % can also be atan2
            freq_curr = (phi./(2*pi)).*srate;
            dampr_curr = log2(p).*(srate);
            dampr_unc(cnt,:) = dampr_curr;%AR_mod(i).dampr;
            freq_unc(cnt,:) = freq_curr;%AR_mod(i).freq;
        end
        % reshape into vectors
        dampr_v_unc = reshape(dampr_unc, 1, []);
        freq_v_unc = reshape(freq_unc, 1, []);
    else
        % conscious
        for i = numel(AR_mod) - ref_num:numel(AR_mod)
            ev = AR_mod(i).imag;
            %a = real(ev); % real part
            %b = imag(ev); % imaginary part
            p = abs(ev); % magnitude
            phi = angle(ev); % can also be atan2 of a b
            freq_curr = (phi./(2*pi)).*srate;
            dampr_curr = log2(p).*(srate);
            dampr_con(i,:) = dampr_curr;%AR_mod(i).dampr;
            freq_con(i,:) = freq_curr;%AR_mod(i).freq;
        end
        % reshape into vectors
        dampr_v_con = reshape(dampr_con, 1, []);
        freq_v_con = reshape(freq_con, 1, []);
        
        %anesthetized
        cnt = 0;
        for i = 1:ref_num
            cnt = cnt + 1;
            ev = AR_mod(i).imag;
            %a = real(ev); % real part
            %b = imag(ev); % imaginary part
            p = abs(ev); % magnitude
            phi = angle(ev); % can also be atan2
            freq_curr = (phi./(2*pi)).*srate;
            dampr_curr = log2(p).*(srate);
            dampr_unc(cnt,:) = dampr_curr;%AR_mod(i).dampr;
            freq_unc(cnt,:) = freq_curr;%AR_mod(i).freq;
        end
        % reshape into vectors
        dampr_v_unc = reshape(dampr_unc, 1, []);
        freq_v_unc = reshape(freq_unc, 1, []);
    end
    
    
    % select only damped modes (damping rate < 0)
    %idx = dampr_v < 0;
    %dampr_v = dampr_v(idx);
    %freq_v = freq_v(idx);
    
    % get rid of negative freq
    idx_unc = find(freq_v_unc > 0);
    freq_v_unc = freq_v_unc(idx_unc);
    dampr_v_unc = dampr_v_unc(idx_unc);
    idx_con = find(freq_v_con > 0);
    freq_v_con = freq_v_con(idx_con);
    dampr_v_con = dampr_v_con(idx_con);
    % get rid of outliers (eigenmode = 0)
    idx_unc = find(dampr_v_unc > -50000);
    freq_v_unc = freq_v_unc(idx_unc);
    dampr_v_unc = dampr_v_unc(idx_unc);
    idx_con = find(dampr_v_con > -50000);
    freq_v_con = freq_v_con(idx_con);
    dampr_v_con = dampr_v_con(idx_con);
    
    % plot
    clf
    scatter(dampr_v_unc, freq_v_unc, 'filled', 'MarkerFaceAlpha', .2, 'MarkerFaceColor', 'r')    
    hold on
    scatter(dampr_v_con, freq_v_con, 'filled', 'MarkerFaceAlpha', .2, 'MarkerFaceColor', 'b')    
    xlabel('Damping Rate (1/s)')
    ylabel('Frequency (Hz)')
    legend('Unconscious', 'Conscious')
    ylim([0, 1024])
    saveas(gca, [save_dir, '/order_', num2str(order), '/images/damping_rate.jpg'], 'jpg');
    
    % plot
    clf
    scatter(dampr_v_unc, freq_v_unc, 'filled', 'MarkerFaceAlpha', .2, 'MarkerFaceColor', 'r')    
    hold on
    scatter(dampr_v_con, freq_v_con, 'filled', 'MarkerFaceAlpha', .2, 'MarkerFaceColor', 'b')    
    xlabel('Damping Rate (1/s)')
    ylabel('Frequency (Hz)')
    legend('Unconscious', 'Conscious')
    ylim([0, 200])
    saveas(gca, [save_dir, '/order_', num2str(order), '/images/damping_rate_zoom.jpg'], 'jpg');
    
    close all
end
% histogram
%clf
%[N,C] = hist3([freq_v_con', dampr_v_con'], [100 10]);
%imagesc(dampr_v_con, freq_v_con, N)
%colorbar
%colormap('hot')