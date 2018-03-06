%% ROC
%
% Fit glm to slow wave, stability, and heart rate data, calculate ROC
% curve
% Currently data was obtained from Eric
%
% @author JStiso March 2017

%% Load data
subjects = [{'HUP125_i'}, {'HUP138_i'}, {'HUP121_i'}, {'HUP119_i'}, {'HUP117_i'}, {'HUP132_i'}];

for s = 1:numel(subjects)
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
    load([data_dir, 'tw', num2str(win), '/Analysis/order_', num2str(order), '/medians']);
    save_dir = [data_dir, 'tw', num2str(win), '/Analysis/'];
    load([data_dir, 'tw', num2str(win), '/Analysis/order_', num2str(order), '/time_s']);
    
    
    % load anethesia start and change of consciousness
    load([data_dir, 'parameters']);
    coc_time = round(parameters.coc/win); % time bin where change in consciousness occurs occurs (arbitrary value right now)
    
    
    %% Make behvioral covariate
    
    oo = false(size(med));
    oo(coc_time:end) = true;
    
    %% Generate Model
    
    mdl = fitglm(med',oo','Distribution','binomial','Link','logit');
    
    %% ROC
    
    scores = mdl.Fitted.Probability;
    [X,Y,T,AUC] = perfcurve(oo',scores,1);
    save([save_dir, '/order_', num2str(order), '/AUC'], 'AUC');
    save([save_dir, '/order_', num2str(order), '/T'], 'T');
    
    clf
    plot(X,Y)
    hold on
    plot(linspace(0,1,numel(X)), linspace(0,1,numel(Y)), 'k')
    xlabel('False positive rate', 'fontsize', 15)
    ylabel('True positive rate', 'fontsize', 15)
    title('ROC for Classification by Logistic Regression', 'fontsize', 20)
    saveas(gca, [save_dir, '/order_', num2str(order), '/images/roc.jpg'], 'jpg')
    
end