
function [ config_mm ] = le_mm_config( )
% This function is a config file w various parameters used throughout these
% analyses

config_mm.comparison = 'moveWait';

% selectMotorSites
config_mm.hfa_range = [70 200];
config_mm.hfa_collapseFreqFlag = true;
config_mm.powConfigNum = 1;
config_mm.fdr_p_thresh = 0.05;
config_mm.num_sig_bins_thresh = 1; % this sets how many fdr-corrected sig. time windows you need to select a motor site


% cluster motor sites based on HFA, tfPow and phase analyses
config_mm.fRange = [3 200];
config_mm.tRange = [0 2000];
config_mm.collapseFreqFlag = false;
config_mm.collapseTimeFlag = true;

config_mm.tfplot_flag = 0; % sets whether to plot individual electrodes or not
config_mm.tfplot_printFlag = 1;


config_mm.phase_tRange = [0 2000];
config_mm.phase_collapseTimeFlag = [0 2000];


config_mm.use_kmeans = true;%false; % either use kmeans or explicit feature to sort elecs
if config_mm.use_kmeans 
    config_mm.explicit_feature_to_sort_by = '';
else
%    config_mm.explicit_feature_to_sort_by = 'theta';
    config_mm.explicit_feature_to_sort_by = 'diff(theta,beta)';
end
% for explicit feature sorting
config_mm.featfRange.theta = [3 8];
config_mm.featfRange.beta = [12 30];


config_mm.kmeans_numKtoEval = 4; % num of kmeans clusters to evaluate
config_mm.kmeans_k = 2; % manually coded k in k means (ideally do this based on jumpClus analysis)
config_mm.kmeans_iters = 100; % number of iterations
config_mm.tf_clim_hfa = [-7 7];
config_mm.tf_clim_lf = [-5 3];
config_mm.tf_clim_phaseZ = [-2 2];