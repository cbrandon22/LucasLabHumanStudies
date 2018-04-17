function [tMat_pow_cat,tStruct_lf,anatAxis,config_pow,insigCentroids,sigCentroids,insigAnat] = le_mm_cat_tfPowPhase(subj,eLbl_list,printFolderLabel,comparison,indvPlot,savePowSubset,clus)
% This function computes power and phase analyses comparing moveWait using
% parameters set in the config file

% Inputs
% subj  =        'HUP001'
% elecLbl_list   cell array of elecLbls (indicates motor sites); output of
% le_mm_selectMotorSites

% Outputs
% tfMat_zPow    z-scored power differences between move and wait
% tfMat_tstats  tstats comparing the move and wait epochs
% tfMat_pvals   p values for various TF bins

% load dirs and config
dirs = le_dirs('motormap');
config_mm = le_mm_config;

% parse inputs
if ~exist('printFolderLabel','var') || isempty(printFolderLabel)
    printFolderLabel = '';
end

if ~exist('comparison','var') || isempty(comparison)
    comparison = config_mm.comparison;
end

if ~exist('savePowSubset','var') || isempty(savePowSubset)
    savePowSubset = false;
end

if ~exist('savePowSubset','var') || isempty(savePowSubset)
    clus = 'noLbl';
end

% load tStructs (low frequency, TF plots)
[tStruct_lf,config_pow] = le_calcTs_wrapper(subj,[],'motormap',[],comparison,...
    config_mm.fRange,config_mm.collapseFreqFlag,config_mm.powConfigNum,eLbl_list);

% load low frequency phase reset (TF plots)
% [phaseStruct] = le_calcPhaseReset_wrapper(subj,[],'motormap',config_mm.comparison,...
%     config_mm.fRange,config_mm.tRange,config_mm.collapseFreqFlag,...
%     config_mm.collapseTimeFlag,config_mm.powConfigNum,eLbl_list);

% make individual TF plots
if config_mm.tfplot_flag
    %le_plotTs(tStruct_hfa,config_pow,config_mm.tfplot_printFlag,printFolderLabel);
    le_plotTs(tStruct_lf,config_pow,config_mm.tfplot_printFlag,printFolderLabel);
end
if savePowSubset
    save(fullfile(dirs.scratch,'POWER',[subj '_' clus '.mat']),'tStruct_lf');
    tMat_pow_cat=[];
    anatAxis=[];
    insigCentroids=[];
    sigCentroids=[];
    insigAnat=[];
    return
end
% concatonate LF TMats for plotting and dimensionality reduction
tMat_pow_cat = squeeze((nanmean(cat(3,tStruct_lf.tMat),2)))';
%zMat_phase_cat = squeeze((nanmean(cat(3,phaseStruct.zMat),2)))';


% sort all cat matrixes by some feature 
if indvPlot
    if config_mm.use_kmeans && length(eLbl_list)>1
        tMat_cat =  tMat_pow_cat; %[tMat_pow_cat zMat_phase_cat];
        
        % do k-means clustering on this matrix to sort these into
        % groups (clusId describes group membership)
        [clusId] = kmeans(tMat_cat,config_mm.kmeans_k,'replicates',config_mm.kmeans_iters);
        
        
        % evaluate num of clusters
        %     [jumps,SUMD] = JumpClus(tMat_cat,config_mm.kmeans_numKtoEval);
        %     swagFig;
        %     plot(1:config_mm.kmeans_numKtoEval,SUMD,'-k','linewidth',3);
        %     swagAxes(gca,20,'Number of clusters','Model error (summed distance)')
        [~,sortIdx] = sort(clusId);
        
    else % sort on some explicit feature
        [sortIdx,exp_feat] = sortByExplicit_local(tMat_pow_cat,...
            tStruct_lf,config_mm,config_pow);
    end
    
    % sort mats
    tMat_pow_cat = tMat_pow_cat(sortIdx,:);
    %zMat_phase_cat = zMat_phase_cat(sortIdx,:);
    eLbl_list = eLbl_list(:,sortIdx);
    % plot matrixes
    
    % get values for x-axis (time bins) and yaxis (freq bins, accounting for ...
    %custom fRange used to generate this tStruct)
    tBins = nanmean(config_pow.timeBins,2)';
    fBins = [config_pow.freqBins];
    fBins = fBins(tStruct_lf(1).fInd); % filter by freq range used for this tStruct
end

 if config_mm.use_kmeans
     tit =    ['k-means k = ' num2str(config_mm.kmeans_k)];
 else
     tit = ['sorted based on ' config_mm.explicit_feature_to_sort_by];
 end

% Use anatomical labels for y axis instead of electrode pairs
[anatStruct] = le_centroids2Anatstruct(subj);
allELbls = {anatStruct.eLbl};
anatAxis = cell(1,length(eLbl_list));
insigAnat = {};
sigCentroids = {};
insigCentroids = {};
for i=1:length(eLbl_list)
    [~,Locb] = ismember(eLbl_list(i),allELbls);
    anatAxis(i) = {anatStruct(Locb).anatAbbr};
    sigAnatNum(i) = {anatStruct(Locb).anatNum};
    sigCentroids{i,1}=subj;
    sigCentroids{i,2} = anatStruct(Locb).X;
    sigCentroids{i,3} = anatStruct(Locb).Y;
    sigCentroids{i,4} = anatStruct(Locb).Z;
    %sigCentroids{i,5} = anatStruct(Locb).anatNum;
end
j=1;
for i=1:length(allELbls)
    if ~ismember(allELbls(i),eLbl_list);
        insigAnat = [insigAnat,{anatStruct(i).anatAbbr}];
        insigCentroids{j,1}=subj;
        insigCentroids{j,2} = anatStruct(i).X;
        insigCentroids{j,3} = anatStruct(i).Y;
        insigCentroids{j,4} = anatStruct(i).Z;
        insigCentroids{j,5} = 0;
        %insigCentroids{j,5} = anatStruct(i).anatNum;
        j=j+1;
    end
end
% allCentroidsAnat = [sigCentroids;insigCentroids];
if indvPlot
    % TMat Power
    swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
    imagesc(1:length(fBins),1:length(eLbl_list),tMat_pow_cat);
    set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
    colormap jet; colorbar;
    swagAxes(gca,20,'Frequency (Hz)','Electrodes',{'t(LF MOVE - LF WAIT)';tit});
    yt = [1:1:length(eLbl_list)];
    set(gca,'ytick',yt,'yticklabel',eLbl_list(yt),'yticklabelrotation',0)
    setFreqYTicks(gca,fBins)
end

% ZMat Phase Reset
% swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
% imagesc(1:length(fBins),1:length(eLbl_list),real(zMat_phase_cat)); 
% set(gca,'ydir','normal','clim',config_mm.tf_clim_phaseZ)
% colormap jet; colorbar;
% swagAxes(gca,20,'Frequency (Hz)','Electrodes',{'z(LF MOVE - LF WAIT)';tit});
% yt = [1:1:length(eLbl_list)];
% set(gca,'ytick',yt,'yticklabel',anatAxis(yt),'yticklabelrotation',0)
% setFreqYTicks(gca,fBins)


%keyboard
function [] = setFreqYTicks(ax,fBins)
yt = [[1:5:length(fBins)] length(fBins)];
if any(diff(yt)==0)
    yt = [1:5:length(fBins)];
end
set(ax,'xtick',yt,...
'xticklabel',round(fBins(yt)))


function [sortIdx,exp_feat] = sortByExplicit_local(tMat_lf_cat,...
    tStruct_lf,config_mm,config_pow)
fBins = [config_pow.freqBins];
fBins = fBins(tStruct_lf(1).fInd); % filter by freq range used for this tStruct


switch config_mm.explicit_feature_to_sort_by
    case {'theta','beta'}
        if strcmp(config_mm.explicit_feature_to_sort_by,'theta')
            fRange = config_mm.featfRange.theta;
        elseif strcmp(config_mm.explicit_feature_to_sort_by,'beta')
            fRange = config_mm.featfRange.theta;
        end
        [~,fInd_start] = min(abs(fRange(1) - fBins));
        [~,fInd_end] = min(abs(fRange(2) - fBins));
        fInd = fInd_start:fInd_end;

        % collapse based on that feature
        exp_feat = nanmean(tMat_lf_cat(:,fInd),2);
        [~,sortIdx] = sort(exp_feat);

case {'diff(theta,beta)'}
        %collapse theta pow
        fRange = config_mm.featfRange.theta;
        [~,fInd_start] = min(abs(fRange(1) - fBins));
        [~,fInd_end] = min(abs(fRange(2) - fBins));
        fInd = fInd_start:fInd_end;
        theta_pow = nanmean(tMat_lf_cat(:,fInd),2);
        
        %collapse beta pow
        fRange = config_mm.featfRange.beta;
        [~,fInd_start] = min(abs(fRange(1) - fBins));
        [~,fInd_end] = min(abs(fRange(2) - fBins));
        fInd = fInd_start:fInd_end;
        beta_pow = nanmean(tMat_lf_cat(:,fInd),2);
        
        exp_feat = [theta_pow beta_pow];
        [~,sortIdx] = sort(exp_feat(:,1) - exp_feat(:,2));

end


% GRAVEYARD 1 - Contains HFA time series mat and Phase struct
% load HFA tStructs for motor sites)
%[tStruct_hfa] = le_calcTs_wrapper(subj,[],'motormap',config_mm.comparison,...
%    config_mm.hfa_range,config_mm.hfa_collapseFreqFlag,config_mm.powConfigNum,eLbl_list);

% load low frequency phase reset (TF plots)
%[phaseStruct] = le_calcPhaseReset_wrapper(subj,[],'motormap',config_mm.comparison,...
%    config_mm.lf_range,config_mm.lf_collapseFreqFlag,config_mm.powConfigNum,eLbl_list);

%tMat_hfa_cat = cat(1,tStruct_hfa.tMat);
%phaseMat_lf_cat = squeeze((nanmean(cat(3,phaseStruct.rbarMat1)-cat(3,phaseStruct.rbarMat2),2)))';

% this matrix contains two pieces of data: 1) a time resolved HFA matrix of
% tstats to capture temporal dynamics of the signal, 
% and 2) a frequency-resolved low frequency matrix 
%tMat_cat = [tMat_hfa_cat tMat_lf_cat]; 

% sort based on sortIdx (feature-based or multi-variate)
%tMat_hfa_cat = tMat_hfa_cat(sortIdx,:);
%phaseMat_lf_cat = phaseMat_lf_cat(sortIdx,:);

% % HFA mat
% swagFig([0.0188    0.0663    0.4297    0.8150]); hold all
% imagesc(tBins,1:length(eLbl_list),tMat_hfa_cat); 
% set(gca,'ydir','normal','clim',config_mm.tf_clim_hfa)
% colormap jet; colorbar
% swagAxes(gca,20,'Time from GO cue (ms)','Electrodes',{'t(HFA MOVE - HFA WAIT)';tit});
% yt = [1:1:length(eLbl_list)];
% set(gca,'ytick',yt,'yticklabel',eLbl_list(yt),'yticklabelrotation',0)
% % label plot and set preferences
% plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [1 1 1])
% 
% % LF Phase mat
% swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
% imagesc(1:length(fBins),1:length(eLbl_list),phaseMat_lf_cat); 
% set(gca,'ydir','normal','clim',[-.15 .15])
% colormap jet; colorbar;
% swagAxes(gca,20,'Frequency (Hz)','Electrodes',{'Phase Non-Uniformity(LF MOVE - LF WAIT)';tit});
% yt = [1:1:length(eLbl_list)];
% set(gca,'ytick',yt,'yticklabel',eLbl_list(yt),'yticklabelrotation',0)
% setFreqYTicks(gca,fBins)
