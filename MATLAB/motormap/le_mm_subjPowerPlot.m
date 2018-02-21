% Plot subject's power spectrum by cluster

subjList = {'HUP084','HUP087','HUP091','HUP092','HUP111'};
plotInsig = false;
plotIndividualClusters = false;
plotAllBipols = true;
plotSigTogether = false;
savePlots = true;

dirs = le_dirs;
config_mm = le_mm_config;
[eInfoStruct] = le_mm_readElecInfo;

for i=1:length(subjList)
    clear subjAnat;
    subj = subjList{i};
    plotDir = fullfile(dirs.scratch,'figs','motormap',subj,'all','moveWait-motor','clustered');
    % Build eLbl_list for each cluster separately
    eLbl_list_clus1 = {};
    eLbl_list_clus2 = {};
    eLbl_list_clus0 = {};
    for j=1:length(eInfoStruct)
        if strcmp(eInfoStruct(j).subj,subj) && eInfoStruct(j).cluster == 1
            [eLbl_list_clus1] = [eLbl_list_clus1, strcat(num2str(eInfoStruct(j).eLbl1),'-',num2str(eInfoStruct(j).eLbl2))];
        elseif strcmp(eInfoStruct(j).subj,subj) && eInfoStruct(j).cluster == 2
            [eLbl_list_clus2] = [eLbl_list_clus2, strcat(num2str(eInfoStruct(j).eLbl1),'-',num2str(eInfoStruct(j).eLbl2))];
        elseif strcmp(eInfoStruct(j).subj,subj) && eInfoStruct(j).cluster == 0
            [eLbl_list_clus0] = [eLbl_list_clus0, strcat(num2str(eInfoStruct(j).eLbl1),'-',num2str(eInfoStruct(j).eLbl2))];
        end
    end
    
    % Get power t stats for each cluster, make plots
    % Cluster 1
    if ~isempty(eLbl_list_clus1)
        [clus_tPow,tStruct,clusAnat,config_pow] = le_mm_cat_tfPowPhase(subj,eLbl_list_clus1,'motor');
        if ~exist('subjAnat')
            tStruct_lf = tStruct;
            subjAnat = clusAnat;
            tMat_pow_cat = clus_tPow;
        else
            tStruct_lf = [tStruct_lf,tStruct];
            subjAnat = [subjAnat,clusAnat];
            tMat_pow_cat = cat(1,tMat_pow_cat,clus_tPow);
        end
        % Plot cluster 1
        if plotIndividualClusters
            tBins = nanmean(config_pow.timeBins,2)';
            fBins = [config_pow.freqBins];
            fBins = fBins(tStruct_lf(1).fInd);
            
            swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
            imagesc(1:length(fBins),1:length(clusAnat),clus_tPow); 
            set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
            colormap jet; colorbar;
            tit = strcat(subj,' Cluster 1 Bipolars');
            swagAxes(gca,10,'Frequency (Hz)','Bipolar Location',{'t(LF Move - LF WAIT)';tit});
            yt = [1:1:length(clusAnat)];
            set(gca,'ytick',yt,'yticklabel',clusAnat(yt),'yticklabelrotation',0)
            le_setFreqYTicks(gca,fBins)
            
            % Save this plot
            if savePlots
               if ~exist(plotDir,'dir'),mkdir(plotDir);end
               cd(plotDir);
               saveas(gca,'cluster1','png');
            end
        end
    end
    
    % Cluster 2
    if ~isempty(eLbl_list_clus2)
        [clus_tPow,tStruct,clusAnat,config_pow] = le_mm_cat_tfPowPhase(subj,eLbl_list_clus2,'motor');
        if ~exist('subjAnat')
            tStruct_lf = tStruct;
            subjAnat = clusAnat;
            tMat_pow_cat = clus_tPow;
        else
            tStruct_lf = [tStruct_lf,tStruct];
            subjAnat = [subjAnat,clusAnat];
            tMat_pow_cat = cat(1,tMat_pow_cat,clus_tPow);
        end
        % Plot cluster 2
        if plotIndividualClusters
            tBins = nanmean(config_pow.timeBins,2)';
            fBins = [config_pow.freqBins];
            fBins = fBins(tStruct_lf(1).fInd);
            
            swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
            imagesc(1:length(fBins),1:length(clusAnat),clus_tPow); 
            set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
            colormap jet; colorbar;
            tit = strcat(subj,' Cluster 2 Bipolars');
            swagAxes(gca,10,'Frequency (Hz)','Bipolar Location',{'t(LF Move - LF WAIT)';tit});
            yt = [1:1:length(clusAnat)];
            set(gca,'ytick',yt,'yticklabel',clusAnat(yt),'yticklabelrotation',0)
            le_setFreqYTicks(gca,fBins)
            
            % Save this plot
            if savePlots
               if ~exist(plotDir,'dir'),mkdir(plotDir);end
               cd(plotDir);
               saveas(gca,'cluster2','png');
            end
        end
    end
    
    % Plot Clusters 1 and 2 together
    if plotSigTogether
        tBins = nanmean(config_pow.timeBins,2)';
        fBins = [config_pow.freqBins];
        fBins = fBins(tStruct_lf(1).fInd);
            
        swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
        imagesc(1:length(fBins),1:length(subjAnat),tMat_pow_cat); 
        set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
        colormap jet; colorbar;
        tit = strcat(subj,' Significant Bipolars');
        swagAxes(gca,10,'Frequency (Hz)','Bipolars',{'t(LF Move - LF WAIT)';tit});
        yt = [1:1:length(subjAnat)];
        set(gca,'ytick',yt,'yticklabel',subjAnat(yt),'yticklabelrotation',0)
        le_setFreqYTicks(gca,fBins)
        
        % Save this plot
        if savePlots
           if ~exist(plotDir,'dir'),mkdir(plotDir);end
           cd(plotDir);
           saveas(gca,'significantBPs','png');
        end
    end
    
    % Insignificant bipolars
    if ~isempty(eLbl_list_clus0)
        [clus_tPow,tStruct,clusAnat,config_pow] = le_mm_cat_tfPowPhase(subj,eLbl_list_clus0,'motor');
        if ~exist('subjAnat')
            tStruct_lf = tStruct;
            subjAnat = clusAnat;
            tMat_pow_cat = clus_tPow;
        else
            tStruct_lf = [tStruct_lf,tStruct];
            subjAnat = [subjAnat,clusAnat];
            tMat_pow_cat = cat(1,tMat_pow_cat,clus_tPow);
        end
        % Plot Insignificant bipolars
        if plotIndividualClusters
            tBins = nanmean(config_pow.timeBins,2)';
            fBins = [config_pow.freqBins];
            fBins = fBins(tStruct_lf(1).fInd);
            
            swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
            imagesc(1:length(fBins),1:length(clusAnat),clus_tPow); 
            set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
            colormap jet; colorbar;
            tit = strcat(subj,' Insignificant Bipolars');
            swagAxes(gca,10,'Frequency (Hz)','Electrodes',{'t(LF Move - LF WAIT)';tit});
            yt = [1:1:length(clusAnat)];
            set(gca,'ytick',yt,'yticklabel',clusAnat(yt),'yticklabelrotation',0)
            le_setFreqYTicks(gca,fBins)
            
            % Save this plot
            if savePlots
               if ~exist(plotDir,'dir'),mkdir(plotDir);end
               cd(plotDir);
               saveas(gca,'insignificantBPs','png');
            end
        end
    end
    
    % Plot all electrodes for patient (from top to bottom of plot:cluster 0,1,2)
    if plotAllBipols
        tBins = nanmean(config_pow.timeBins,2)';
        fBins = [config_pow.freqBins];
        fBins = fBins(tStruct_lf(1).fInd);
            
        swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
        imagesc(1:length(fBins),1:length(subjAnat),tMat_pow_cat); 
        set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
        colormap jet; colorbar;
        tit = strcat(subj,' All Bipolars');
        swagAxes(gca,10,'Frequency (Hz)','Bipolars',{'t(LF Move - LF WAIT)';tit});
        yt = [1:1:length(subjAnat)];
        set(gca,'ytick',yt,'yticklabel',subjAnat(yt),'yticklabelrotation',0)
        le_setFreqYTicks(gca,fBins)
        
        % Save this plot
        if savePlots
           if ~exist(plotDir,'dir'),mkdir(plotDir);end
           cd(plotDir);
           saveas(gca,'allBPs','png');
        end
    end
end