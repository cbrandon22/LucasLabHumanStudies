subjList = {'HUP084','HUP087','HUP091','HUP092','HUP111'};

dirs = le_dirs;
config_mm = le_mm_config;
[eInfoStruct] = le_mm_readElecInfo;

for i=1:length(subjList)
    subj = subjList{i};
    
    % lGenerate anatStruct
    [anatStruct] = le_centroids2Anatstruct(subj);
    
    % load tStruct and configPow
    [tStruct_hfa,config_pow] = le_calcTs_wrapper(subj,[],'motormap',config_mm.comparison,...
        config_mm.hfa_range,config_mm.hfa_collapseFreqFlag,config_mm.powConfigNum);
    
    % concatonate p values from each elec to a single matrix
    if size(tStruct_hfa(1).pMat,1)>1
        error ('tStruct needs to be computed for a single HFA band...pMat needs to be a row vec for each elec')
    end
    pMat = cat(1,tStruct_hfa.pMat);
    
    % apply FDR correction on all p-values
    %pMat_sig_idx = pMat <= 0.001;
    pMat_sig_idx = pMat <= max(fdr(pMat,config_mm.fdr_p_thresh));   %0.001;
    
    % identify motor sites as those that have a sig p - value during the
    % post-movement interval
    
    % first, identify post-movement time bins
    tIdx_postMove = find(mean(config_pow.timeBins,2) > 0,1);
    
    % clean pMat by removing prestim values
    pMat = pMat(:,tIdx_postMove:end);
    pMat_sig_idx = pMat_sig_idx(:,tIdx_postMove:end);
    
    % identify motor electrodes
    motor_idx = sum(pMat_sig_idx,2)>config_mm.num_sig_bins_thresh;
    
    % get eLbls for this patient (same order as pMat)
    subj_eLbls = cell(length(anatStruct),1);
    for ii=1:length(anatStruct)
        subj_eLbls(ii) = {anatStruct(ii).eLbl};
    end
    % get eLbls organized by cluster
    eLbl_list_clus1 = {};
    eLbl_list_clus2 = {};
    eLbl_list_clus0 = {};
    for j=1:length(eInfoStruct)
        if strcmp(eInfoStruct(j).subj,subj) && eInfoStruct(j).cluster == 1
            [eLbl_list_clus1] = [eLbl_list_clus1; strcat(num2str(eInfoStruct(j).eLbl1),'-',num2str(eInfoStruct(j).eLbl2))];
        elseif strcmp(eInfoStruct(j).subj,subj) && eInfoStruct(j).cluster == 2
            [eLbl_list_clus2] = [eLbl_list_clus2; strcat(num2str(eInfoStruct(j).eLbl1),'-',num2str(eInfoStruct(j).eLbl2))];
        elseif strcmp(eInfoStruct(j).subj,subj) && eInfoStruct(j).cluster == 0
            [eLbl_list_clus0] = [eLbl_list_clus0; strcat(num2str(eInfoStruct(j).eLbl1),'-',num2str(eInfoStruct(j).eLbl2))];
        end
    end
    % sort pMats by cluster
    clus1pMat = zeros(length(eLbl_list_clus1),size(pMat,2));
    for j=1:length(eLbl_list_clus1)
        [~,eInd] = ismember(eLbl_list_clus1(j),subj_eLbls);
        clus1pMat(j,:) = pMat(eInd,:);
    end
    clus2pMat = zeros(length(eLbl_list_clus2),size(pMat,2));
    for j=1:length(eLbl_list_clus2)
        [~,eInd] = ismember(eLbl_list_clus2(j),subj_eLbls);
        clus2pMat(j,:) = pMat(eInd,:);
    end
    clus0pMat = zeros(length(eLbl_list_clus0),size(pMat,2));
    for j=1:length(eLbl_list_clus0)
        [~,eInd] = ismember(eLbl_list_clus0(j),subj_eLbls);
        clus0pMat(j,:) = pMat(eInd,:);
    end
    % learn to plot
    keyboard
    for j=1:size(clus1pMat,1)
        plot(mean(clus1pMat),'b')
        hold all
    end
    for j=1:size(clus2pMat,1)
        plot(mean(clus2pMat),'r')
        hold all
    end
    for j=1:size(clus0pMat,1)
        plot(mean(clus0pMat),'k')
        hold all
    end
end