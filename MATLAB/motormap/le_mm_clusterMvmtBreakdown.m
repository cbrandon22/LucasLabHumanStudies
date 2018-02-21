% This script breaks down each of the 3 clusters created from all motor electrodes
% (move vs wait) and creates color plots for each of the movement types
% This script requires you to change the way le_mm_cat_tfPowPhase chooses a
% comparison. when le_mm_cat_tfPowPhase calls the calcTs wrapper (currently line 25), it gets
% comparison from the config file. Change this to get config from the
% function input variable before running.

dirs = le_dirs;
cd (dirs.scratch);
load('moveWaitClusterStruct.mat');
config_mm = le_mm_config;

% Break struct into clusters
clusSeparationInd = [];
for i=2:length(tStruct_lf)
    if tStruct_lf(i).clusID ~= tStruct_lf(i-1).clusID
        clusSeparationInd = [clusSeparationInd, i-1];
    end
end

allDiff = [];
centroids = {};

for clusId=1:length(clusSeparationInd)+1
    % Separate into cluster groups
    if clusId==1
        tStruct_clus = tStruct_lf(1:clusSeparationInd(1));
    elseif clusId==2
        tStruct_clus = tStruct_lf((clusSeparationInd(1)+1):clusSeparationInd(2));
    elseif clusId==3
        tStruct_clus = tStruct_lf((clusSeparationInd(2)+1):length(tStruct_lf));
    end
    
    % build eLbl_list for each subject, call cat_tPowPhase, concatenate
    i=1;
    firstSubj = true;
    while i<=length(tStruct_clus)
        subj = tStruct_clus(i).subj;
        eLbl_list = [];
        while i<=length(tStruct_clus) && strcmp(subj,tStruct_clus(i).subj)
            eLbl_list = [eLbl_list, {tStruct_clus(i).eLbl}];
            i=i+1;
        end    
        [m_subj_tPow,m_tStruct,subjAnat,config_pow,~,sigCentroids] = le_mm_cat_tfPowPhase(subj,eLbl_list,'motor','mouthWait');
        %[l_subj_tPow,l_tStruct] = le_mm_cat_tfPowPhase(subj,eLbl_list,'motor','leftWait');
        %[r_subj_tPow,r_tStruct] = le_mm_cat_tfPowPhase(subj,eLbl_list,'motor','rightWait');
        [h_subj_tPow,h_tStruct] = le_mm_cat_tfPowPhase(subj,eLbl_list,'motor','handWait');
        if firstSubj
            firstSubj = false;
            mouth_tStruct = m_tStruct;
            %left_tStruct = l_tStruct;
            %right_tStruct = r_tStruct;
            hand_tStruct = h_tStruct;
            clusCentroids = sigCentroids;
            anatAxis = subjAnat;
            mouth_tMat_pow_cat = m_subj_tPow;
            %left_tMat_pow_cat = l_subj_tPow;
            %right_tMat_pow_cat = r_subj_tPow;
            hand_tMat_pow_cat = h_subj_tPow;
        else
            mouth_tStruct = [mouth_tStruct,m_tStruct];
            %left_tStruct = [left_tStruct,l_tStruct];
            %right_tStruct = [right_tStruct,r_tStruct];
            hand_tStruct = [hand_tStruct,h_tStruct];
            anatAxis = [anatAxis,subjAnat];
            mouth_tMat_pow_cat = cat(1,mouth_tMat_pow_cat,m_subj_tPow);
            %left_tMat_pow_cat = cat(1,left_tMat_pow_cat,l_subj_tPow);
            %right_tMat_pow_cat = cat(1,right_tMat_pow_cat,r_subj_tPow);
            hand_tMat_pow_cat = cat(1,hand_tMat_pow_cat,h_subj_tPow);
            clusCentroids = [clusCentroids;sigCentroids];
        end
    end
    
    % k means sort
%     mouth_tMat_cat =  [mouth_tMat_pow_cat];
%     %left_tMat_cat =  [left_tMat_pow_cat];
%     %right_tMat_cat =  [right_tMat_pow_cat];
%     hand_tMat_cat =  [hand_tMat_pow_cat];
%     
%     [m_clusId] = kmeans(mouth_tMat_cat,config_mm.kmeans_k,'replicates',config_mm.kmeans_iters);
%     %[l_clusId] = kmeans(left_tMat_cat,config_mm.kmeans_k,'replicates',config_mm.kmeans_iters);
%     %[r_clusId] = kmeans(right_tMat_cat,config_mm.kmeans_k,'replicates',config_mm.kmeans_iters);
%     [h_clusId] = kmeans(hand_tMat_cat,config_mm.kmeans_k,'replicates',config_mm.kmeans_iters);
%     
%     [~,m_sortIdx] = sort(m_clusId);
%     %[~,l_sortIdx] = sort(l_clusId);
%     %[~,r_sortIdx] = sort(r_clusId);
%     [~,h_sortIdx] = sort(h_clusId);
%     
%     % sort mats
%     mouth_tMat_pow_cat = mouth_tMat_pow_cat(m_sortIdx,:);
%     %left_tMat_pow_cat = left_tMat_pow_cat(l_sortIdx,:);
%     %right_tMat_pow_cat = right_tMat_pow_cat(r_sortIdx,:);
%     hand_tMat_pow_cat = hand_tMat_pow_cat(h_sortIdx,:);
%     m_anatAxis = anatAxis(:,m_sortIdx);
%     %l_anatAxis = anatAxis(:,l_sortIdx);
%     %r_anatAxis = anatAxis(:,r_sortIdx);
%     h_anatAxis = anatAxis(:,h_sortIdx);
%     %tStruct_lf = tStruct_lf(:,sortIdx);

    
    % make plots for this cluster
    tBins = nanmean(config_pow.timeBins,2)';
    fBins = [config_pow.freqBins];
    fBins = fBins(tStruct_lf(1).fInd);
    
    % mouth
%     swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
%     imagesc(1:length(fBins),1:length(anatAxis),mouth_tMat_pow_cat); 
%     set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
%     colormap jet; colorbar;
%     tit = strcat('Cluster ', num2str(clusId),' Mouth');
%     swagAxes(gca,10,'Frequency (Hz)','Electrodes',{'t(LF MOUTH - LF WAIT)';tit});
%     yt = [1:1:length(anatAxis)];
%     set(gca,'ytick',yt,'yticklabel',anatAxis(yt),'yticklabelrotation',0)
%     le_setFreqYTicks(gca,fBins)
%     
    % left
%     swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
%     imagesc(1:length(fBins),1:length(anatAxis),left_tMat_pow_cat); 
%     set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
%     colormap jet; colorbar;
%     tit = strcat('Cluster ', num2str(clusId),' Left');
%     swagAxes(gca,10,'Frequency (Hz)','Electrodes',{'t(LF LEFT - LF WAIT)';tit});
%     yt = [1:1:length(anatAxis)];
%     set(gca,'ytick',yt,'yticklabel',anatAxis(yt),'yticklabelrotation',0)
%     le_setFreqYTicks(gca,fBins)
    
    % right
%     swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
%     imagesc(1:length(fBins),1:length(anatAxis),right_tMat_pow_cat); 
%     set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
%     colormap jet; colorbar;
%     tit = strcat('Cluster ', num2str(clusId),' Right');
%     swagAxes(gca,10,'Frequency (Hz)','Electrodes',{'t(LF RIGHT - LF WAIT)';tit});
%     yt = [1:1:length(anatAxis)];
%     set(gca,'ytick',yt,'yticklabel',anatAxis(yt),'yticklabelrotation',0)
%     le_setFreqYTicks(gca,fBins)  
%     
    % hands
%     swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
%     imagesc(1:length(fBins),1:length(anatAxis),hand_tMat_pow_cat); 
%     set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
%     colormap jet; colorbar;
%     tit = strcat('Cluster ', num2str(clusId),' Hands');
%     swagAxes(gca,10,'Frequency (Hz)','Electrodes',{'t(LF HANDS - LF WAIT)';tit});
%     yt = [1:1:length(anatAxis)];
%     set(gca,'ytick',yt,'yticklabel',anatAxis(yt),'yticklabelrotation',0)
%     le_setFreqYTicks(gca,fBins)  
    
    % Calculate difference between hands and mouth t-stats, take abs value,
    % sum across rows, sort on difference
    tMatDiff = sum(abs(mouth_tMat_pow_cat - hand_tMat_pow_cat),2);
    [sortedDiff,sortInd] = sort(tMatDiff);
    diffSortHands = hand_tMat_pow_cat(sortInd,:);
    diffSortMouth = mouth_tMat_pow_cat(sortInd,:);
    diffAnatAxis = anatAxis(sortInd);
    clusCentroids = clusCentroids(sortInd,:);
    
    for ii=1:length(diffAnatAxis)
        diffAnatAxis(ii) = strcat(diffAnatAxis(ii),' (',num2str(sortedDiff(ii)),')');
    end
    % make difference sorted plots
    swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
    imagesc(1:length(fBins),1:length(diffAnatAxis),diffSortHands); 
    set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
    colormap jet; colorbar;
    tit = strcat('Cluster ', num2str(clusId),' Hands');
    swagAxes(gca,10,'Frequency (Hz)','Electrodes',{'t(LF HANDS - LF WAIT)';tit});
    yt = [1:1:length(diffAnatAxis)];
    set(gca,'ytick',yt,'yticklabel',diffAnatAxis(yt),'yticklabelrotation',0)
    le_setFreqYTicks(gca,fBins) 
    
    swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
    imagesc(1:length(fBins),1:length(diffAnatAxis),diffSortMouth); 
    set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
    colormap jet; colorbar;
    tit = strcat('Cluster ', num2str(clusId),' Mouth');
    swagAxes(gca,10,'Frequency (Hz)','Electrodes',{'t(LF Mouth - LF WAIT)';tit});
    yt = [1:1:length(diffAnatAxis)];
    set(gca,'ytick',yt,'yticklabel',diffAnatAxis(yt),'yticklabelrotation',0)
    le_setFreqYTicks(gca,fBins)
    
    % all difference values for each cluster
    allDiff = [allDiff;sortedDiff];
    
    %Select only bp's in upper quartile of difference
    upperQuartile = 75.4216; % had to determine earlier using complete allDiff
    mouth_specific_tMat = [];
    hand_specific_tMat = [];
    non_specific_tMat = [];
    mouthAnat = {};
    handAnat = {};
    nonSpecificAnat = {};
    addClusId = clusId*10;
    for jj=1:length(sortedDiff)
        diff = sum(diffSortMouth(jj,26:length(fBins)) - diffSortHands(jj,26:length(fBins)));
        if diff > 0 && sum(diffSortMouth(jj,26:length(fBins)))>25 && jj>find(sortedDiff>upperQuartile,1)% mouth electrode
            mouth_specific_tMat = [mouth_specific_tMat; diffSortMouth(jj,:)];
            clusCentroids(jj,5) = {1+addClusId};
            mouthAnat = [mouthAnat,diffAnatAxis(jj)];
        elseif diff < 0 && sum(diffSortHands(jj,26:length(fBins)))>25 && jj>find(sortedDiff>upperQuartile,1)% hand electrode
            hand_specific_tMat = [hand_specific_tMat; diffSortHands(jj,:)];
            clusCentroids(jj,5) = {2+addClusId};
            handAnat = [handAnat,diffAnatAxis(jj)];
        elseif sum(diffSortMouth(jj,26:length(fBins)))>25 && sum(diffSortHands(jj,26:length(fBins)))>25
            non_specific_tMat = [non_specific_tMat;squeeze(nanmean(tStruct_clus(jj).tMat,2))']; 
            clusCentroids(jj,5) = {3+addClusId};
            nonSpecificAnat = [nonSpecificAnat,diffAnatAxis(jj)];
        end
    end
    if ~isempty(mouth_specific_tMat) || ~isempty(hand_specific_tMat)
        saveCentroidInd = [];
        for cent=1:length(clusCentroids)
            if ~isempty(clusCentroids{cent,5})
                saveCentroidInd = [saveCentroidInd,cent];
            end
        end
        clusCentroids = clusCentroids(saveCentroidInd,:);
        if clusId ==1
            centroids = clusCentroids;
        else
            centroids = [centroids;clusCentroids];
        end
    end
    % Plot hand, mouth, and nonspecific tMats
    swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
    imagesc(1:length(fBins),1:length(handAnat),hand_specific_tMat); 
    set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
    colormap jet; colorbar;
    tit = strcat('Cluster ', num2str(clusId),' hf hand');
    swagAxes(gca,10,'Frequency (Hz)','Electrodes',{'t(LF HANDS - LF WAIT)';tit});
    yt = [1:1:length(handAnat)];
    set(gca,'ytick',yt,'yticklabel',handAnat(yt),'yticklabelrotation',0)
    le_setFreqYTicks(gca,fBins)
    
    swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
    imagesc(1:length(fBins),1:length(mouthAnat),mouth_specific_tMat); 
    set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
    colormap jet; colorbar;
    tit = strcat('Cluster ', num2str(clusId),' hf mouth');
    swagAxes(gca,10,'Frequency (Hz)','Electrodes',{'t(LF HANDS - LF WAIT)';tit});
    yt = [1:1:length(mouthAnat)];
    set(gca,'ytick',yt,'yticklabel',mouthAnat(yt),'yticklabelrotation',0)
    le_setFreqYTicks(gca,fBins)
    
    swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
    imagesc(1:length(fBins),1:length(nonSpecificAnat),non_specific_tMat); 
    set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
    colormap jet; colorbar;
    tit = strcat('Cluster ', num2str(clusId),' hf hands and mouth');
    swagAxes(gca,10,'Frequency (Hz)','Electrodes',{'t(LF MOVE - LF WAIT)';tit});
    yt = [1:1:length(nonSpecificAnat)];
    set(gca,'ytick',yt,'yticklabel',nonSpecificAnat(yt),'yticklabelrotation',0)
    le_setFreqYTicks(gca,fBins)
   
end
fclose('all');

