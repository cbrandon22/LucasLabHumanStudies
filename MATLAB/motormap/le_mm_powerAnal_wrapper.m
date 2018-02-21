% This script wraps around le_mm_selectMotorSites and le_mm_cat_tfPowPhase
% to create a stacked tpower plot.

subjList = {'HUP084','HUP087','HUP091','HUP092','HUP111'};%'HUP089','HUP090'
dirs = le_dirs('motormap');
config_mm = le_mm_config;

for i=1:length(subjList)
    subj = subjList{i};
    [eLbl_list_motor,eLbl_list_nonmotor] = le_mm_selectMotorSites (subj);
    if ~isempty(eLbl_list_motor)
        [subj_tPow,tStruct,subjAnat,config_pow,subjInsigCentroids,subjSigCentroids,subjInsigAnat] = le_mm_cat_tfPowPhase(subj,eLbl_list_motor,'motor',[],1);
        if i==1
           tStruct_lf = tStruct;
           anatAxis = subjAnat;
           insigAnat = subjInsigAnat;
           tMat_pow_cat = subj_tPow;
           InsigCentroids = subjInsigCentroids;
           sigCentroids = subjSigCentroids;
        else
            tStruct_lf = [tStruct_lf,tStruct];
            anatAxis = [anatAxis,subjAnat];
            insigAnat = [insigAnat,subjInsigAnat];
            tMat_pow_cat = cat(1,tMat_pow_cat,subj_tPow);
            InsigCentroids = [InsigCentroids;subjInsigCentroids];
            sigCentroids = [sigCentroids;subjSigCentroids];
        end
    end
    countPrimInsig = 0;
    for iii=1:length(subjInsigAnat)
        if ismember(subjInsigAnat{iii},{'L-PRE','L-POST','R-PRE','R-POST'})
            countPrimInsig = countPrimInsig+1;
        end
    end
end
countPrimInsig = 0;
for iii=1:length(insigAnat)
    if ismember(insigAnat{iii},{'L-PRE','L-POST','R-PRE','R-POST'})
        countPrimInsig = countPrimInsig+1;
    end
end

if config_mm.use_kmeans
   tMat_cat =  tMat_pow_cat; %[tMat_pow_cat zMat_phase_cat]; 

    % do k-means clustering on this matrix to sort these into
    % groups (clusId describes group membership)
    [clusId] = kmeans(tMat_cat,config_mm.kmeans_k,'replicates',config_mm.kmeans_iters);


    % evaluate num of clusters
    %eva = evalclusters(tMat_cat,'kmeans','CalinskiHarabasz','KList',[1:6])
%     [jumps,SUMD] = JumpClus(tMat_cat,config_mm.kmeans_numKtoEval);
%     swagFig;
%     plot(1:config_mm.kmeans_numKtoEval,SUMD,'-k','linewidth',3);
%     swagAxes(gca,20,'Number of clusters','Model error (summed distance)')
     [sortedClusId,sortIdx] = sort(clusId);

else % sort on some explicit feature
    [sortIdx,exp_feat] = le_sortByExplicit(tMat_pow_cat,...
    tStruct_lf,config_mm,config_pow);
end 
   
% sort mats
tMat_pow_cat = tMat_pow_cat(sortIdx,:);
%zMat_phase_cat = zMat_phase_cat(sortIdx,:);
anatAxis = anatAxis(:,sortIdx);
tStruct_lf = tStruct_lf(:,sortIdx);
sigCentroids = sigCentroids(sortIdx,:);


if config_mm.use_kmeans
    tit =    ['k-means k = ' num2str(config_mm.kmeans_k)];
else
    tit = ['sorted based on ' config_mm.explicit_feature_to_sort_by];
end

tBins = nanmean(config_pow.timeBins,2)';
fBins = [config_pow.freqBins];
fBins = fBins(tStruct_lf(1).fInd);

clusSeparationInd = [];
for i=2:length(sortedClusId)
    if sortedClusId(i) ~= sortedClusId(i-1)
        clusSeparationInd = [clusSeparationInd, i-1];
    end
end

tStruct_clus1 = tStruct_lf(1:clusSeparationInd(1));
tStruct_clus2 = tStruct_lf((clusSeparationInd(1)+1):length(tStruct_lf));

clus1Centroids = sigCentroids(1:clusSeparationInd(1),:);
clus2Centroids = sigCentroids((clusSeparationInd(1)+1):length(sigCentroids),:);

tMat_clus1 = tMat_pow_cat(1:clusSeparationInd(1),:);
tMat_clus2 = tMat_pow_cat((clusSeparationInd(1)+1):length(tMat_pow_cat),:);

% Cell array of eLbls
eLbls = cell(1,length(tStruct_lf));
for i=1: length(tStruct_lf)
    eLbls{i} = tStruct_lf(i).eLbl;
end

fclose('all');

keyboard
%%
% Create list of mono eLbls and elecLbls 
for i=1:length(tStruct_lf)
    tStruct_lf(i).clusID = sortedClusId(i);
end
elecLblcell = cell(length(tStruct_lf),2);
eLblcell = cell(length(tStruct_lf),2);
for i=1:length(tStruct_lf)
    if i==1
        subj = tStruct_lf(i).subj;
        subjAnatStruct = le_centroids2Anatstruct(subj);
        subjeLbls = cell(length(subjAnatStruct),1);
        for j=1:length(subjAnatStruct)
            subjeLbls(j) = {subjAnatStruct(j).eLbl};
        end
    end
    [~,locb] = ismember(tStruct_lf(i).eLbl,subjeLbls);
    elecLblcell(i,:)=strsplit(subjAnatStruct(locb).elecLbl,'-');
    if i<length(tStruct_lf)&&~strcmp(tStruct_lf(i).subj,tStruct_lf(i+1).subj)
        subj= tStruct_lf(i+1).subj;
        subjAnatStruct = le_centroids2Anatstruct(subj);
        subjeLbls = cell(length(subjAnatStruct),1);
        for j=1:length(subjAnatStruct)
            subjeLbls(j) = {subjAnatStruct(j).eLbl};
        end
    end
    eLblcell(i,:)=strsplit(tStruct_lf(i).eLbl,'-');
end
% Same for insignificant electrodes
elecLblcell = cell(length(InsigCentroids),2);
eLblcell = cell(length(InsigCentroids),2);
anatcell = cell(length(InsigCentroids),1);
for i=1:length(InsigCentroids)
    if i==1
        subj = InsigCentroids{i,1};
        subjAnatStruct = le_centroids2Anatstruct(subj);
        subjxyz = cell(length(subjAnatStruct),1);
        for j=1:length(subjAnatStruct)
            subjxyz{j} = strcat(num2str(subjAnatStruct(j).X,'%.4f'),',',num2str(subjAnatStruct(j).Y,'%.4f'),',',num2str(subjAnatStruct(j).Z,'%.4f'));
        end
    end
    xyzString = strcat(num2str(InsigCentroids{i,2},'%.4f'),',',num2str(InsigCentroids{i,3},'%.4f'),',',num2str(InsigCentroids{i,4},'%.4f'));
    [~,locb] = ismember(xyzString,subjxyz);
    elecLblcell(i,:)=strsplit(subjAnatStruct(locb).elecLbl,'-');
    eLblcell(i,:)=strsplit(subjAnatStruct(locb).eLbl,'-');
    anatcell{i}= subjAnatStruct(locb).anatAbbr;
    if i<length(InsigCentroids)&&~strcmp(InsigCentroids{i,1},InsigCentroids{i+1,1})
        subj= InsigCentroids{i+1,1};
        subjAnatStruct = le_centroids2Anatstruct(subj);
        subjxyz = cell(length(subjAnatStruct),1);
        for j=1:length(subjAnatStruct)
            subjxyz{j} = strcat(num2str(subjAnatStruct(j).X,'%.4f'),',',num2str(subjAnatStruct(j).Y,'%.4f'),',',num2str(subjAnatStruct(j).Z,'%.4f'));
        end
    end
end
%%
% Subclustering based on low frequency range only(less than 30.5284 Hz fbin(26))
tMat_clus = tMat_clus1; %pick cluster to subcluster
subAnatAxis = anatAxis(1:clusSeparationInd(1)); %clus1
%subAnatAxis = anatAxis((clusSeparationInd(1)+1):clusSeparationInd(2)); %clus2
%subAnatAxis = anatAxis((clusSeparationInd(2)+1):length(anatAxis)); %clus3
[subClusId] = kmeans(tMat_clus(:,1:26),3,'replicates',config_mm.kmeans_iters); %k=2
[sortedSubClusId,sortSubIdx] = sort(subClusId);
subClusSeparationInd = [];
for i=2:length(sortedSubClusId)
    if sortedSubClusId(i) ~= sortedSubClusId(i-1)
        subClusSeparationInd = [subClusSeparationInd, i-1];
    end
end
tMat_clus = tMat_clus(sortSubIdx,:);
subAnatAxis = subAnatAxis(sortSubIdx);
%subClusCentroids = clus1Centroids(sortSubIdx,:);
%subClusCentroids = clus2Centroids(sortSubIdx,:);
%subClusCentroids = clus3Centroids(sortSubIdx,:);
%plotSubClusID = sortedSubClusId + 10; %clus1
%plotSubClusID = sortedSubClusId + 30; %clus3

%plot subclustered k=2
%all electrodes
swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
imagesc(1:length(fBins),1:length(subAnatAxis),tMat_clus); 
set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
colormap jet; colorbar;
tit = 'Cluster 3 subclustered (k=3)';
swagAxes(gca,10,'Frequency (Hz)','Electrodes',{'t(LF LF MOVE - LF WAIT)';tit});
yt = [1:1:length(subAnatAxis)];
set(gca,'ytick',yt,'yticklabel',subAnatAxis(yt),'yticklabelrotation',0)
le_setFreqYTicks(gca,fBins)
%subClus1
swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
imagesc(1:length(fBins),1:length(subAnatAxis(1:subClusSeparationInd(1))),tMat_clus(1:subClusSeparationInd(1),:)); 
set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
colormap jet; colorbar;
tit = 'Cluster 3 subcluster1 (k=3)';
swagAxes(gca,10,'Frequency (Hz)','Electrodes',{'t(LF INSTRUCT - LF WAIT)';tit});
yt = [1:1:length(subAnatAxis(1:subClusSeparationInd(1)))];
set(gca,'ytick',yt,'yticklabel',subAnatAxis(yt),'yticklabelrotation',0)
le_setFreqYTicks(gca,fBins)
%subClus2 (k=2)
swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
imagesc(1:length(fBins),1:length(subAnatAxis((subClusSeparationInd(1)+1):length(subAnatAxis))),tMat_clus((subClusSeparationInd(1)+1):length(tMat_clus),:)); 
set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
colormap jet; colorbar;
tit = 'Cluster 1 subcluster2 (k=2)';
swagAxes(gca,10,'Frequency (Hz)','Electrodes',{'t(LF INSTRUCT - LF WAIT)';tit});
yt = [1:1:length(subAnatAxis((subClusSeparationInd(1)+1):length(subAnatAxis)))];
set(gca,'ytick',yt,'yticklabel',subAnatAxis(yt+(subClusSeparationInd(1))),'yticklabelrotation',0)
le_setFreqYTicks(gca,fBins)

%subClus2 (k=3)
swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
imagesc(1:length(fBins),1:length(subAnatAxis((subClusSeparationInd(1)+1):subClusSeparationInd(2))),tMat_clus((subClusSeparationInd(1)+1):subClusSeparationInd(2),:)); 
set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
colormap jet; colorbar;
tit = 'Cluster 3 subcluster2 (k=3)';
swagAxes(gca,10,'Frequency (Hz)','Electrodes',{'t(LF INSTRUCT - LF WAIT)';tit});
yt = [1:1:length(subAnatAxis((subClusSeparationInd(1)+1):subClusSeparationInd(2)))];
set(gca,'ytick',yt,'yticklabel',subAnatAxis(yt+(subClusSeparationInd(1))),'yticklabelrotation',0)
le_setFreqYTicks(gca,fBins)

%subClus3 (k=3)
swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
imagesc(1:length(fBins),1:length(subAnatAxis((subClusSeparationInd(2)+1):length(subAnatAxis))),tMat_clus((subClusSeparationInd(2)+1):length(tMat_clus),:)); 
set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
colormap jet; colorbar;
tit = 'Cluster 3 subcluster3 (k=3)';
swagAxes(gca,10,'Frequency (Hz)','Electrodes',{'t(LF INSTRUCT - LF WAIT)';tit});
yt = [1:1:length(subAnatAxis((subClusSeparationInd(2)+1):length(subAnatAxis)))];
set(gca,'ytick',yt,'yticklabel',subAnatAxis(yt+(subClusSeparationInd(2))),'yticklabelrotation',0)
le_setFreqYTicks(gca,fBins)
%%

% Color Plots (Frequency and t-value for selected electrodes)
% First cluster only
swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
imagesc(1:length(fBins),1:length(anatAxis(1:clusSeparationInd(1))),tMat_pow_cat(1:clusSeparationInd(1),:)); 
set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
colormap jet; colorbar;
tit = 'Cluster 1';
swagAxes(gca,10,'Frequency (Hz)','Electrode Pairs',{'t(LF Move - LF WAIT)';tit});
yt = [1:1:length(anatAxis(1:clusSeparationInd(1)))];
set(gca,'ytick',yt,'yticklabel',eLbls(yt),'yticklabelrotation',0)
le_setFreqYTicks(gca,fBins)

% Second cluster only
swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
imagesc(1:length(fBins),1:length(anatAxis((clusSeparationInd(1)+1):length(anatAxis))),tMat_pow_cat((clusSeparationInd(1)+1):length(tMat_pow_cat),:)); 
set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
colormap jet; colorbar;
tit = 'Cluster 2';
swagAxes(gca,10,'Frequency (Hz)','Electrode Pairs',{'t(LF Move - LF WAIT)';tit});
yt = [1:1:length(anatAxis((clusSeparationInd(1)+1):length(anatAxis)))];
set(gca,'ytick',yt,'yticklabel',eLbls(yt+clusSeparationInd(1)),'yticklabelrotation',0)
le_setFreqYTicks(gca,fBins)

% Third cluster only
% swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
% imagesc(1:length(fBins),1:length(anatAxis(:,(clusSeparationInd(2)+1):clusSeparationInd(3))),tMat_pow_cat((clusSeparationInd(2)+1):clusSeparationInd(3),:)); 
% set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
% colormap jet; colorbar;
% swagAxes(gca,10,'Frequency (Hz)','Electrodes',{'t(LF MOVE - LF WAIT)';tit});
% yt = [1:1:length(anatAxis(:,(clusSeparationInd(2)+1):clusSeparationInd(3)))];
% set(gca,'ytick',yt,'yticklabel',anatAxis(yt),'yticklabelrotation',0)
% le_setFreqYTicks(gca,fBins)

%Final cluster only (k=3)
% swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
% imagesc(1:length(fBins),1:length(anatAxis((clusSeparationInd(2)+1):length(anatAxis))),tMat_pow_cat((clusSeparationInd(2)+1):length(tMat_pow_cat),:)); 
% set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
% colormap jet; colorbar;
% tit = 'Cluster 3';
% swagAxes(gca,10,'Frequency (Hz)','Electrodes',{'t(LF INSTRUCT - LF WAIT)';tit});
% yt = [1:1:length(anatAxis((clusSeparationInd(2)+1):length(anatAxis)))];
% set(gca,'ytick',yt,'yticklabel',anatAxis(yt+clusSeparationInd(2)),'yticklabelrotation',0)
% le_setFreqYTicks(gca,fBins)

% Final cluster only (k=4)
% swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
% imagesc(1:length(fBins),1:length(anatAxis(:,(clusSeparationInd(3)+1):length(anatAxis))),tMat_pow_cat((clusSeparationInd(3)+1):length(tMat_pow_cat),:)); 
% set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
% colormap jet; colorbar;
% swagAxes(gca,10,'Frequency (Hz)','Electrodes',{'t(LF MOVE - LF WAIT)';tit});
% yt = [1:1:length(anatAxis(:,(clusSeparationInd(3)+1):length(anatAxis)))];
% set(gca,'ytick',yt,'yticklabel',anatAxis(yt),'yticklabelrotation',0)
% le_setFreqYTicks(gca,fBins)

% All clusters together
swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
imagesc(1:length(fBins),1:length(tStruct_lf),tMat_pow_cat); 
set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
colormap jet; colorbar;
tit = 'Motor-Related Bipolar Pairs';
swagAxes(gca,10,'Frequency (Hz)','Electrode Pairs',{'t(LF LF Move - LF WAIT)';tit});
yt = [1:1:length(tStruct_lf)];
set(gca,'ytick',yt,'yticklabel',eLbls(yt),'yticklabelrotation',0)
le_setFreqYTicks(gca,fBins)
clusSeparation = refline(0,clusSeparationInd(1)+0.5);
clusSeparation.Color = 'k';
clusSeparation.LineWidth = 2;

% Plot average +/-SEM t-values for clusters
figure
plotSEM(fBins,nanmean(tMat_clus2),SEM(tMat_clus2),'b')
%plotSEM(fBins,nanmean(tMat_clus1),SEM(tMat_clus1),'r') 
hold all
%plotSEM(fBins,nanmean(tMat_clus2),SEM(tMat_clus2),'b')
plotSEM(fBins,nanmean(tMat_clus1),SEM(tMat_clus1),'r')
%plotSEM(fBins,nanmean(tMat_clus3),SEM(tMat_clus3),'r') 
legend('Cluster1','Cluster2','Location','northeastoutside')

plot(fBins,nanmean(tMat_clus1),'--k','linewidth',2)
plot(fBins,nanmean(tMat_clus2),'--k','linewidth',2)
%plot(fBins,nanmean(tMat_clus3),'--k','linewidth',2)
tit = 'Cluster Average Power Difference Move-Wait';
swagAxes(gca,20,'Frequency (Hz)','Average Normalized Power',tit)

% Shuffle cluster
clusAnat = {tStruct_clus2.anatAx};
shuffleInd = randperm(length(tMat_clus2));
shuffleClus = tMat_clus2(shuffleInd,:);
shuffleAnat = clusAnat(shuffleInd);

swagFig([0.4492    0.0663    0.4266    0.8150]); hold all
imagesc(1:length(fBins),1:length(clusAnat),tMat_clus2); 
set(gca,'ydir','normal','clim',config_mm.tf_clim_lf)
colormap jet; colorbar;
tit = 'Cluster 3';
swagAxes(gca,10,'Frequency (Hz)','Electrodes',{'t(LF INSTRUCT - LF WAIT)';tit});
yt = [1:1:length(clusAnat)];
set(gca,'ytick',yt,'yticklabel',clusAnat,'yticklabelrotation',0)
le_setFreqYTicks(gca,fBins)