subj='HUP111';
dirs = le_dirs;

[eLbl_list,eLbl_list_nonmotor] = le_mm_selectMotorSites (subj);
config_mm = le_mm_config;
[tStruct_lf,config_pow] = le_calcTs_wrapper(subj,[],'motormap',config_mm.comparison,...
    config_mm.fRange,config_mm.collapseFreqFlag,config_mm.powConfigNum,eLbl_list);

powDiff=struct();
[anatStruct] = le_centroids2Anatstruct(subj);

for i=1:length(eLbl_list);  
    powDiffTemp=tStruct_lf(i).pow1-tStruct_lf(i).pow2; 
    powDiffTemp=squeeze(mean(powDiffTemp,2)); % collapse across time
    powDiff(i).pow=powDiffTemp;
    
    [~,Locb] = ismember(eLbl_list(i),{anatStruct.eLbl});
    powDiff(i).eLbl = eLbl_list(i);
    powDiff(i).anatAbbr = {anatStruct(Locb).anatAbbr};
    clear powDiffTemp;  
end

%figure; pcolor([1:length(tStruct_lf(1).fInd)],[1:length(eLbl_list)], real(powDiff(:,:))); colorbar; shading flat
%caxis ([-0.5 0.5])

save([dirs.scratch '/powDiff/' subj '_powDiff.mat'], 'powDiff')
