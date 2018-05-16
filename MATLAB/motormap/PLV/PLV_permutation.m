% Permutation test for PLV for conservative statistical significance
% estimate. See Weaver, K. E. et al 2016 for similar technique
%
% We get a pre/post cue distribution of move/wait PLV (a max and min 
% distribution for a total of 8 distributions) by shuffling the time-series
% data in windows of size $shuffleWin. From each permutation, we select
% $threshValueCount PLV's (at each high/low extreme) to build each
% distribution. We then compare this to the real PLVs, selecting those
% above the 95th percentile (alpha 0.05) as significant.

%Parameters:
shuffleWin = 250; % window size to shuffle in each permutation (ms)
iterations = 20; % number of permutations
threshValueCount = 50; % number of extreme PLV's to save (e.g. top/bottom 50 values)
trimPLV = 512; % samples to trim from each end of PLV
alpha = 0.5; % significance level to compare with permutation distribution

if ~exist('clus_info','var'),defineClusters_k2;end
subjList = {'HUP084','HUP087','HUP091','HUP092','HUP111'};
for sub=1:length(subjList)
    for jj=1:length(allClus)
          freqStruct.labels = {'lowFreq250_highFreq50'};
          freqStruct.range1 = [3,33];
          freqStruct.range2 = [70,100];
          freqStruct.order1 = [250];
          freqStruct.order2 = [50];
        for thisFreq=1:length(freqStruct.labels)
            subj=subjList{sub};
            currClus=allClus{jj};
            dirs = le_dirs('motormap');
            offset=-3; %seconds
            baseTimeInd1=475; %samples
            baseTimeInd2=1525; %samples
            
            % create currStruct
            eval(sprintf('currClusStruct=%sStruct;',currClus));
            clusSubjs = {currClusStruct.subj};
            currStructInd = strcmp(subj,clusSubjs);
            currStruct = currClusStruct(currStructInd);
            for i=1:length(currStruct)
                datMove(i,:,:)=currStruct(i).datMove;
                datWait(i,:,:)=currStruct(i).datWait;
            end
            currDatMove=permute(datMove, [1 3 2]);
            currDatWait=permute(datWait, [1 3 2]);
            
            filtSpecT.range=freqStruct.range1(thisFreq,:);
            filtSpecT.order=freqStruct.order1(thisFreq);
            filtSpecG.range=freqStruct.range2(thisFreq,:);
            filtSpecG.order=freqStruct.order2(thisFreq);
            
            % Permutation loop
            for iter = 1:iterations
                
            end
            
            
        end  
    end 
end
% split into windows

% Permutation loop to shuffle windows

    %Select highest/lowest PLV's from each loop
    
clearvars -except 'allClus' 'clus0Struct' 'clus1Struct' 'clus2Struct' 'clus_info' 'clusHStruct' 'subjVec' 'useBipolar';