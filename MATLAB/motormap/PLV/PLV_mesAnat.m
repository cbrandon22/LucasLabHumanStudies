if ~exist('clus_info','var'),defineClusters_k2;end
task = 'motormap';
%freqRangeList = [{'delta600','theta250'};{'alpha200','beta100'};...
        %{'lowFreq250','highFreq50'};{'broadtheta250','broadtheta300'}];
freqRangeList = {'lowFreq250','highFreq50'};
catClus = [{'clus1','clus2','clus0','clush'}]; %clusters to concatenate
catClusLbl = 'all'; % label for titles
subjPlots = 0;
for jj=1:size(freqRangeList,1)
    iterations = 3000; %number of samples for bootstrap CI
    normPlv = 1; %Use raw/normalized plv
    anatPlot = 1; %Breakdown clusters by anatomical region
    roi2plot = 'Association'; %Which ROI to plot 'Association'
    
    trimPLV = 512; %ms to trim from beginning and end of PLV
    sampleWin = 512; %window size (samples)
    sampleStep = 256; %window step size (samples)
    cueInd = 1525-trimPLV; %Index of wait or go cue
    
    dirs = le_dirs(task);
    freqRange = freqRangeList(jj,:);
    plvDir = fullfile(dirs.scratch,'PHASE/k2/btwRegions',[freqRange{1} '_' freqRange{2}]);
    saveFigs = 1;
    figDir = fullfile(dirs.scratch,'figs/PLV/mes/anat/hedgesg',[freqRange{1} '_' freqRange{2}]);
    if ~exist(figDir,'dir'),mkdir(figDir);end

    cd(plvDir);
    subjList = {clus_info.subj};
    subjList = unique(subjList);
    mSubjVec={};
    nmSubjVec={};
    btwSubjVec={};
    m_hfMoveDat=[];
    m_hfWaitDat=[];
    m_lfMoveDat=[];
    m_lfWaitDat=[];
    nm_hfMoveDat=[];
    nm_hfWaitDat=[];
    nm_lfMoveDat=[];
    nm_lfWaitDat=[];
    btw_hfMoveDat=[];
    btw_hfWaitDat=[];
    btw_lfMoveDat=[];
    btw_lfWaitDat=[];
    %Loop through subjects/clusters to combine data
    for i=1:length(subjList)
        subj = subjList{i};
        for ii=1:length(allClus)
            %load this subj+clus
            cd(plvDir);
            clus = allClus{ii};
            if ismember(clus,catClus)
                if exist([subj '_clusPLV_combinedROI_' clus '.mat'],'file')==2
                    load([subj '_clusPLV_combinedROI_' clus '.mat']);
                end
                %add this cluster data
                if exist('motor_normGammaMove','var')
                    m_hfMoveDat = [m_hfMoveDat;motor_normGammaMove];
                    m_hfWaitDat = [m_hfWaitDat;motor_normGammaWait];
                    m_lfMoveDat = [m_lfMoveDat;motor_normThetaMove];
                    m_lfWaitDat = [m_lfWaitDat;motor_normThetaWait];
                    for j=1:size(motor_normGammaMove,1)
                        mSubjVec = [mSubjVec;subj];
                    end
                end
                if exist('nonMotor_normGammaMove','var')
                    nm_hfMoveDat = [nm_hfMoveDat;nonMotor_normGammaMove];
                    nm_hfWaitDat = [nm_hfWaitDat;nonMotor_normGammaWait];
                    nm_lfMoveDat = [nm_lfMoveDat;nonMotor_normThetaMove];
                    nm_lfWaitDat = [nm_lfWaitDat;nonMotor_normThetaWait];
                    for j=1:size(nonMotor_normGammaMove,1)
                        nmSubjVec = [nmSubjVec;subj];
                    end
                end
                if exist('btw_normGammaMove','var')
                    btw_hfMoveDat = [btw_hfMoveDat;btw_normGammaMove];
                    btw_hfWaitDat = [btw_hfWaitDat;btw_normGammaWait];
                    btw_lfMoveDat = [btw_lfMoveDat;btw_normThetaMove];
                    btw_lfWaitDat = [btw_lfWaitDat;btw_normThetaWait];
                    for j=1:size(btw_normGammaMove,1)
                        btwSubjVec = [btwSubjVec;subj];
                    end
                end
                clear nonMotor* motor* btw_Gamma* btw_norm* btw_Theta*
            end
        end
        %make subj anat PLV plots
        if subjPlots
            tt = t(trimPLV:(length(t)-2*trimPLV));
            %% Motor high frequency
            tit = [subj ' Motor High Frequency'];
            moveDat = mean(m_hfMoveDat,1)';
            waitDat = mean(m_hfWaitDat,1)';
            moveDat = moveDat(trimPLV:(length(moveDat)-2*trimPLV));
            waitDat = waitDat(trimPLV:(length(waitDat)-2*trimPLV));
            [m_hfMoveWait,m_hfMoveMove,m_hfWaitWait,m_hfWin] = PLV_mm_mesPlot(moveDat,waitDat,tt,sampleWin,sampleStep,cueInd,tit,iterations,saveFigs,figDir);
            %% Motor low frequency
            tit = [subj ' Motor Low Frequency'];
            moveDat = mean(m_lfMoveDat,1)';
            waitDat = mean(m_lfWaitDat,1)';
            moveDat = moveDat(trimPLV:(length(moveDat)-2*trimPLV));
            waitDat = waitDat(trimPLV:(length(waitDat)-2*trimPLV));
            [m_lfMoveWait,m_lfMoveMove,m_lfWaitWait,m_lfWin] = PLV_mm_mesPlot(moveDat,waitDat,tt,sampleWin,sampleStep,cueInd,tit,iterations,saveFigs,figDir);
            %% Nonmotor high frequency
            tit = [subj ' Nonmotor High Frequency'];
            moveDat = mean(nm_hfMoveDat,1)';
            waitDat = mean(nm_hfWaitDat,1)';
            moveDat = moveDat(trimPLV:(length(moveDat)-2*trimPLV));
            waitDat = waitDat(trimPLV:(length(waitDat)-2*trimPLV));
            [nm_hfMoveWait,nm_hfMoveMove,nm_hfWaitWait,nm_hfWin] = PLV_mm_mesPlot(moveDat,waitDat,tt,sampleWin,sampleStep,cueInd,tit,iterations,saveFigs,figDir);
            %% Nonmotor low frequency
            tit = [subj ' Motor Low Frequency'];
            moveDat = mean(nm_lfMoveDat,1)';
            waitDat = mean(nm_lfWaitDat,1)';
            moveDat = moveDat(trimPLV:(length(moveDat)-2*trimPLV));
            waitDat = waitDat(trimPLV:(length(waitDat)-2*trimPLV));
            [nm_lfMoveWait,nm_lfMoveMove,nm_lfWaitWait,nm_lfWin] = PLV_mm_mesPlot(moveDat,waitDat,tt,sampleWin,sampleStep,cueInd,tit,iterations,saveFigs,figDir);
            %% Between Regions high frequency
            tit = [subj ' Between Regions High Frequency'];
            moveDat = mean(btw_hfMoveDat,1)';
            waitDat = mean(btw_hfWaitDat,1)';
            moveDat = moveDat(trimPLV:(length(moveDat)-2*trimPLV));
            waitDat = waitDat(trimPLV:(length(waitDat)-2*trimPLV));
            [btw_hfMoveWait,btw_hfMoveMove,btw_hfWaitWait,btw_hfWin] = PLV_mm_mesPlot(moveDat,waitDat,tt,sampleWin,sampleStep,cueInd,tit,iterations,saveFigs,figDir);
            %% Between Regions low frequency
            tit = [subj ' Between Regions Low Frequency'];
            moveDat = mean(btw_lfMoveDat,1)';
            waitDat = mean(btw_lfWaitDat,1)';
            moveDat = moveDat(trimPLV:(length(moveDat)-2*trimPLV));
            waitDat = waitDat(trimPLV:(length(waitDat)-2*trimPLV));
            [btw_lfMoveWait,btw_lfMoveMove,btw_lfWaitWait,btw_lfWin] = PLV_mm_mesPlot(moveDat,waitDat,tt,sampleWin,sampleStep,cueInd,tit,iterations,saveFigs,figDir);
        end
    end
    %make Across-Subjs anat PLV plots
    t = t(trimPLV:(length(t)-2*trimPLV));
    %% Motor high frequency
    tit = [catClusLbl ' Primary Cortex High Frequency'];
    tempMove = m_hfMoveDat(:,trimPLV:(size(m_hfMoveDat,2)-2*trimPLV));
    tempWait = m_hfWaitDat(:,trimPLV:(size(m_hfWaitDat,2)-2*trimPLV));
    m_hfMoveMeans = [mean(tempMove(:,1:cueInd),2),mean(tempMove(:,cueInd+1:end),2)];
    m_hfWaitMeans = [mean(tempWait(:,1:cueInd),2),mean(tempWait(:,cueInd+1:end),2)];
    moveDat = mean(m_hfMoveDat,1)';
    waitDat = mean(m_hfWaitDat,1)';
    moveDat = moveDat(trimPLV:(length(moveDat)-2*trimPLV));
    waitDat = waitDat(trimPLV:(length(waitDat)-2*trimPLV));
    [m_hfMoveWait,m_hfMoveMove,m_hfWaitWait,m_hfWin] = PLV_mm_mesPlot(moveDat,waitDat,t,sampleWin,sampleStep,cueInd,tit,iterations);
    if saveFigs
        print(gcf,fullfile(figDir,tit),'-dpng');close
    end
    %% Motor low frequency
    tit = [catClusLbl ' Primary Cortex Low Frequency'];
    tempMove = m_lfMoveDat(:,trimPLV:(size(m_lfMoveDat,2)-2*trimPLV));
    tempWait = m_lfWaitDat(:,trimPLV:(size(m_lfWaitDat,2)-2*trimPLV));
    m_lfMoveMeans = [mean(tempMove(:,1:cueInd),2),mean(tempMove(:,cueInd+1:end),2)];
    m_lfWaitMeans = [mean(tempWait(:,1:cueInd),2),mean(tempWait(:,cueInd+1:end),2)];
    moveDat = mean(m_lfMoveDat,1)';
    waitDat = mean(m_lfWaitDat,1)';
    moveDat = moveDat(trimPLV:(length(moveDat)-2*trimPLV));
    waitDat = waitDat(trimPLV:(length(waitDat)-2*trimPLV));
    [m_lfMoveWait,m_lfMoveMove,m_lfWaitWait,m_lfWin] = PLV_mm_mesPlot(moveDat,waitDat,t,sampleWin,sampleStep,cueInd,tit,iterations);
    if saveFigs
        print(gcf,fullfile(figDir,tit),'-dpng');close
    end
    %% Nonmotor high frequency
    tit = [catClusLbl ' Association Cortex High Frequency'];
    tempMove = nm_hfMoveDat(:,trimPLV:(size(nm_hfMoveDat,2)-2*trimPLV));
    tempWait = nm_hfWaitDat(:,trimPLV:(size(nm_hfWaitDat,2)-2*trimPLV));
    nm_hfMoveMeans = [mean(tempMove(:,1:cueInd),2),mean(tempMove(:,cueInd+1:end),2)];
    nm_hfWaitMeans = [mean(tempWait(:,1:cueInd),2),mean(tempWait(:,cueInd+1:end),2)];
    moveDat = mean(nm_hfMoveDat,1)';
    waitDat = mean(nm_hfWaitDat,1)';
    moveDat = moveDat(trimPLV:(length(moveDat)-2*trimPLV));
    waitDat = waitDat(trimPLV:(length(waitDat)-2*trimPLV));
    [nm_hfMoveWait,nm_hfMoveMove,nm_hfWaitWait,nm_hfWin] = PLV_mm_mesPlot(moveDat,waitDat,t,sampleWin,sampleStep,cueInd,tit,iterations);
    if saveFigs
        print(gcf,fullfile(figDir,tit),'-dpng');close
    end
    %% Nonmotor low frequency
    tit = [catClusLbl ' Association Cortex Low Frequency'];
    tempMove = nm_lfMoveDat(:,trimPLV:(size(nm_lfMoveDat,2)-2*trimPLV));
    tempWait = nm_lfMoveDat(:,trimPLV:(size(nm_lfMoveDat,2)-2*trimPLV));
    nm_lfMoveMeans = [mean(tempMove(:,1:cueInd),2),mean(tempMove(:,cueInd+1:end),2)];
    nm_lfWaitMeans = [mean(tempWait(:,1:cueInd),2),mean(tempWait(:,cueInd+1:end),2)];
    moveDat = mean(nm_lfMoveDat,1)';
    waitDat = mean(nm_lfWaitDat,1)';
    moveDat = moveDat(trimPLV:(length(moveDat)-2*trimPLV));
    waitDat = waitDat(trimPLV:(length(waitDat)-2*trimPLV));
    [nm_lfMoveWait,nm_lfMoveMove,nm_lfWaitWait,nm_lfWin] = PLV_mm_mesPlot(moveDat,waitDat,t,sampleWin,sampleStep,cueInd,tit,iterations);
    if saveFigs
        print(gcf,fullfile(figDir,tit),'-dpng');close
    end
    %% Between Regions high frequency
    tit = [catClusLbl ' Between Regions High Frequency'];
    tempMove = btw_hfMoveDat(:,trimPLV:(size(btw_hfMoveDat,2)-2*trimPLV));
    tempWait = btw_hfWaitDat(:,trimPLV:(size(btw_hfWaitDat,2)-2*trimPLV));
    btw_hfMoveMeans = [mean(tempMove(:,1:cueInd),2),mean(tempMove(:,cueInd+1:end),2)];
    btw_hfWaitMeans = [mean(tempWait(:,1:cueInd),2),mean(tempWait(:,cueInd+1:end),2)];
    moveDat = mean(btw_hfMoveDat,1)';
    waitDat = mean(btw_hfWaitDat,1)';
    moveDat = moveDat(trimPLV:(length(moveDat)-2*trimPLV));
    waitDat = waitDat(trimPLV:(length(waitDat)-2*trimPLV));
    [btw_hfMoveWait,btw_hfMoveMove,btw_hfWaitWait,btw_hfWin] = PLV_mm_mesPlot(moveDat,waitDat,t,sampleWin,sampleStep,cueInd,tit,iterations);
    if saveFigs
        print(gcf,fullfile(figDir,tit),'-dpng');close
    end
    %% Between Regions low frequency
    tit = [catClusLbl ' Between Regions Cortex Low Frequency'];
    tempMove = btw_lfMoveDat(:,trimPLV:(size(btw_lfMoveDat,2)-2*trimPLV));
    tempWait = btw_lfWaitDat(:,trimPLV:(size(btw_lfWaitDat,2)-2*trimPLV));
    btw_lfMoveMeans = [mean(tempMove(:,1:cueInd),2),mean(tempMove(:,cueInd+1:end),2)];
    btw_lfWaitMeans = [mean(tempWait(:,1:cueInd),2),mean(tempWait(:,cueInd+1:end),2)];
    moveDat = mean(btw_lfMoveDat,1)';
    waitDat = mean(btw_lfWaitDat,1)';
    moveDat = moveDat(trimPLV:(length(moveDat)-2*trimPLV));
    waitDat = waitDat(trimPLV:(length(waitDat)-2*trimPLV));
    [btw_lfMoveWait,btw_lfMoveMove,btw_lfWaitWait,btw_lfWin] = PLV_mm_mesPlot(moveDat,waitDat,t,sampleWin,sampleStep,cueInd,tit,iterations);
    if saveFigs
        print(gcf,fullfile(figDir,tit),'-dpng');close
    end
end
cd(figDir)
save([catClusLbl 'Motor'], 'm_*', 'mSubjVec', 't')
save([catClusLbl 'Nonmotor'], 'nm_*','nmSubjVec', 't')
save([catClusLbl 'Between'], 'btw_*','btwSubjVec', 't')
clearvars -except clus1Struct clus2Struct clus0Struct clusHStruct useBipolar clus_info subjVec allClus subjList sub