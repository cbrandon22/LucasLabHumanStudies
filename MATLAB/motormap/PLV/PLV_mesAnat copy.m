if ~exist('clus_info','var'),defineClusters_k2;end
task = 'motormap';
%freqRangeList = [{'delta600','theta250'};{'alpha200','beta100'};...
        %{'lowFreq250','highFreq50'};{'broadtheta250','broadtheta300'}];
freqRangeList = {'lowFreq250','highFreq50'};
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
    plvDir = fullfile(dirs.scratch,'PHASE/k2/norm1',[freqRange{1} '_' freqRange{2}]);
    saveFigs = 1;
    figDir = fullfile(dirs.scratch,'figs/PLV/mes/anat/hedgesg',[freqRange{1} '_' freqRange{2}]);
    if ~exist(figDir,'dir'),mkdir(figDir);end

    cd(plvDir);
    subjList = {clus_info.subj};
    subjList = unique(subjList);
    mSubjVec={};
    nmSubjVec={};
    
    %Loop through subjects/clusters to combine data
    for i=1:length(subjList)
        subj = subjList{i};
        m_hfMoveDat=[];
        m_hfWaitDat=[];
        m_lfMoveDat=[];
        m_lfWaitDat=[];
        nm_hfMoveDat=[];
        nm_hfWaitDat=[];
        nm_lfMoveDat=[];
        nm_lfWaitDat=[];
        for ii=1:length(allClus)
            %load this subj+clus
            cd(plvDir);
            clus = allClus{ii};
            if ~strcmp(clus,'clusH')
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
                clear nonMotor* motor*
            end
        end
        %make subj anat PLV plots
        if subjPlots
            tt = t(trimPLV:(length(t)-trimPLV));
            %% Motor high frequency
            tit = [subj ' Motor High Frequency'];
            moveDat = mean(m_hfMoveDat,1)';
            waitDat = mean(m_hfWaitDat,1)';
            moveDat = moveDat(trimPLV:(length(moveDat)-trimPLV));
            waitDat = waitDat(trimPLV:(length(waitDat)-trimPLV));
            [m_hfMoveWait,m_hfMoveMove,m_hfWaitWait] = PLV_mm_mesPlot(moveDat,waitDat,tt,sampleWin,sampleStep,cueInd,tit,iterations,saveFigs,figDir);
            %% Motor low frequency
            tit = [subj ' Motor Low Frequency'];
            moveDat = mean(m_lfMoveDat,1)';
            waitDat = mean(m_lfWaitDat,1)';
            moveDat = moveDat(trimPLV:(length(moveDat)-trimPLV));
            waitDat = waitDat(trimPLV:(length(waitDat)-trimPLV));
            [m_lfMoveWait,m_lfMoveMove,m_lfWaitWait] = PLV_mm_mesPlot(moveDat,waitDat,tt,sampleWin,sampleStep,cueInd,tit,iterations,saveFigs,figDir);
            %% Nonmotor high frequency
            tit = [subj ' Nonmotor High Frequency'];
            moveDat = mean(nm_hfMoveDat,1)';
            waitDat = mean(nm_hfWaitDat,1)';
            moveDat = moveDat(trimPLV:(length(moveDat)-trimPLV));
            waitDat = waitDat(trimPLV:(length(waitDat)-trimPLV));
            [nm_hfMoveWait,nm_hfMoveMove,nm_hfWaitWait] = PLV_mm_mesPlot(moveDat,waitDat,tt,sampleWin,sampleStep,cueInd,tit,iterations,saveFigs,figDir);
            %% Nonmotor low frequency
            tit = [subj ' Motor Low Frequency'];
            moveDat = mean(nm_lfMoveDat,1)';
            waitDat = mean(nm_lfWaitDat,1)';
            moveDat = moveDat(trimPLV:(length(moveDat)-trimPLV));
            waitDat = waitDat(trimPLV:(length(waitDat)-trimPLV));
            [nm_lfMoveWait,nm_lfMoveMove,nm_lfWaitWait] = PLV_mm_mesPlot(moveDat,waitDat,tt,sampleWin,sampleStep,cueInd,tit,iterations,saveFigs,figDir);
        end
    end
    %make Across-Subjs anat PLV plots
    t = t(trimPLV:(length(t)-trimPLV));
    %% Motor low frequency
    tit = 'All Motor Low Frequency';
    moveDat = mean(m_lfMoveDat,1)';
    waitDat = mean(m_lfWaitDat,1)';
    moveDat = moveDat(trimPLV:(length(moveDat)-trimPLV));
    waitDat = waitDat(trimPLV:(length(waitDat)-trimPLV));
    [m_lfMoveWait,m_lfMoveMove,m_lfWaitWait] = PLV_mm_mesPlot(moveDat,waitDat,t,sampleWin,sampleStep,cueInd,tit,iterations,saveFigs,figDir);
    %% Motor high frequency
    tit = 'All Motor High Frequency';
    moveDat = mean(m_hfMoveDat,1)';
    waitDat = mean(m_hfWaitDat,1)';
    moveDat = moveDat(trimPLV:(length(moveDat)-trimPLV));
    waitDat = waitDat(trimPLV:(length(waitDat)-trimPLV));
    try
        [m_hfMoveWait,m_hfMoveMove,m_hfWaitWait] = PLV_mm_mesPlot(moveDat,waitDat,t,sampleWin,sampleStep,cueInd,tit,iterations,saveFigs,figDir);
    catch ex
        ex.stack
    end
    %% Motor low frequency
%     tit = 'All Motor Low Frequency';
%     moveDat = mean(m_lfMoveDat,1)';
%     waitDat = mean(m_lfWaitDat,1)';
%     moveDat = moveDat(trimPLV:(length(moveDat)-trimPLV));
%     waitDat = waitDat(trimPLV:(length(waitDat)-trimPLV));
%     [m_lfMoveWait,m_lfMoveMove,m_lfWaitWait] = PLV_mm_mesPlot(moveDat,waitDat,t,sampleWin,sampleStep,cueInd,tit,iterations,saveFigs,figDir);
%     clear moveDat waitDat
    %% Nonmotor high frequency
    tit = 'All Nonmotor High Frequency';
    moveDat = mean(nm_hfMoveDat,1)';
    waitDat = mean(nm_hfWaitDat,1)';
    moveDat = moveDat(trimPLV:(length(moveDat)-trimPLV));
    waitDat = waitDat(trimPLV:(length(waitDat)-trimPLV));
    [nm_hfMoveWait,nm_hfMoveMove,nm_hfWaitWait] = PLV_mm_mesPlot(moveDat,waitDat,t,sampleWin,sampleStep,cueInd,tit,iterations,saveFigs,figDir);
    clear moveDat waitDat
    %% Nonmotor low frequency
    tit = 'All Motor Low Frequency';
    moveDat = mean(nm_lfMoveDat,1)';
    waitDat = mean(nm_lfWaitDat,1)';
    moveDat = moveDat(trimPLV:(length(moveDat)-trimPLV));
    waitDat = waitDat(trimPLV:(length(waitDat)-trimPLV));
    [nm_lfMoveWait,nm_lfMoveMove,nm_lfWaitWait] = PLV_mm_mesPlot(moveDat,waitDat,t,sampleWin,sampleStep,cueInd,tit,iterations,saveFigs,figDir);
end
clearvars -except clus1Struct clus2Struct clus0Struct clusHStruct useBipolar clus_info subjVec allClus subjList sub