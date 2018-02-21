% Perform 2 sample permutation test on motormap PLV
if ~exist('clus_info','var'),defineClusters_k2;end
freqRangeList = [{'delta600','theta250'};{'alpha200','beta100'};...
        {'lowFreq250','highFreq50'};{'broadtheta250','broadtheta300'}];
summaryStruct = struct();
for jj=1:size(freqRangeList,1)
    iterations = 1000; %number of samples for bootstrap
    shuffles = []; %default shuffle
    doWilcox = 1; %non-parametric
    doSig = 1;
    
    baseTimeInd1=475; %samples
    baseTimeInd2=1525; %samples
    
    dirs = le_dirs;
    freqRange = freqRangeList(jj,:)
    plvDir = fullfile(dirs.scratch,'PHASE/k2/',[freqRange{1} '_' freqRange{2}]);
    statDir = fullfile(dirs.scratch,'PHASE/stats/2sec',[freqRange{1} '_' freqRange{2}]);
    sigStatDir = fullfile(dirs.scratch,'PHASE/stats/2secSig',[freqRange{1} '_' freqRange{2}]);
    saveStats = 1;
    if ~exist(statDir,'dir'),mkdir(statDir);end
    if ~exist(sigStatDir,'dir'),mkdir(sigStatDir);end

    cd(plvDir);
    subjList = {clus_info.subj};
    subjList = unique(subjList);
    
    allclus1_mNormGammaMove = [];
    allclus1_mNormGammaWait = [];
    allclus1_mNormThetaMove = [];
    allclus1_mNormThetaWait = [];
    allclus2_mNormGammaMove = [];
    allclus2_mNormGammaWait = [];
    allclus2_mNormThetaMove = [];
    allclus2_mNormThetaWait = [];
    
    for i=1:length(subjList)
        for ii=1:length(allClus)
            %load this subj+clus
            cd(plvDir);
            subj = subjList{i};
            clus = allClus{ii};
            load([subj '_clusPLV_combinedROI_' clus '.mat']);
            
            if exist('motor_normGammaMove','var')
                %bootstrap setup
                gammaMoveDat = mean(motor_GammaMove,1)';
                gammaWaitDat = mean(motor_GammaWait,1)';
                gammaMoveBase = gammaMoveDat(baseTimeInd1:baseTimeInd2);
                gammaWaitBase = gammaWaitDat(baseTimeInd1:baseTimeInd2);
                gammaMoveData = gammaMoveDat(baseTimeInd2:end-512);
                gammaWaitData = gammaWaitDat(baseTimeInd2:end-512);
                disp([subj ': ' clus ' ' freqRange{2} ' move']);
                %[p_boot,p_sig,shuffles] = bootstrap(iterations, gammaMoveBase, gammaMoveData, shuffles, doWilcox, doSig);
                [rs_boot,rs_stat,shuffles] = bootstrap(iterations, gammaMoveBase, gammaMoveData, shuffles, doWilcox, doSig);
                n = length(find(abs(cell2mat({rs_boot.zval}))>=abs(rs_stat.zval)));
                p_val = (n+1)/(iterations+1);
                if size(summaryStruct,1) < 2;
                    summaryStruct.subj = subj;
                    summaryStruct.clus = clus;
                    summaryStruct.freq = freqRange{2};
                    summaryStruct.stim = 'move';
                    summaryStruct.z_boot = {rs_boot.zval};
                    summaryStruct.z_dat = rs_stat.zval;
                    summaryStruct.p_value = p_val;
                else
                    currStruct.subj = subj;
                    currStruct.clus = clus;
                    currStruct.freq = freqRange{2};
                    currStruct.stim = 'move';
                    currStruct.z_boot = {rs_boot.zval};
                    currStruct.z_dat = rs_stat.zval;
                    currStruct.p_value = p_val;
                    summaryStruct = [summaryStruct;currStruct];
                end
                if saveStats
                    cd(statDir);
                    save([subj '_' clus '_' freqRange{2} '_moveBoot.mat'], 'p_val','rs_boot', 'rs_stat','shuffles')
                end
                disp([subj ': ' clus ' ' freqRange{2} ' wait']);
                %[p_boot,p_sig,shuffles] = bootstrap(iterations, gammaWaitBase, gammaWaitData, shuffles, doWilcox, doSig);
                [rs_boot,rs_stat,shuffles] = bootstrap(iterations, gammaWaitBase, gammaWaitData, shuffles, doWilcox, doSig);
                n = length(find(abs(cell2mat({rs_boot.zval}))>=abs(rs_stat.zval)));
                p_val = (n+1)/(iterations+1);
                currStruct.subj = subj;
                currStruct.clus = clus;
                currStruct.freq = freqRange{2};
                currStruct.stim = 'wait';
                currStruct.z_boot = {rs_boot.zval};
                currStruct.z_dat = rs_stat.zval;
                currStruct.p_value = p_val;
                summaryStruct = [summaryStruct;currStruct];
                if saveStats
                    cd(statDir);
                    save([subj '_' clus '_' freqRange{2} '_waitBoot.mat'], 'p_val','rs_boot', 'rs_stat','shuffles')
                end
                if ii ==1
                    allclus1_mNormGammaMove = [allclus1_mNormGammaMove; motor_GammaMove];
                    allclus1_mNormGammaWait = [allclus1_mNormGammaWait; motor_GammaWait];
                else
                    allclus2_mNormGammaMove = [allclus2_mNormGammaMove; motor_GammaMove];
                    allclus2_mNormGammaWait = [allclus2_mNormGammaWait; motor_GammaWait];
                end
            end
            
            if exist('motor_normThetaMove','var')
                thetaMoveDat = mean(motor_ThetaMove,1)';
                thetaWaitDat = mean(motor_ThetaWait,1)';
                thetaMoveBase = thetaMoveDat(baseTimeInd1:baseTimeInd2);
                thetaWaitBase = thetaWaitDat(baseTimeInd1:baseTimeInd2);
                thetaMoveData = thetaMoveDat(baseTimeInd2:end-512);
                thetaWaitData = thetaWaitDat(baseTimeInd2:end-512);
                disp([subj ': ' clus ' ' freqRange{1} ' move']);
                %[p_boot,p_sig,shuffles] = bootstrap(iterations, thetaMoveBase, thetaMoveData, shuffles, doWilcox, doSig);
                [rs_boot,rs_stat,shuffles] = bootstrap(iterations, thetaMoveBase, thetaMoveData, shuffles, doWilcox, doSig);
                n = length(find(abs(cell2mat({rs_boot.zval}))>=abs(rs_stat.zval)));
                p_val = (n+1)/(iterations+1);
                currStruct.subj = subj;
                currStruct.clus = clus;
                currStruct.freq = freqRange{1};
                currStruct.stim = 'move';
                currStruct.z_boot = {rs_boot.zval};
                currStruct.z_dat = rs_stat.zval;
                currStruct.p_value = p_val;
                summaryStruct = [summaryStruct;currStruct];
                if saveStats
                    cd(statDir);
                    save([subj '_' clus '_' freqRange{1} '_moveBoot.mat'], 'p_val','rs_boot', 'rs_stat','shuffles')
                end
                disp([subj ': ' clus ' ' freqRange{1} ' wait']);
                %[p_boot,p_sig,shuffles] = bootstrap(iterations, thetaWaitBase, thetaWaitData, shuffles, doWilcox, doSig);
                [rs_boot,rs_stat,shuffles] = bootstrap(iterations, thetaWaitBase, thetaWaitData, shuffles, doWilcox, doSig);
                n = length(find(abs(cell2mat({rs_boot.zval}))>=abs(rs_stat.zval)));
                p_val = (n+1)/(iterations+1);
                currStruct.subj = subj;
                currStruct.clus = clus;
                currStruct.freq = freqRange{1};
                currStruct.stim = 'wait';
                currStruct.z_boot = {rs_boot.zval};
                currStruct.z_dat = rs_stat.zval;
                currStruct.p_value = p_val;
                summaryStruct = [summaryStruct;currStruct];
                if saveStats
                    cd(statDir);
                    save([subj '_' clus '_' freqRange{1} '_waitBoot.mat'], 'p_val','rs_boot', 'rs_stat','shuffles')
                end
                if ii ==1
                    allclus1_mNormThetaMove = [allclus1_mNormThetaMove; motor_ThetaMove];
                    allclus1_mNormThetaWait = [allclus1_mNormThetaWait; motor_ThetaWait];
                else
                    allclus2_mNormThetaMove = [allclus2_mNormThetaMove; motor_ThetaMove];
                    allclus2_mNormThetaWait = [allclus2_mNormGammaWait; motor_ThetaWait];
                end
            end
            clearvars -except clus1Struct clus2Struct useBipolar clus_info subjVec allClus i ii plvDir subjList saveStats statDir sigStatDir allclus* t freqRange* base* iterations doWilcox doSig shuffles pVec summaryStruct
        end
    end
    
    %All Cluster 1
    gammaMoveDat = mean(allclus1_mNormGammaMove,1)';
    gammaWaitDat = mean(allclus1_mNormGammaWait,1)';
    gammaMoveBase = gammaMoveDat(baseTimeInd1:baseTimeInd2);
    gammaWaitBase = gammaWaitDat(baseTimeInd1:baseTimeInd2);
    gammaMoveData = gammaMoveDat(baseTimeInd2:end-512);
    gammaWaitData = gammaWaitDat(baseTimeInd2:end-512);
    disp(['allClus1 ' freqRange{2} ' move']);
    %[p_boot,p_sig,shuffles] = bootstrap(iterations, gammaMoveBase, gammaMoveData, shuffles, doWilcox, doSig);
    [rs_boot,rs_stat,shuffles] = bootstrap(iterations, gammaMoveBase, gammaMoveData, shuffles, doWilcox, doSig);
    n = length(find(abs(cell2mat({rs_boot.zval}))>=abs(rs_stat.zval)));
    p_val = (n+1)/(iterations+1);
    currStruct.subj = 'allsubj';
    currStruct.clus = 'allclus1';
    currStruct.freq = freqRange{2};
    currStruct.stim = 'move';
    currStruct.z_boot = {rs_boot.zval};
    currStruct.z_dat = rs_stat.zval;
    currStruct.p_value = p_val;
    summaryStruct = [summaryStruct;currStruct];
    if saveStats
        cd(statDir);
        save(['allClus1_' freqRange{2} '_moveBoot.mat'], 'p_val','rs_boot', 'rs_stat','shuffles')
    end
    disp(['allClus1 ' freqRange{2} ' wait']);
    %[p_boot,p_sig,shuffles] = bootstrap(iterations, gammaWaitBase, gammaWaitData, shuffles, doWilcox, doSig);
    [rs_boot,rs_stat,shuffles] = bootstrap(iterations, gammaWaitBase, gammaWaitData, shuffles, doWilcox, doSig);
    n = length(find(abs(cell2mat({rs_boot.zval}))>=abs(rs_stat.zval)));
    p_val = (n+1)/(iterations+1);
    currStruct.subj = 'allsubj';
    currStruct.clus = 'allclus1';
    currStruct.freq = freqRange{2};
    currStruct.stim = 'wait';
    currStruct.z_boot = {rs_boot.zval};
    currStruct.z_dat = rs_stat.zval;
    currStruct.p_value = p_val;
    summaryStruct = [summaryStruct;currStruct];
    if saveStats
        cd(statDir);
        save(['allClus1_' freqRange{2} '_waitBoot.mat'], 'p_val','rs_boot', 'rs_stat','shuffles')
    end
    
    thetaMoveDat = mean(allclus1_mNormThetaMove,1)';
    thetaWaitDat = mean(allclus1_mNormThetaWait,1)';
    thetaMoveBase = thetaMoveDat(baseTimeInd1:baseTimeInd2);
    thetaWaitBase = thetaWaitDat(baseTimeInd1:baseTimeInd2);
    thetaMoveData = thetaMoveDat(baseTimeInd2:end-512);
    thetaWaitData = thetaWaitDat(baseTimeInd2:end-512);
    disp(['allClus1 ' freqRange{1} ' move']);
    %[p_boot,p_sig,shuffles] = bootstrap(iterations, thetaMoveBase, thetaMoveData, shuffles, doWilcox, doSig);
    [rs_boot,rs_stat,shuffles] = bootstrap(iterations, thetaMoveBase, thetaMoveData, shuffles, doWilcox, doSig);
    n = length(find(abs(cell2mat({rs_boot.zval}))>=abs(rs_stat.zval)));
    p_val = (n+1)/(iterations+1);
    currStruct.subj = 'allsubj';
    currStruct.clus = 'allclus1';
    currStruct.freq = freqRange{1};
    currStruct.stim = 'move';
    currStruct.z_boot = {rs_boot.zval};
    currStruct.z_dat = rs_stat.zval;
    currStruct.p_value = p_val;
    summaryStruct = [summaryStruct;currStruct];
    if saveStats
        cd(statDir);
        save(['allClus1_' freqRange{1} '_moveBoot.mat'], 'p_val','rs_boot', 'rs_stat','shuffles')
    end
    disp(['allClus1 ' freqRange{1} ' wait']);
    %[p_boot,p_sig,shuffles] = bootstrap(iterations, thetaWaitBase, thetaWaitData, shuffles, doWilcox, doSig);
    [rs_boot,rs_stat,shuffles] = bootstrap(iterations, thetaWaitBase, thetaWaitData, shuffles, doWilcox, doSig);
    n = length(find(abs(cell2mat({rs_boot.zval}))>=abs(rs_stat.zval)));
    p_val = (n+1)/(iterations+1);
    currStruct.subj = 'allsubj';
    currStruct.clus = 'allclus1';
    currStruct.freq = freqRange{1};
    currStruct.stim = 'wait';
    currStruct.z_boot = {rs_boot.zval};
    currStruct.z_dat = rs_stat.zval;
    currStruct.p_value = p_val;
    summaryStruct = [summaryStruct;currStruct];
    if saveStats
        cd(statDir);
        save(['allClus1_' freqRange{1} '_waitBoot.mat'], 'p_val','rs_boot', 'rs_stat','shuffles')
    end
    
    %All cluster 2
    gammaMoveDat = mean(allclus2_mNormGammaMove,1)';
    gammaWaitDat = mean(allclus2_mNormGammaWait,1)';
    gammaMoveBase = gammaMoveDat(baseTimeInd1:baseTimeInd2);
    gammaWaitBase = gammaWaitDat(baseTimeInd1:baseTimeInd2);
    gammaMoveData = gammaMoveDat(baseTimeInd2:end-512);
    gammaWaitData = gammaWaitDat(baseTimeInd2:end-512);
    disp(['allClus2 ' freqRange{2} ' move']);
    %[p_boot,p_sig,shuffles] = bootstrap(iterations, gammaMoveBase, gammaMoveData, shuffles, doWilcox, doSig);
    [rs_boot,rs_stat,shuffles] = bootstrap(iterations, gammaMoveBase, gammaMoveData, shuffles, doWilcox, doSig);
    n = length(find(abs(cell2mat({rs_boot.zval}))>=abs(rs_stat.zval)));
    p_val = (n+1)/(iterations+1);
    currStruct.subj = 'allsubj';
    currStruct.clus = 'allclus2';
    currStruct.freq = freqRange{2};
    currStruct.stim = 'move';
    currStruct.z_boot = {rs_boot.zval};
    currStruct.z_dat = rs_stat.zval;
    currStruct.p_value = p_val;
    summaryStruct = [summaryStruct;currStruct];
    if saveStats
        cd(statDir);
        save(['allClus2_' freqRange{2} '_moveBoot.mat'], 'p_val','rs_boot', 'rs_stat','shuffles')
    end
    disp(['allClus2 ' freqRange{2} ' wait']);
    %[p_boot,p_sig,shuffles] = bootstrap(iterations, gammaWaitBase, gammaWaitData, shuffles, doWilcox, doSig);
    [rs_boot,rs_stat,shuffles] = bootstrap(iterations, gammaWaitBase, gammaWaitData, shuffles, doWilcox, doSig);
    n = length(find(abs(cell2mat({rs_boot.zval}))>=abs(rs_stat.zval)));
    p_val = (n+1)/(iterations+1);
    currStruct.subj = 'allsubj';
    currStruct.clus = 'allclus2';
    currStruct.freq = freqRange{2};
    currStruct.stim = 'wait';
    currStruct.z_boot = {rs_boot.zval};
    currStruct.z_dat = rs_stat.zval;
    currStruct.p_value = p_val;
    summaryStruct = [summaryStruct;currStruct];
    if saveStats
        cd(statDir);
        save(['allClus2_' freqRange{2} '_waitBoot.mat'], 'p_val','rs_boot', 'rs_stat','shuffles')
    end
    
    thetaMoveDat = mean(allclus2_mNormThetaMove,1)';
    thetaWaitDat = mean(allclus2_mNormThetaWait,1)';
    thetaMoveBase = thetaMoveDat(baseTimeInd1:baseTimeInd2);
    thetaWaitBase = thetaWaitDat(baseTimeInd1:baseTimeInd2);
    thetaMoveData = thetaMoveDat(baseTimeInd2:end-512);
    thetaWaitData = thetaWaitDat(baseTimeInd2:end-512);
    disp(['allClus2 ' freqRange{1} ' move']);
    %[p_boot,p_sig,shuffles] = bootstrap(iterations, thetaMoveBase, thetaMoveData, shuffles, doWilcox, doSig);
    [rs_boot,rs_stat,shuffles] = bootstrap(iterations, thetaMoveBase, thetaMoveData, shuffles, doWilcox, doSig);
    n = length(find(abs(cell2mat({rs_boot.zval}))>=abs(rs_stat.zval)));
    p_val = (n+1)/(iterations+1);
    currStruct.subj = 'allsubj';
    currStruct.clus = 'allclus2';
    currStruct.freq = freqRange{1};
    currStruct.stim = 'move';
    currStruct.z_boot = {rs_boot.zval};
    currStruct.z_dat = rs_stat.zval;
    currStruct.p_value = p_val;
    summaryStruct = [summaryStruct;currStruct];
    if saveStats
        cd(statDir);
        save(['allClus2_' freqRange{1} '_moveBoot.mat'], 'p_val','rs_boot', 'rs_stat','shuffles')
    end
    disp(['allClus2 ' freqRange{1} ' wait']);
    %[p_boot,p_sig,shuffles] = bootstrap(iterations, thetaWaitBase, thetaWaitData, shuffles, doWilcox, doSig);
    [rs_boot,rs_stat,shuffles] = bootstrap(iterations, thetaWaitBase, thetaWaitData, shuffles, doWilcox, doSig);
    n = length(find(abs(cell2mat({rs_boot.zval}))>=abs(rs_stat.zval)));
    p_val = (n+1)/(iterations+1);
    currStruct.subj = 'allsubj';
    currStruct.clus = 'allclus2';
    currStruct.freq = freqRange{1};
    currStruct.stim = 'wait';
    currStruct.z_boot = {rs_boot.zval};
    currStruct.z_dat = rs_stat.zval;
    currStruct.p_value = p_val;
    summaryStruct = [summaryStruct;currStruct];
    if saveStats
        cd(statDir);
        save(['allClus2_' freqRange{1} '_waitBoot.mat'], 'p_val','rs_boot', 'rs_stat','shuffles')
    end
    
    disp('Done')
    clearvars -except clus1Struct clus2Struct useBipolar clus_info subjVec allClus freqRangeList pVec summaryStruct
end