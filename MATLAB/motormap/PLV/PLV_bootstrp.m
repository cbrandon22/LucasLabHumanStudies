% Perform 2 sample bootstrap t-test on motormap PLV
if ~exist('clus_info','var'),defineClusters_k2;end
freqRangeList = [{'alpha200','beta125'};{'broadtheta250','broadtheta350'};...
        {'delta800','theta350'};{'lowFreq150','highFreq25'};{'lowFreq250','highFreq50'}];
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
                gammaMoveDat = mean(motor_normGammaMove,1)';
                gammaWaitDat = mean(motor_normGammaWait,1)';
                gammaMoveBase = gammaMoveDat(baseTimeInd1:baseTimeInd2);
                gammaWaitBase = gammaWaitDat(baseTimeInd1:baseTimeInd2);
                gammaMoveData = gammaMoveDat(baseTimeInd2:end);
                gammaWaitData = gammaWaitDat(baseTimeInd2:end);
                disp([subj ': ' clus ' ' freqRange{2} ' move']);
                m = bootstrp(iterations, @mean, gammaMoveBase, gammaMoveData);
                figure;
                [fi,xi] = ksdensity(m);
                plot(xi,fi);
                if saveStats
                    cd(statDir);
                    save([subj '_' clus '_' freqRange{2} '_moveBoot.mat'], 'p_boot', 'p_sig','shuffles')
                    if p_sig
                        cd(sigStatDir)
                        save([subj '_' clus '_' freqRange{2} '_moveBoot.mat'], 'p_boot', 'p_sig','shuffles')
                    end
                end
                disp([subj ': ' clus ' ' freqRange{2} ' wait']);
                [p_boot,p_sig,shuffles] = bootstrap(iterations, gammaWaitBase, gammaWaitData, shuffles, doWilcox, doSig);
                if saveStats
                    cd(statDir);
                    save([subj '_' clus '_' freqRange{2} '_waitBoot.mat'], 'p_boot', 'p_sig','shuffles')
                    if p_sig
                        cd(sigStatDir)
                        save([subj '_' clus '_' freqRange{2} '_waitBoot.mat'], 'p_boot', 'p_sig','shuffles')
                    end
                end
                if ii ==1
                    allclus1_mNormGammaMove = [allclus1_mNormGammaMove; motor_normGammaMove];
                    allclus1_mNormGammaWait = [allclus1_mNormGammaWait; motor_normGammaWait];
                else
                    allclus2_mNormGammaMove = [allclus2_mNormGammaMove; motor_normGammaMove];
                    allclus2_mNormGammaWait = [allclus2_mNormGammaWait; motor_normGammaWait];
                end
            end
            
            if exist('motor_normThetaMove','var')
                thetaMoveDat = mean(motor_normThetaMove,1)';
                thetaWaitDat = mean(motor_normThetaWait,1)';
                thetaMoveBase = thetaMoveDat(baseTimeInd1:baseTimeInd2);
                thetaWaitBase = thetaWaitDat(baseTimeInd1:baseTimeInd2);
                thetaMoveData = thetaMoveDat(baseTimeInd2:end);
                thetaWaitData = thetaWaitDat(baseTimeInd2:end);
                disp([subj ': ' clus ' ' freqRange{1} ' move']);
                [p_boot,p_sig,shuffles] = bootstrap(iterations, thetaMoveBase, thetaMoveData, shuffles, doWilcox, doSig);
                if saveStats
                    cd(statDir);
                    save([subj '_' clus '_' freqRange{1} '_moveBoot.mat'], 'p_boot', 'p_sig','shuffles')
                    if p_sig
                        cd(sigStatDir)
                        save([subj '_' clus '_' freqRange{1} '_moveBoot.mat'], 'p_boot', 'p_sig','shuffles')
                    end
                end
                disp([subj ': ' clus ' ' freqRange{1} ' wait']);
                [p_boot,p_sig,shuffles] = bootstrap(iterations, thetaWaitBase, thetaWaitData, shuffles, doWilcox, doSig);
                if saveStats
                    cd(statDir);
                    save([subj '_' clus '_' freqRange{1} '_waitBoot.mat'], 'p_boot', 'p_sig','shuffles')
                    if p_sig
                        cd(sigStatDir)
                        save([subj '_' clus '_' freqRange{1} '_waitBoot.mat'], 'p_boot', 'p_sig','shuffles')
                    end
                end
                if ii ==1
                    allclus1_mNormThetaMove = [allclus1_mNormThetaMove; motor_normThetaMove];
                    allclus1_mNormThetaWait = [allclus1_mNormThetaWait; motor_normThetaWait];
                else
                    allclus2_mNormThetaMove = [allclus2_mNormThetaMove; motor_normThetaMove];
                    allclus2_mNormThetaWait = [allclus2_mNormGammaWait; motor_normThetaWait];
                end
            end
            clearvars -except clus1Struct clus2Struct useBipolar clus_info subjVec allClus i ii plvDir subjList saveStats statDir sigStatDir allclus* t freqRange* base* iterations doWilcox doSig shuffles
        end
    end
    
    %All Cluster 1
    gammaMoveDat = mean(allclus1_mNormGammaMove,1)';
    gammaWaitDat = mean(allclus1_mNormGammaWait,1)';
    gammaMoveBase = gammaMoveDat(baseTimeInd1:baseTimeInd2);
    gammaWaitBase = gammaWaitDat(baseTimeInd1:baseTimeInd2);
    gammaMoveData = gammaMoveDat(baseTimeInd2:end);
    gammaWaitData = gammaWaitDat(baseTimeInd2:end);
    disp(['allClus1 ' freqRange{2} ' move']);
    [p_boot,p_sig,shuffles] = bootstrap(iterations, gammaMoveBase, gammaMoveData, shuffles, doWilcox, doSig);
    if saveStats
        cd(statDir);
        save(['allClus1_' freqRange{2} '_moveBoot.mat'], 'p_boot', 'p_sig','shuffles')
        if p_sig
            cd(sigStatDir)
            save(['allClus1_' freqRange{2} '_moveBoot.mat'], 'p_boot', 'p_sig','shuffles')
        end
    end
    disp(['allClus1 ' freqRange{2} ' wait']);
    [p_boot,p_sig,shuffles] = bootstrap(iterations, gammaWaitBase, gammaWaitData, shuffles, doWilcox, doSig);
    if saveStats
        cd(statDir);
        save(['allClus1_' freqRange{2} '_waitBoot.mat'], 'p_boot', 'p_sig','shuffles')
        if p_sig
            cd(sigStatDir)
            save(['allClus1_' freqRange{2} '_waitBoot.mat'], 'p_boot', 'p_sig','shuffles')
        end
    end
    
    thetaMoveDat = mean(allclus1_mNormThetaMove,1)';
    thetaWaitDat = mean(allclus1_mNormThetaWait,1)';
    thetaMoveBase = thetaMoveDat(baseTimeInd1:baseTimeInd2);
    thetaWaitBase = thetaWaitDat(baseTimeInd1:baseTimeInd2);
    thetaMoveData = thetaMoveDat(baseTimeInd2:end);
    thetaWaitData = thetaWaitDat(baseTimeInd2:end);
    disp(['allClus1 ' freqRange{1} ' move']);
    [p_boot,p_sig,shuffles] = bootstrap(iterations, thetaMoveBase, thetaMoveData, shuffles, doWilcox, doSig);
    if saveStats
        cd(statDir);
        save(['allClus1_' freqRange{1} '_moveBoot.mat'], 'p_boot', 'p_sig','shuffles')
        if p_sig
            cd(sigStatDir)
            save(['allClus1_' freqRange{1} '_moveBoot.mat'], 'p_boot', 'p_sig','shuffles')
        end
    end
    disp(['allClus1 ' freqRange{1} ' wait']);
    [p_boot,p_sig,shuffles] = bootstrap(iterations, thetaWaitBase, thetaWaitData, shuffles, doWilcox, doSig);
    if saveStats
        cd(statDir);
        save(['allClus1_' freqRange{1} '_waitBoot.mat'], 'p_boot', 'p_sig','shuffles')
        if p_sig
            cd(sigStatDir)
            save(['allClus1_' freqRange{1} '_waitBoot.mat'], 'p_boot', 'p_sig','shuffles')
        end
    end
    
    %All cluster 2
    gammaMoveDat = mean(allclus2_mNormGammaMove,1)';
    gammaWaitDat = mean(allclus2_mNormGammaWait,1)';
    gammaMoveBase = gammaMoveDat(baseTimeInd1:baseTimeInd2);
    gammaWaitBase = gammaWaitDat(baseTimeInd1:baseTimeInd2);
    gammaMoveData = gammaMoveDat(baseTimeInd2:end);
    gammaWaitData = gammaWaitDat(baseTimeInd2:end);
    disp(['allClus2 ' freqRange{2} ' move']);
    [p_boot,p_sig,shuffles] = bootstrap(iterations, gammaMoveBase, gammaMoveData, shuffles, doWilcox, doSig);
    if saveStats
        cd(statDir);
        save(['allClus2_' freqRange{2} '_moveBoot.mat'], 'p_boot', 'p_sig','shuffles')
        if p_sig
            cd(sigStatDir)
            save(['allClus2_' freqRange{2} '_moveBoot.mat'], 'p_boot', 'p_sig','shuffles')
        end
    end
    disp(['allClus2 ' freqRange{2} ' wait']);
    [p_boot,p_sig,shuffles] = bootstrap(iterations, gammaWaitBase, gammaWaitData, shuffles, doWilcox, doSig);
    if saveStats
        cd(statDir);
        save(['allClus2_' freqRange{2} '_waitBoot.mat'], 'p_boot', 'p_sig','shuffles')
        if p_sig
            cd(sigStatDir)
            save(['allClus2_' freqRange{2} '_waitBoot.mat'], 'p_boot', 'p_sig','shuffles')
        end
    end
    
    thetaMoveDat = mean(allclus2_mNormThetaMove,1)';
    thetaWaitDat = mean(allclus2_mNormThetaWait,1)';
    thetaMoveBase = thetaMoveDat(baseTimeInd1:baseTimeInd2);
    thetaWaitBase = thetaWaitDat(baseTimeInd1:baseTimeInd2);
    thetaMoveData = thetaMoveDat(baseTimeInd2:end);
    thetaWaitData = thetaWaitDat(baseTimeInd2:end);
    disp(['allClus2 ' freqRange{1} ' move']);
    [p_boot,p_sig,shuffles] = bootstrap(iterations, thetaMoveBase, thetaMoveData, shuffles, doWilcox, doSig);
    if saveStats
        cd(statDir);
        save(['allClus2_' freqRange{1} '_moveBoot.mat'], 'p_boot', 'p_sig','shuffles')
        if p_sig
            cd(sigStatDir)
            save(['allClus2_' freqRange{1} '_moveBoot.mat'], 'p_boot', 'p_sig','shuffles')
        end
    end
    disp(['allClus2 ' freqRange{1} ' wait']);
    [p_boot,p_sig,shuffles] = bootstrap(iterations, thetaWaitBase, thetaWaitData, shuffles, doWilcox, doSig);
    if saveStats
        cd(statDir);
        save(['allClus2_' freqRange{1} '_waitBoot.mat'], 'p_boot', 'p_sig','shuffles')
        if p_sig
            cd(sigStatDir)
            save(['allClus2_' freqRange{1} '_waitBoot.mat'], 'p_boot', 'p_sig','shuffles')
        end
    end
    
    disp('Done')
    clearvars -except clus1Struct clus2Struct useBipolar clus_info subjVec allClus freqRangeList
end