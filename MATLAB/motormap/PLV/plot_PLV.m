
if ~exist('clus_info','var'),defineClusters_k2;end
freqRangeList = [{'delta600','theta250'};{'alpha200','beta100'};...
    {'lowFreq250','highFreq50'};{'broadtheta250','broadtheta300'}];
for jj=1:size(freqRangeList,1)
    dirs = le_dirs;
    freqRange = freqRangeList(jj,:)
    plvDir = fullfile(dirs.scratch,'PHASE/k2/',[freqRange{1} '_' freqRange{2}]);
    saveFigs = 1;
    figDir = fullfile(dirs.scratch,'figs/motormap/PLV',[freqRange{1} '_' freqRange{2}]);
    if ~exist(figDir,'dir'),mkdir(figDir);end
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
                %plot
                hold on
                title([subj ' - ' clus ' ' freqRange{2}])
                plot(t,fgsmooth(mean(motor_GammaMove,1),10),'g')
                plot(t,fgsmooth(mean(motor_GammaWait,1),10),'r')
                xlim([-2 6.25])
                ylim([0 1])
                plot([0 0], [-7 7], 'b-', 'linewidth', 2)
                plot([-3 7], [0 0], 'k--', 'linewidth', 2)
                plot([0 5], [-1 -1], 'k-', 'linewidth' ,2)
                xlabel('Time (s)')
                ylabel('Raw PLV')
                legend Move Wait
                if saveFigs
                    cd(figDir);
                    print(gcf,[subj ' - ' clus ' ' freqRange{2}],'-dpng');close
                end
                %             if ii ==1
                %                 allclus1_mNormGammaMove = [allclus1_mNormGammaMove; motor_GammaMove];
                %                 allclus1_mNormGammaWait = [allclus1_mNormGammaWait; motor_GammaWait];
                %             else
                %                 allclus2_mNormGammaMove = [allclus2_mNormGammaMove; motor_GammaMove];
                %                 allclus2_mNormGammaWait = [allclus2_mNormGammaWait; motor_GammaWait];
                %             end
            end
            
            if exist('motor_normThetaMove','var')
                hold on
                title([subj ' - ' clus ' ' freqRange{1}])
                plot(t,fgsmooth(mean(motor_ThetaMove,1),10),'g')
                plot(t,fgsmooth(mean(motor_ThetaWait,1),10),'r')
                xlim([-2 6.25])
                ylim([0 1])
                plot([0 0], [-7 7], 'b-', 'linewidth', 2)
                plot([-3 7], [0 0], 'k--', 'linewidth', 2)
                plot([0 5], [-1 -1], 'k-', 'linewidth' ,2)
                xlabel('Time (s)')
                ylabel('Raw PLV')
                legend Move Wait
                if saveFigs
                    cd(figDir);
                    print(gcf,[subj ' - ' clus ' ' freqRange{1}],'-dpng');close
                end
                %             if ii ==1
                %                 allclus1_mNormThetaMove = [allclus1_mNormThetaMove; motor_ThetaMove];
                %                 allclus1_mNormThetaWait = [allclus1_mNormThetaWait; motor_ThetaWait];
                %             else
                %                 allclus2_mNormThetaMove = [allclus2_mNormThetaMove; motor_ThetaMove];
                %                 allclus2_mNormThetaWait = [allclus2_mNormGammaWait; motor_ThetaWait];
                %             end
            end
            clearvars -except clus1Struct clus2Struct useBipolar clus_info subjVec allClus i ii plvDir subjList saveFigs figDir allclus* t freqRange freqRangeList
        end
    end
    
    % %plot
    % hold on
    % title(['all clus1 ' freqRange{2}])
    % plot(t,fgsmooth(mean(allclus1_mNormGammaMove),10),'g')
    % plot(t,fgsmooth(mean(allclus1_mNormGammaWait),10),'r')
    % xlim([-2 6.25])
    % ylim([0 1])
    % plot([0 0], [-7 7], 'b-', 'linewidth', 2)
    % plot([-3 7], [0 0], 'k--', 'linewidth', 2)
    % plot([0 5], [-1 -1], 'k-', 'linewidth' ,2)
    % xlabel('Time (s)')
    % ylabel('Change in PLV')
    % legend Move Wait
    % if saveFigs
    %     cd(figDir);
    %     print(gcf,['all_clus1 ' freqRange{2}],'-dpng');close
    % end
    %
    % hold on
    % title(['all clus1 ' freqRange{1}])
    % plot(t,fgsmooth(mean(allclus1_mNormThetaMove),10),'g')
    % plot(t,fgsmooth(mean(allclus1_mNormThetaWait),10),'r')
    % xlim([-2 6.25])
    % ylim([0 1])
    % plot([0 0], [-7 7], 'b-', 'linewidth', 2)
    % plot([-3 7], [0 0], 'k--', 'linewidth', 2)
    % plot([0 5], [-1 -1], 'k-', 'linewidth' ,2)
    % xlabel('Time (s)')
    % ylabel('Change in PLV')
    % legend Move Wait
    % if saveFigs
    %     cd(figDir);
    %     print(gcf,['all_clus1 ' freqRange{1}],'-dpng');close
    % end
    %
    % hold on
    % title(['all clus2 ' freqRange{2}])
    % plot(t,fgsmooth(mean(allclus2_mNormGammaMove),10),'g')
    % plot(t,fgsmooth(mean(allclus2_mNormGammaWait),10),'r')
    % xlim([-2 6.25])
    % ylim([0 1])
    % plot([0 0], [-7 7], 'b-', 'linewidth', 2)
    % plot([-3 7], [0 0], 'k--', 'linewidth', 2)
    % plot([0 5], [-1 -1], 'k-', 'linewidth' ,2)
    % xlabel('Time (s)')
    % ylabel('Change in PLV')
    % legend Move Wait
    % if saveFigs
    %     cd(figDir);
    %     print(gcf,['all_clus2 ' freqRange{2}],'-dpng');close
    % end
    %
    % hold on
    % title(['all clus2 ' freqRange{1}])
    % plot(t,fgsmooth(mean(allclus2_mNormThetaMove),10),'g')
    % plot(t,fgsmooth(mean(allclus2_mNormThetaWait),10),'r')
    % xlim([-2 6.25])
    % ylim([0 1])
    % plot([0 0], [-7 7], 'b-', 'linewidth', 2)
    % plot([-3 7], [0 0], 'k--', 'linewidth', 2)
    % plot([0 5], [-1 -1], 'k-', 'linewidth' ,2)
    % xlabel('Time (s)')
    % ylabel('Change in PLV')
    % legend Move Wait
    % if saveFigs
    %     cd(figDir);
    %     print(gcf,['all_clus2 ' freqRange{1}],'-dpng');close
    % end
    
    clearvars -except clus1Struct clus2Struct useBipolar clus_info subjVec allClus freqRangeList
end
disp('done')