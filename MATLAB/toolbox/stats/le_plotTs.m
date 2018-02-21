function le_plotTs(tStruct,config_pow,task,printFlag,printFolderLbl)
% This function takes a tStruct and power configuration and creates
% a normalized power plot for each electrode. If tStruct contains
% information for multiple electrodes, it will loop through those
% electrodes when plotting. KEY POINT: It makes differents plots based 
% on the dimensionality of tMat and powMats. If 2d data is given, it will
% plotSEM row vectors showing mean +/- SEM. If 3d data is given, it will 
% make a TF plot.


% Inputs:
% tStruct       struct created by le_calcTs.m
% config_pow    power configuration (also output from le_calcTs.m)
% printFlag     if set to 1, it will print the figures to disk; default = 0
% printFolderLbl %optional - string describing elec group that will print
                 % figs in a separate folder


if ~exist('printFlag','var') || isempty(printFlag)
    printFlag = 0;
end
if ~exist('printFolderLbl','var') || isempty(printFolderLbl)
    printFolderLbl = '';
end

% set dirs
dirs = le_dirs(task);

% deduce collapseFreqFlag
collapseFreqFlag = size(tStruct(1).tMat,1) == 1;

% load anatStruct (must go to subject's tal dir to get that anatStruct)
if strcmp(tStruct(1).task, 'motormap')
    [anatStruct] = le_centroids2Anatstruct(tStruct(1).subj);
elseif strcmp(tStruct(1).task, 'CCDT')
    jacFile = fullfile(dirs.data,'eeg',tStruct(1).task,tStruct(1).subj,'docs','jacksheet.txt');
    fid = fopen(jacFile,'r');
    JAC=textscan(fid,'%d%s','delimiter','\t');
    JAC{:,2} = strtrim(JAC{:,2});
    fclose(fid);
    channels = JAC{1};
    eLbls = JAC{2};
    excludeChan = {'EKG1','EKG2'};
    [lia,locb] = ismember(excludeChan,JAC{2});
    if sum(lia)>0
        channels(locb) = [];
        eLbls(locb) = [];
    end
    anatStruct.eLbl = eLbls;
end
    

% get values for x-axis (time bins) and yaxis (freq bins, accounting for ...
%custom fRange used to generate this tStruct)
tBins = nanmean(config_pow.timeBins,2)';
fBins = [config_pow.freqBins];
fBins = fBins(tStruct(1).fInd); % filter by freq range used for this tStruct



% make both plots
for i=1:length(tStruct)
    
    if collapseFreqFlag
        % make figure here
        swagFig

        plot(tBins,nanmean(tStruct(i).pow1,1),'r');
        hold all
        plot(tBins,nanmean(tStruct(i).pow2,1),'b');
        tit = ['ch ' tStruct(i).eLbl];
        %tit = anatStruct(~cellfun(@isempty,strfind({anatStruct.eLbl},tStruct(i).eLbl))).elecLbl;
        swagAxes(gca,20,'Time from GO cue (ms)','Normalized Power',tit)
        legend(tStruct(i).retLbl1,tStruct(i).retLbl2,'Location','northeastoutside')


        plotSEM(tBins,nanmean(tStruct(i).pow1,1),SEM(tStruct(i).pow1),'r');
        hold all
        plotSEM(tBins,nanmean(tStruct(i).pow2,1),SEM(tStruct(i).pow2),'b');
        plot(tBins,nanmean(tStruct(i).pow1,1),'--k','linewidth',2);
        plot(tBins,nanmean(tStruct(i).pow2,1),'--k','linewidth',2);

        % label plot and set preferences
        plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [0 0 0])

        
    else % power was not collapsed, so make a TF plot
        %
         swagFig([0.0109    0.5487    0.9109    0.3312])
        
        % Z POW 1
        ax(1) = subplot(1,3,1); hold all
        imagesc(tBins,1:length(fBins),nanmean(tStruct(i).pow1,3))
        set(gca,'ydir','normal','clim',[-1 1])
        setFreqYTicks(gca,fBins)
        swagAxes(gca,20,'Time from GO cue (ms)','Frequency (Hz)',['zPow ' tStruct(i).retLbl1])
        % plot line
        plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [1 1 1])
        c = colorbar('eastoutside');colormap('jet')
        
        % Z POW 2
        ax(2) = subplot(1,3,2); hold all
        imagesc(tBins,1:length(fBins),nanmean(tStruct(i).pow2,3))
        set(gca,'ydir','normal','clim',[-1 1])
        setFreqYTicks(gca,fBins)
        swagAxes(gca,20,'Time from GO cue (ms)','Frequency (Hz)',['zPow ' tStruct(i).retLbl2])
        % plot line
        plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [1 1 1])
        c = colorbar('eastoutside');colormap('jet')
        
        % t(Z POW 1 - zPow 2)
        ax(3) = subplot(1,3,3); hold all
        imagesc(tBins,1:length(fBins),nanmean(tStruct(i).tMat,3))
        set(gca,'ydir','normal','clim',[-5 3])
        setFreqYTicks(gca,fBins)
        swagAxes(gca,20,'Time from GO cue (ms)','Frequency (Hz)',['t(' tStruct(i).retLbl1 '-' tStruct(i).retLbl2 ')'])
        % plot line
        plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [1 1 1])
        c = colorbar('eastoutside');colormap('jet')
        
       
        
    end
    % parse printFlag
    if printFlag ==1
       if collapseFreqFlag
         freqLbl = [num2str(round(fBins(1))) '-to-' num2str(round(fBins(end)))];
       else
         freqLbl = ['TF-' num2str(round(fBins(1))) '-to-' num2str(round(fBins(end)))];
       end 
        
       printDir = fullfile(dirs.scratch,'figs',tStruct(1).task,...
           tStruct(1).subj,tStruct(1).sessLbl,[tStruct(1).comparison '-' printFolderLbl],freqLbl);
       cd_mkdir(printDir)
       print(gcf,tStruct(i).eLbl,'-dpng');close
    end
          
end


% Plot Y ticks for TF plots
function [] = setFreqYTicks(ax,fBins)
yt = [[1:5:length(fBins)] length(fBins)];
if any(diff(yt)==0)
    yt = [1:5:length(fBins)];
end
set(ax,'ytick',yt,...
'yticklabel',round(fBins(yt)))

