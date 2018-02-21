function oddball_erp
% This function plots ERPs from the oddball task
%Inputs:
subj = 'HUP142_i';
saveFigs = 1;
elec = 68;
targEvents = 'TARGETHF';
backEvents = 'BACKGROUNDHF'; %Only used if subtractBackground is selected
alignment = 'startSound';
signalStartMS = -400;
signalDuration = 1200;

dirs = le_dirs('oddball');
saveFigsDir = fullfile(dirs.scratch,'figs',subj,'erp');
lfpDir = fullfile(dirs.data,'eeg',subj,'lfp.noreref');
evFile = fullfile(dirs.data,'eeg',subj,'behavioral','session_0/events.mat');

%Settings
zScoreRandBase = 0; % Plots z-score voltages using random baseline
subtractBackground = 0; %Subtract background tone signal from target

load(evFile,'events');
comparison = [targEvents backEvents];
[retain1,retain2] = le_events2retain(events,'oddball',alignment,comparison);

for i=1:length(events)
    if isempty(events(i).response)
        events(i).response = NaN;
    end
end
goCorrect = strcmp({events.type},'GO') & cell2mat({events.response})==1;

lfp1 = an_getlfp_ms_wrapper(elec,events(retain1),signalDuration,signalStartMS);
lfp2 = an_getlfp_ms_wrapper(elec,events(retain2),signalDuration,signalStartMS);
t = signalStartMS:1:signalStartMS+signalDuration-1;

if zScoreRandBase
    baseEv=get_baseline_random(events,3,1,[],'lfp');
    baseLfp = an_getlfp_ms_wrapper(elec,baseEv,signalDuration,signalStartMS,100);
    baseMean = nanmean(baseLfp(:,100:end-101),1);
    baseStd = std(baseLfp(:,100:end-101),1);
    zlfp1 = (lfp1-baseMean)./baseStd;
    zlfp2 = (lfp2-baseMean)./baseStd;
    
    if subtractBackground
        pLfp = zlfp1-zlfp2;
        %figure;
        hold all;
        eplotit({pLfp},t,'test','normV',0)
        plot(t,pLfp);
        plot([0 0],ylim,'k')
        title([subj ' ' num2str(elec) ' Normalized, Subtracted']);
%         if saveFigs
%             if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
%             cd(saveFigsDir);
%             print(gcf,[elec '_NormSub'],'-dpng');close
%         end
    else
        figure;
        hold all
        plot(t,zlfp1);
        %plotSEM(t,zlfp1,std(zlfp1));
        plot(t,zlfp2);
        %plotSEM(t,zlfp2,std(zlfp2),[0.99 0.01 0.01]);
        plot([0 0],ylim,'k')
        title([subj ' ' num2str(elec) ' Normalized']);
        legend('Target','Background');
    end
    return
end
if subtractBackground
    figure;
    hold all
    pLfp = lfp1-lfp2;
    plot(t,pLfp);
    title([subj ' ' num2str(elec) ' Raw, Subtracted']);
else
    hold all
    p=eplotit({lfp1,lfp2},t,'test','normV',0);
    legend([p.h{1},p.h{2}],'Target','Background');
%     plot(t,lfp1);
%     plot(t,lfp2);
%     plot([0 0],ylim,'k')
%     title([subj ' ' num2str(elec) ' Raw']);
%     legend('Target','Background');
end