function oddball_erp
% This function plots ERPs from the oddball task

%Inputs:
subj = 'HUP150_i';
%Settings:
zScoreRandBase = 0; % Plots z-score voltages using random baseline
subtractBackground = 0; %Subtract background tone signal from target
averageElec = 0; %averages all electrodes
stackedPlot = 0; %create stacked plot using all electrodes
plotThreeD = 0; %creates 3D plot
heatmapping = 1; %creates heatmap plot
saveFigs = 0; %save plots in specified folder
docar = 0;
dirs = le_dirs('oddball');
saveFigsDir = fullfile(dirs.data,'scratch','figs',subj,'erp'); %folder to save plots
%144 badChannels = [26,27,28,33,37,29,30,31,32,86];
%142 badChannels = [52,91,92,74,94,95,97,98,99,21,22,24,30,31,32,33,43,44,46,47,117,118,119,120,121,122,123,124,125,126];
%142 badChannels = [52,91,92,74,94,95,97,98,99,21,22,24,30,31,32,33,43,44,46,47];
%149 badChannels = [64,24];
%147 badChannels = [46,30,6];
badChannels = [];

targEvents = 'TARGETHF';
backEvents = 'BACKGROUNDHF'; %Only used if subtractBackground is selected
alignment = 'startSound';
signalStartMS = -400;
signalDuration = 1350;

lfpDir = fullfile(dirs.data,'eeg',subj,'lfp.noreref');
evFile = fullfile(dirs.data,'eeg',subj,'behavioral','session_0/events.mat');

load(evFile,'events');
comparison = [targEvents backEvents];
[retain1,retain2] = le_events2retain(events,'oddball',alignment,comparison);

for i=1:length(events)
    if isempty(events(i).response)
        events(i).response = NaN;
    end
end
goCorrect = strcmp({events.type},'GO') & cell2mat({events.response})==1;

%Read electrodeMap
xlfile = fullfile(dirs.data,'eeg',subj,'docs/electrodeMap');
[~,~,xlcells] = xlsread(xlfile);
emptyChannels = [];
lbls = {};
for i=1:length(xlcells)
    if isnan(xlcells{i,3})
        emptyChannels = [emptyChannels, i];
    elseif i<129
        lbls = [lbls;{i}, xlcells{i,3}];
    end
end
channels = cell2mat(lbls(:,1));
%For averages using all electrodes
if averageElec
    lfp1mat = NaN(length(channels),signalDuration);
    lfp2mat = NaN(length(channels),signalDuration);
    for elec = 1:length(channels)
        if ismember(channels(elec),badChannels),continue;end
        lfp1mat((elec),:) = nanmean(an_getlfp_ms_wrapper(channels(elec),events(retain1),signalDuration,signalStartMS),1);
        %noise = mtmlinenoise(lfp1mat(elec,:),3,1000,1000,60:60:300);
        %lfp1mat(elec,:) = lfp1mat(elec,:)-noise;
        lfp2mat((elec),:) = nanmean(an_getlfp_ms_wrapper(channels(elec),events(retain2),signalDuration,signalStartMS),1);
        %noise = mtmlinenoise(lfp2mat(elec,:),3,1000,1000,60:60:300);
        %lfp2mat(elec,:) = lfp2mat(elec,:)-noise;
    end
    avglfp1 = nanmean(lfp1mat((1:length(channels)),:));
    avglfp2 = nanmean(lfp2mat((1:length(channels)),:));
    t = signalStartMS:1:signalStartMS+signalDuration-1;
    hold all;
    plot(t,avglfp1);
    plot(t,avglfp2);
    if saveFigs
        if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
        cd(saveFigsDir);
        print(gcf,['Averages'],'-dpng');
        close;
    end
    return
end
%for stacked plot
if stackedPlot
    pLfpmat = NaN(length(channels),signalDuration);
    for elec = 1:length(channels)
    if ismember(channels(elec),badChannels),continue;end
    lfp1 = nanmean(an_getlfp_ms_wrapper(channels(elec),events(retain1),signalDuration,signalStartMS),1);
    lfp2 = nanmean(an_getlfp_ms_wrapper(channels(elec),events(retain2),signalDuration,signalStartMS),1);
    t = signalStartMS:1:signalStartMS+signalDuration-1;
    label = char(lbls(channels(elec),2));
    baseEv=get_baseline_random(events,3,1,[],'lfp');
    baseLfp = an_getlfp_ms_wrapper(channels(elec),baseEv,signalDuration,signalStartMS,100);
    baseMean = nanmean(baseLfp(:,100:end-101),1);
    baseStd = std(baseLfp(:,100:end-101),1);
    zlfp1 = (lfp1-baseMean)./baseStd;
    zlfp2 = (lfp2-baseMean)./baseStd;
    pLfpmat((elec),:) = (zlfp1-zlfp2)+channels(elec);
    end
    hold all;
    plot(t,(pLfpmat((1:length(channels)),:)));
    plot([0 0],ylim,'k');
    title([subj ' ' label ' Normalized, Subtracted'],'Interpreter','none');
    if saveFigs
        if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
        cd(saveFigsDir);
        print(gcf,['Stacked'],'-dpng');
        close;
    end
end
%for 3D plot
if plotThreeD
    pLfpmat = NaN(length(channels),signalDuration);
    for elec = 1:length(channels)
    if ismember(channels(elec),badChannels),continue;end
    lfp1 = nanmean(an_getlfp_ms_wrapper(channels(elec),events(retain1),signalDuration,signalStartMS),1);
    lfp2 = nanmean(an_getlfp_ms_wrapper(channels(elec),events(retain2),signalDuration,signalStartMS),1);
    t = signalStartMS:1:signalStartMS+signalDuration-1;
    label = char(lbls(channels(elec),2));
    baseEv=get_baseline_random(events,3,1,[],'lfp');
    baseLfp = an_getlfp_ms_wrapper(channels(elec),baseEv,signalDuration,signalStartMS,100);
    baseMean = nanmean(baseLfp(:,100:end-101),1);
    baseStd = std(baseLfp(:,100:end-101),1);
    zlfp1 = (lfp1-baseMean)./baseStd;
    zlfp2 = (lfp2-baseMean)./baseStd;
    pLfpmat((elec),:) = (zlfp1-zlfp2);
    end
    h=waterfall(t,channels,pLfpmat);
    set( h, 'LineWidth', 4 );
    hidden off;
    rotate3d;
    xlabel('Time');
    ylabel('Channels');
    zlabel('Z-score');
    return
end
%for heatmap plot
if heatmapping
    if docar
        lfp1mat=NaN(length(channels),signalDuration);
        lfp2mat=NaN(length(channels),signalDuration);
        t = signalStartMS:1:signalStartMS+signalDuration-1;
        for elec = 1:length(channels)
            if ismember(channels(elec),badChannels),continue;end
            lfp1 = nanmean(an_getlfp_ms_wrapper(channels(elec),events(retain1),signalDuration,signalStartMS),1);
            lfp2 = nanmean(an_getlfp_ms_wrapper(channels(elec),events(retain2),signalDuration,signalStartMS),1);
            lfp1mat((elec),:) = (lfp1);  
            lfp2mat((elec),:) = (lfp2);
        end
        if zScoreRandBase
            lfpmat = lfp1mat-lfp2mat;
            car = nanmean(lfpmat);
            pLfpmat = lfpmat-car;
            zLfpmat = (pLfpmat)./std(car);
        h=waterfall(t,channels,pLfpmat);
        set( h, 'LineWidth', 4 );
        hidden off;
        view([0, 90]);
        xlabel('Time');
        set(gca,'ytick',[1:126],'yticklabel',lbls(1:126,2),'FontSize',4)
        ylim([0 length(channels)])
            return
        end
        carLfp1=nanmean(lfp1mat);
        carLfp2=nanmean(lfp2mat);
        rerefLfp1mat=lfp1mat-carLfp1;
        rerefLfp2mat=lfp2mat-carLfp2;
        pLfpmat=rerefLfp1mat-rerefLfp2mat;
        h=waterfall(t,channels,pLfpmat);
        set( h, 'LineWidth', 4 );
        hidden off;
        view([0, 90]);
        xlabel('Time');
        set(gca,'ytick',[1:126],'yticklabel',lbls(1:126,2),'FontSize',4)
        ylim([0 length(channels)])
        return
    end
    pLfpmat = NaN(length(channels),signalDuration);
    for elec = 1:length(channels)
    if ismember(channels(elec),badChannels),continue;end
    lfp1 = nanmean(an_getlfp_ms_wrapper(channels(elec),events(retain1),signalDuration,signalStartMS),1);
    lfp2 = nanmean(an_getlfp_ms_wrapper(channels(elec),events(retain2),signalDuration,signalStartMS),1);
    t = signalStartMS:1:signalStartMS+signalDuration-1;
    label = char(lbls(channels(elec),2));
    baseEv=get_baseline_random(events,3,1,[],'lfp');
    baseLfp = an_getlfp_ms_wrapper(channels(elec),baseEv,signalDuration,signalStartMS,100);
    baseMean = nanmean(baseLfp(:,100:end-101),1);
    baseStd = std(baseLfp(:,100:end-101),1);
    zlfp1 = (lfp1-baseMean)./baseStd;
    zlfp2 = (lfp2-baseMean)./baseStd;
    pLfpmat((elec),:) = (zlfp1-zlfp2);
    end
    h=waterfall(t,channels,pLfpmat);
    set( h, 'LineWidth', 4 );
    hidden off;
    view([0, 90]);
    xlabel('Time');
    set(gca,'ytick',[1:126],'yticklabel',lbls(1:126,2),'FontSize',4)
    ylim([0 length(channels)])
    return
end
%for individual electrodes
for elec = 1:length(channels)
    if ismember(channels(elec),badChannels),continue;end
    lfp1 = nanmean(an_getlfp_ms_wrapper(channels(elec),events(retain1),signalDuration,signalStartMS),1);
    lfp2 = nanmean(an_getlfp_ms_wrapper(channels(elec),events(retain2),signalDuration,signalStartMS),1);
    t = signalStartMS:1:signalStartMS+signalDuration-1;
    label = char(lbls(channels(elec),2));
    if zScoreRandBase
        baseEv=get_baseline_random(events,3,1,[],'lfp');
        baseLfp = an_getlfp_ms_wrapper(channels(elec),baseEv,signalDuration,signalStartMS,100);
        baseMean = nanmean(baseLfp(:,100:end-101),1);
        baseStd = std(baseLfp(:,100:end-101),1);
        zlfp1 = (lfp1-baseMean)./baseStd;
        zlfp2 = (lfp2-baseMean)./baseStd;
        
        if subtractBackground
            pLfp = zlfp1-zlfp2;
            hold all;
%             eval(sprintf('%s = plot(t,pLfp);',label));
            plot(t,pLfp);
            plot([0 0],ylim,'k');
            title([subj ' ' label ' Normalized, Subtracted'],'Interpreter','none');
            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                  print(gcf,['Norm_Sub_' label],'-dpng');
                  close;
            end
        else
            hold all
            plot(t,zlfp1);
            plot(t,zlfp2);
            plot([0 0],ylim,'k');
            title([subj ' ' label ' Normalized'],'Interpreter','none');
            legend('Target','Background');
            if saveFigs
                if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
                cd(saveFigsDir);
                print(gcf,['Norm_' label],'-dpng');
                close;
            end
        end
        continue
    end
    if subtractBackground
        pLfp = lfp1-lfp2;      
         plot(t,pLfp);
         title([subj ' ' label ' Raw, Subtracted'],'Interpreter','none');
        if saveFigs
            if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
            cd(saveFigsDir);
            print(gcf,['Sub_' label],'-dpng');
            close;
        end
    else
        lfp1 = an_getlfp_ms_wrapper(channels(elec),events(retain1),signalDuration,signalStartMS);
        lfp2 = an_getlfp_ms_wrapper(channels(elec),events(retain2),signalDuration,signalStartMS);
        p=eplotit({lfp1 lfp2},t,label,'normV',0);
%         hold all;
%         plot(t,lfp1);
%         plot(t,lfp2);
         title([subj ' CSC' num2str(channels(elec)) ' ' label ' Raw'],'Interpreter','none');
         legend([p.h{1},p.h{2}],'Target','Background');
        if saveFigs
            if ~exist(saveFigsDir,'dir'),mkdir(saveFigsDir);end
            cd(saveFigsDir);
            print(gcf,[label],'-dpng');
            close;
        end
    end
end
keyboard;