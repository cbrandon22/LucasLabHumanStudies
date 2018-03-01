function ram_sta_cage(varargin)
% function ram_sta_cage(varargin)
%   Stimulus-triggered average across trials - cage experiment.
%
%   DR 05/2015

% parameters
ddir = 'F:\'; % directory
subj = 'ditka'; % subject
tank = 'ditka_cage_051815'; % tank 

% data
fs = 10000; % sampling rate (Hz)
scale = 3.3/2^12/1000*1e6; % amplitude scale factor to convert to uV
ISIr = [28 32]; % inter-stimulus interval range (s)
dur = 5; % train duration (ms)
thresh = 200; % stimulus artifact threshold
twin = [-100 200]; % peri-stimulus window (ms)
debg = 1; % debug plots? (0 or 1)
try 
    load([ddir subj filesep 'processed' filesep tank '.mat'],'-mat');
    size(STA); size(trainstat); size(tsta); size(t);
catch
    cd([ddir subj '\' tank]);
    Nfl = length(dir('*_*.txt'));
    DAT = [];
    for ii = 0:Nfl-1
        fnm = [tank(end-5:end-2) '_' num2str(ii) '.txt'];
        fid = fopen(fnm,'r');
        DAT = [DAT; fread(fid,inf,'uint16')];
        fclose(fid);
    end

    % stimulus trains
    [b,a] = butter(2,300/(fs/2),'high');
    fDAT = filtfilt(b,a,DAT);
    ind = find(fDAT(1:end-1)<thresh & fDAT(2:end)>thresh);
    ind2 = find(diff(ind)>ISIr(1)*fs & diff(ind)<ISIr(2)*fs);
    itrig = ind(ind2+1)-2; % indices of train onset
    trainstat = zeros(length(itrig),2); % train number of pulses and mean frequency
    for ii = 1:length(itrig)
        iii = find(ind>=itrig(ii) & ind<=itrig(ii)+dur/1000*fs);
        trainstat(ii,1) = length(iii);
        trainstat(ii,2) = fs/mean(diff(ind(iii)));
    end
    
    % stim-triggered sweeps
    swin = round(twin(1)/1000*fs):round(twin(2)/1000*fs);
    t = swin/fs*1000;
    STA = zeros(length(itrig),length(t));
    for ii = 1:length(itrig)
        STA(ii,:) = DAT(itrig(ii)+swin)';
        if debg
            figure('Name','sta cage','NumberTitle','off','Units','normalized','Position',[1/3 1/4 1/3 1/2]);
            plot(t,STA(ii,:),'k');
            axis tight; axis square; set(gca,'Box','off'); xlabel('ms'); title([num2str(trainstat(ii,1)) ' pulses @ ' num2str(trainstat(ii,2)) ' Hz']);
            pause;
            close;
        end
    end
    tsta = itrig/fs/3600; % train onset (hours)
    
    % save
    save([ddir subj filesep 'processed' filesep tank],'STA','tsta','trainstat','t','-append');
end

% pre-plot
STA = scale*(STA-median(STA(:,t<0),2)*ones(1,length(t))); % center and scale
% ind = find(tsta>1); % restrict to certain hours of recording
% ind = find(tsta>0.2 & tsta<0.65);
% trainstat = trainstat(ind,:);
% STA = STA(ind,:);

% plot
sf = unique(trainstat(:,1));
for ii = 1:length(sf)
    csta = STA(trainstat(:,1)==sf(ii),:);
    chz = trainstat(trainstat(:,1)==sf(ii),2);
    if size(csta,1)<5, continue; end
    figure('Name','sta cage','NumberTitle','off','Units','normalized','Position',[1/3 1/4 1/3 1/2]);
    plot(t,mean(csta),'k'); 
    axis tight; axis square; set(gca,'Box','off'); xlabel('ms'); ylabel('\muV'); title([num2str(sf(ii)) ' pulses @ ' num2str(mean(chz),'%4.1f') ' Hz']);
    text(max(get(gca,'XLim')),max(get(gca,'YLim')),['N = ' num2str(size(csta,1))],'HorizontalAlignment','right','VerticalAlignment','top');
end
