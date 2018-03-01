function fidel_entrain
% function fidel_entrain
%
%   DR 05/2015

% parameters
ddir = 'G:\'; % directory
subj = 'fidel'; % subject
ch = 11; % recording channel

% % MS3-4, EC11, sedated
% sdat = {'fidel_stim_DT1_052615',4;...% [tank,block]
%         'fidel_stim_DT1_052615',5;...
%         'fidel_stim_DT1_052615',6;...
%         'fidel_stim_DT1_052615',7;...
%         'fidel_stim_DT1_052615',8}; 

% % MS3-4, EC11, awake
% sdat = {'fidel_CCDT_DT1_052815',3;...% [tank,block]
%         'fidel_CCDT_DT1_052815',4;...
%         'fidel_CCDT_DT1_052815',5;...
%         'fidel_CCDT_DT1_052815',6;...
%         'fidel_CCDT_DT1_052815',7;...
%         'fidel_CCDT_DT1_052915',3;...
%         'fidel_CCDT_DT1_052915',5};    

% % MS2-3, EC11, awake
% sdat = {'fidel_CCDT_DT1_060315',5;
%         'fidel_CCDT_DT1_060315',2;%12
%         'fidel_CCDT_DT1_060315',6;%11
%         'fidel_CCDT_DT1_060315',10;
%         'fidel_CCDT_DT1_060315',7;
%         'fidel_CCDT_DT1_060315',4;
%         'fidel_CCDT_DT1_060315',8}; 
    
% MS1-2, EC11, awake
sdat = {%'fidel_CCDT_DT1_060915',8; % 30Hz
        'fidel_CCDT_DT1_060815',7; % 40Hz
        'fidel_CCDT_DT1_060915',7; % 50Hz
        'fidel_CCDT_DT1_060815',8; % 60Hz
        'fidel_CCDT_DT1_060915',6; % 70Hz
        'fidel_CCDT_DT1_060815',9};% 80Hz 
    
pout = [10 90]; % outlier percentiles

% statistics
f = zeros(size(sdat,1),1); % stim frequency
m = zeros(size(sdat,1),2); % mean [amp,freq]
v = zeros(size(sdat,1),2); % std [amp,freq]
for ii = 1:size(sdat,1)
    try load([ddir subj '\processed\' sdat{ii,1} '_' num2str(sdat{ii,2}) '.mat'],'-mat');
    catch disp('run ''ram_epstats'' first'); return; end
    f(ii) = sfreq;
    epstats = entrain_stats{ch};

    % remove outliers
    pct = prctile(epstats,pout); 
    amp = epstats(:,1);
    amp(amp<pct(1,1) | amp>pct(2,1)) = [];
    freq = epstats(:,2);
    freq(freq<pct(1,2) | freq>pct(2,2)) = [];
    
    % stats
    m(ii,:) = [mean(amp), mean(freq)];
    v(ii,:) = [std(amp), std(freq)];
end

% plot
figure('Name',subj,'NumberTitle','off','Units','normalized','Position',[1/3 1/8 1/3 3/4],'Color','w');
subplot(2,1,1), hold on;
errorbar(f,m(:,1),v(:,1),'ko','MarkerFaceColor','k');
axis square; set(gca,'Box','off','XLim',[20 100],'YLim',[0 300]);
ylabel('response amplitude (\muV)');
xlabel('stimulus frequency (Hz)');
subplot(2,1,2), hold on;
errorbar(f,m(:,2),v(:,2),'ko','MarkerFaceColor','k');
axis square; set(gca,'Box','off','XLim',[20 100],'YLim',[20 100]);
line(get(gca,'XLim'),get(gca,'YLim'),'Color','k');
ylabel('response frequency (Hz)');
xlabel('stimulus frequency (Hz)');