function ditka_entrain
% function ditka_entrain
%
%   DR 04/2015

% parameters
ddir = 'D:\'; % directory
subj = 'ditka'; % subject
ch = 12; % recording channel

% MS4-5
sdat = {'ditka_stim_DT1_033015',4; % 50Hz
        'ditka_stim_DT1_033015',7; % 40Hz
        'ditka_stim_DT1_033015',8; % 80Hz
        'ditka_stim_DT1_040215',2; % 30Hz
        'ditka_stim_DT1_040215',3; % 60Hz
        'ditka_stim_DT1_040615',6};% 70Hz
    
% EC4-5, HC6
% sdat = {'ditka_stim_DT1_040215',11;...% [tank,block]
%         'ditka_stim_DT1_040615',4;...
%         'ditka_stim_DT1_040215',10;...
%         'ditka_stim_DT1_040615',2;...
%         'ditka_stim_DT1_040215',12;...
%         'ditka_stim_DT1_040815',7;...
%         'ditka_stim_DT1_040815',8}; 

% % MS3-4, HC6
% sdat = {'ditka_stim_DT1_041315',4;...% [tank,block]
%         'ditka_stim_DT1_041315',5;...
%         'ditka_stim_DT1_041315',6;...
%         'ditka_stim_DT1_041315',7;...
%         'ditka_stim_DT1_041315',8;...
%         'ditka_stim_DT1_041315',9};     
    
% % MS5-6, HC6
% sdat = {'ditka_stim_DT1_041515',2;...% [tank,block]
%         'ditka_stim_DT1_041515',3;...
%         'ditka_stim_DT1_041515',8;...
%         'ditka_stim_DT1_041515',5;...
%         'ditka_stim_DT1_041515',6;...
%         'ditka_stim_DT1_041515',7}; 

% MS2-3, HC6
% sdat = {'ditka_stim_DT1_050115',2;...% [tank,block]
%         'ditka_stim_DT1_050115',3;...
%         'ditka_stim_DT1_050115',5}; 
    
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
axis square; set(gca,'Box','off','XLim',[10 90],'YLim',[0 700]);
ylabel('response amplitude (\muV)');
xlabel('stimulus frequency (Hz)');
subplot(2,1,2), hold on;
errorbar(f,m(:,2),v(:,2),'ko','MarkerFaceColor','k');
axis square; set(gca,'Box','off','XLim',[10 90],'YLim',[10 90]);
line(get(gca,'XLim'),get(gca,'YLim'),'Color','k');
ylabel('response frequency (Hz)');
xlabel('stimulus frequency (Hz)');