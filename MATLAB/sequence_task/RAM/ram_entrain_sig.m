function ram_entrain_sig
% function ram_entrain_sig
%
%   DR 10/2015

% parameters
ddir = 'D:\'; % directory
subj = 'ditka'; % subject
ch = 9:16; % channels
pout = [10 90]; % outlier percentiles

% DITKA MS4-5
sdat = {'ditka_stim_DT1_033015',4; % 50Hz
        'ditka_stim_DT1_033015',7; % 40Hz
        'ditka_stim_DT1_033015',8; % 80Hz
        'ditka_stim_DT1_040215',2; % 30Hz
        'ditka_stim_DT1_040215',3; % 60Hz
        'ditka_stim_DT1_040615',6};% 70Hz
    
% % FIDEL MS1-2
% sdat = {%'fidel_CCDT_DT1_060915',8; % 30Hz
%         'fidel_CCDT_DT1_060815',7; % 40Hz
%         'fidel_CCDT_DT1_060915',7; % 50Hz
%         'fidel_CCDT_DT1_060815',8; % 60Hz
%         'fidel_CCDT_DT1_060915',6; % 70Hz
%         'fidel_CCDT_DT1_060815',9};% 80Hz    
    
% statistics
f = zeros(size(sdat,1),1); % stim frequency
p = zeros(size(sdat,1),length(ch)); % pval
a = zeros(size(sdat,1),length(ch)); % amp
for ii = 1:size(sdat,1)
    try load([ddir subj '\processed\' sdat{ii,1} '_' num2str(sdat{ii,2}) '.mat'],'-mat');
    catch disp('run ''ram_epstats'' first'); return; end
    f(ii) = sfreq;
    for jj = 1:length(ch)
        epstats = entrain_stats{ch(jj)};
        
        % remove outliers
        pct = prctile(epstats,pout);
        amp = epstats(:,1);
        amp(amp<pct(1,1) | amp>pct(2,1)) = [];
        freq = epstats(:,2);
        freq(freq<pct(1,2) | freq>pct(2,2)) = [];
        
        % stats
        [~,p(ii,jj)] = ttest(freq,f(ii));
        a(ii,jj) = mean(amp);
    end
end
[f,indf] = sort(f);
p = p(indf,:);    
a = a(indf,:);

% plot
psig = 0.05/numel(p);
a = a/max(a(:));
figure('Name',subj,'NumberTitle','off','Units','normalized','Position',[1/3 1/4 1/3 1/2],'Color','w');
for ii = 1:size(p,1)
    for jj = 1:size(p,2)
        if p(ii,jj)>psig
            rectangle('Position',[jj-a(ii,jj)/2 ii-a(ii,jj)/2 a(ii,jj) a(ii,jj)],'Curvature',[1 1],'LineWidth',1,'FaceColor','k','EdgeColor','k');
        end
    end
end
set(gca,'Box','off','XLim',[0 size(p,2)+1],'XTick',1:size(p,2),'XTickLabel',ch,...
    'YLim',[0 size(p,1)+1],'YTick',1:size(p,1),'YTickLabel',f); axis equal;
xlabel('channel'); ylabel('frequency (Hz)');
