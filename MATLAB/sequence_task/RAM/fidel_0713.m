function fidel_0713
% function fidel_0713

% parameters
p.ddir = 'E:\'; % directory
p.subj = 'fidel'; % subject
p.chrec = 3; % rec channel
p.band = [8 12]; % frequency band (Hz)
p.dns = 10; % downsample factor
p.rln = 0; % remove line noise? (0 or 1)
tkbk = {'fidel_CCDT_DT1_061615',3;  % 0 Hz
        'fidel_CCDT_DT1_061715',4;  % 6 Hz
        'fidel_CCDT_DT1_061615',11; % 8 Hz       
        'fidel_CCDT_DT1_061615',12; % 10 Hz
        'fidel_CCDT_DT1_061715',6;  % 12 Hz
        'fidel_CCDT_DT1_061715',8;  % 25 Hz
        'fidel_CCDT_DT1_061715',9;  % 50 Hz
        'fidel_CCDT_DT1_061715',10};% 100 Hz   
prewin = [-305 -5]; % pre-stim window (ms)
pstwin = [5 305]; % post-stim window (ms)

% compute band power pre and post stimulus
N = size(tkbk,1);
m = zeros(N,1);
v = zeros(N,1);
psig = zeros(N,1);
fstm = zeros(N,1);
for ii = 1:N
    p.tank = tkbk{ii,1}; % tank
    p.block = tkbk{ii,2}; % block
    p.twin = prewin; % time window around train (ms)
    p.bnm = 'alphapre'; % band name
    ram_pow(p);
    p.twin = pstwin; % time window around train (ms)
    p.bnm = 'alphapst'; % band name
    ram_pow(p);
    load([p.ddir p.subj '\processed\' tkbk{ii,1} '_' num2str(tkbk{ii,2}) '.mat'],'-mat');
    m(ii) = median(alphapst-alphapre);
    v(ii) = 1.96*mad(alphapst-alphapre,1)/sqrt(50);
    psig(ii) = ranksum(alphapre,alphapst);
    fstm(ii) = sfreq;
    clear alphapre alphapst
end    

% plot
figure('Name','fidel 0713','NumberTitle','off','Units','normalized','Position',[1/4 1/8 1/2 3/4],'Color','w');
bar(1:N,m,'k'); hold on;
errorbar(1:N,m,v,'k.');
for ii = 1:N
    if psig(ii)<0.05
        text(ii,m(ii)+v(ii),'\bf^*','FontSize',14,'HorizontalAlignment','center','VerticalAlignment','bottom');
    end
end
set(gca,'Box','off','XTick',1:N,'XTickLabel',fstm); axis square;
xlabel('stimulus train frequency (Hz)');
ylabel(['change in ' num2str(p.band(1)) '-' num2str(p.band(2)) ' Hz power (\muV^2)']);