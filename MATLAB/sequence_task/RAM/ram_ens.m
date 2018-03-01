function ram_ens
% function ram_ens
%   Analyze files from Medtronic ENS system.
%   09/2015 system test:
%   row =  1      2-9,   10-13, 14,  15,  16,  17
%   dat = [samp#, HC8-1, EC8-5, EC3, EC1, EC4, EC2] 
%
%   DR 09/2015

% parameters
ddir = 'E:\fidel\fidel_ENS_091515\'; % data directory
dtank = 'Fidel_OL_EC12.txt'; % data tank
twin = [-300 -10]; % time window around pulse/burst (ms)
fs = 1000; % sampling rate (Hz)

% load data
cd(ddir);
dat = csvread(dtank);

% find stim triggers
thresh = -1e4;
chthr = 10;
itrain = find(dat(1:end-1,chthr)>thresh & dat(2:end,chthr)<thresh);
itrain(itrain<10*fs) = []; % remove any stimuli within first 10-sec (startup transient)
swin = round(twin(1)/1000*fs):round(twin(2)/1000*fs);
indwin = itrain*ones(1,length(swin)) + ones(length(itrain),1)*swin;
indwin(indwin(:,1)<1,:) = [];
indwin(indwin(:,end)>length(dat),:) = [];
t = swin/fs*1000;
N = size(indwin,1);

% convert to uV
dat = dat*0.2; % 0.217

% % sta plot
% figure('Name',dtank,'NumberTitle','off','Units','normalized','Position',[1/3 1/8 1/3 3/4],'Color','w');
% clr = colormap(jet(7)); ind = 1;
% for ii = 8:-1:2 % HC1-8    
% % for ii = [15 17 14 16 13:-1:10] % EC1-8
%     cdat = dat(:,ii);
%     datwin = cdat(indwin);
%     datwin = datwin - median(datwin(:,t<-1),2)*ones(1,length(t)); % remove dc offset
%     plot(t,median(datwin),'Color',clr(ind,:)); hold on;
%     ind = ind+1;
% end
% axis square; set(gca,'Box','off','XLim',twin);

% % sta + ci plot
% figure('Name',dtank,'NumberTitle','off','Units','normalized','Position',[1/4 1/8 1/2 3/4],'Color','w');
% chrec = 6; % rec channel
% cdat = dat(:,chrec);
% datwin = cdat(indwin);
% datwin = datwin - median(datwin(:,t<-1),2)*ones(1,length(t)); % remove dc offset
% m = median(datwin);
% v = 1.96*std(datwin)/sqrt(N)/1.6;
% patch([t fliplr(t)],[m+v fliplr(m-v)],'k','FaceColor',[0.8 0.8 0.8],'EdgeColor','none'); hold on;
% plot(t,m,'Color','k');
% axis square; set(gca,'Box','off','XLim',twin);

% % stx plot
% chrec = 6; % rec channel
% cdat = dat(:,chrec);
% datwin = cdat(indwin);
% scale = 200; % uV
% figure('Name',dtank,'NumberTitle','off','Units','normalized','Position',[1/3 1/8 1/3 3/4],'Color','w');
% for ii = 1:N
%     plot(t,datwin(ii,:)/scale+ii,'k','Tag',num2str(ii)); hold on;
% end
% set(gca,'Box','off','XLim',twin,'YLim',[0 N(1)+1]); xlabel('time (ms)'); ylabel('stimulus number');
% line(twin(end)*ones(1,2),[0 10],'Color','k');
% text(twin(end),5,[num2str(scale*10) '\muV'],'HorizontalAlignment','left','VerticalAlignment','middle');

% freq plot
chrec = 6; % rec channel
frq = 5:100; % frequencies (Hz)
cdat = dat(:,chrec);
datwin = cdat(indwin);
bp = zeros(N,1);
for ii = 1:N
    bp(ii) = bandpower(datwin(ii,:),fs,[6 12]); 
end
P = zeros(N,length(frq));
for ii = 1:N
    x = datwin(ii,:);
    x = x - median(x);
    P(ii,:) = pmtm(x,3,frq,fs);
end
cln = 2;
epgroup = ones(N,1);
epgroup([11 15 18 22 24 27 31 33 34 39 40 43 48 51 55 60]) = 2;
m = zeros(cln,length(frq));
v = zeros(cln,length(frq));
for ii = 1:cln
    ind = find(epgroup==ii);
    m(ii,:) = mean(P(ind,:));
    v(ii,:) = 1.96*std(P(ind,:))/sqrt(length(ind));
    v(ii,frq<20) = v(ii,frq<20)/1.6;
end
figure('Name',dtank,'NumberTitle','off','Units','normalized','Position',[1/3 1/8 1/3 3/4],'Color','w');
clr = ['k','g'];
for ii = 1:cln
    patch([frq,fliplr(frq)],[m(ii,:)+v(ii,:),fliplr(m(ii,:)-v(ii,:))],'k','FaceColor',[0.8 0.8 0.8],'EdgeColor','none'); hold on;
end
for ii = 1:cln
    plot(frq,m(ii,:),'Color',clr(ii),'Tag',num2str(ii));
end
set(gca,'Box','off','XScale','log','YScale','log','XLim',[frq(1) frq(end)]); axis tight;
xlabel('Hz'); ylabel('\muV^2/Hz'); title([num2str(twin(1)) ' to ' num2str(twin(2)) ' ms']); axis square;
[h,p,ci,stats] = ttest2(bp(find(epgroup==1)),bp(find(epgroup==2)));
disp(['t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(p)]);
