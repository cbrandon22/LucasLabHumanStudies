function ram_sta(varargin)
% function ram_sta(varargin)
%   Stimulus-triggered average across trials.
%
%   DR 01/2015

% parameters
ddir = 'G:\'; % directory
subj = 'fidel'; % subject
tank = 'fidel_CCDT_DT1_060915'; % tank
block = 7; % block
chrec = 9:16; % rec channel(s)
chstim = 1; % stim channel(s)
twin = [-50 450]; % time window around train (ms)
dns = 10; % downsample factor
rln = 0; % remove line noise? (0 or 1)
clus = 0; % plot EP clusters? (0 or 1)
pci = 0; % plot 95% confidence intervals? (0 or 1)
if nargin
    v2struct(varargin{1});
end

% window
try 
    load([ddir subj '\processed\' tank '_' num2str(block) '.mat'],'-mat');
%     itrain = itrain(126:150);
%     itrain = itrain + round((spuls-1)/sfreq*fs); % last pulse of each train (sample)
catch
    cd([ddir subj '\' tank '\Block-' num2str(block)]);
    fid = fopen([tank '_Block-' num2str(block) '_xRAW_ch25.sev'],'r');
    fseek(fid,40,'bof');
    dat = fread(fid,inf,'single');
    fclose(fid);
    fs = 24414.06;
    itrain = find(dat(1:end-1)<0.5 & dat(2:end)>0.5); % trigger for each pulse train
    itrain = itrain + round(0.5/1000*fs);
end
swin = round(twin(1)/1000*fs):round(twin(2)/1000*fs);
indwin = itrain*ones(1,length(swin)) + ones(length(itrain),1)*swin;
indwin(indwin(:,1)<1,:) = [];
t = swin/fs*1000;

% stim channel plot
cd([ddir subj '\' tank '\Block-' num2str(block)]);
figure('Name',tank,'NumberTitle','off','Units','normalized','Position',[1/3 1/4 1/3 1/2],'Color','w');
% subplot(4,1,1), hold on; clr = [0 0 0; 0.7 0.7 0.7];
% for ii = 1:length(chstim)
%     fid = fopen([tank '_Block-' num2str(block) '_xRAW_ch' num2str(chstim(ii)+16) '.sev'],'r');
%     fseek(fid,40,'bof');
%     dat = fread(fid,inf,'single')/1e6; % uA/V
%     fclose(fid);
%     datwin = dat(indwin);
%     m = mean(datwin);
%     plot(t,m,'Color',clr(ii,:));
% end
% axis tight; set(gca,'Box','off','XLim',[t(1) t(end)]);
% text(max(get(gca,'XLim')),max(get(gca,'YLim')),['N = ' num2str(size(indwin,1))],'HorizontalAlignment','right','VerticalAlignment','top');
% ylabel('\muA or V');
 
% rec channel plot
indwin = downsample(indwin',dns)'; % downsample
t = downsample(swin/fs*1000,dns);
Ns = size(indwin,2);
if length(chrec)>1 | ~clus | ~exist('epgroup')
    m = zeros(length(chrec),Ns);
    v = zeros(length(chrec),Ns);
    for ii = 1:length(chrec)
        fid = fopen([tank '_Block-' num2str(block) '_xRAW_ch' num2str(chrec(ii)) '.sev'],'r');
        fseek(fid,40,'bof');
        dat = fread(fid,inf,'single'); % uV
        fclose(fid);
        if rln
            ln = mtmlinenoise(dat,3,round(fs),round(fs),60:60:180);
            dat = dat - ln;
        end
%         [b,a] = butter(2,5/(fs/2),'high');
%         dat = filter(b,a,dat);
        datwin = dat(indwin); 
        datwin = datwin - median(datwin(:,t<-1),2)*ones(1,Ns); % remove dc offset
        m(ii,:) = mean(datwin); 
        v(ii,:) = 1.96*std(datwin)/sqrt(size(datwin,1));
    end
    subplot(4,1,1:4), hold on; if length(chrec)>1, clr = colormap(jet(length(chrec))); else clr = [0 0 0]; end
    if pci
        for ii = 1:length(chrec)
            patch([t,fliplr(t)],[m(ii,:)+v(ii,:),fliplr(m(ii,:)-v(ii,:))],'k','FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
        end
    end
    for ii = 1:length(chrec)
        plot(t,m(ii,:),'Color',clr(ii,:),'Tag',num2str(chrec(ii)));
    end
else
    fid = fopen([tank '_Block-' num2str(block) '_xRAW_ch' num2str(chrec) '.sev'],'r');
    fseek(fid,40,'bof');
    dat = fread(fid,inf,'single'); % uV
    fclose(fid);
    if rln
        ln = mtmlinenoise(dat,3,round(fs),round(fs),60:60:180);
        dat = dat - ln;
    end
%     [b,a] = butter(2,1/(fs/2),'high');
%     dat = filter(b,a,dat);
    datwin = dat(indwin);
    datwin = datwin - median(datwin,2)*ones(1,Ns); % remove dc offset
    cln = length(unique(epgroup));
    m = zeros(cln,size(indwin,2));
    v = zeros(cln,size(indwin,2));
    for ii = 1:cln
        ind = find(epgroup==ii);
        m(ii,:) = mean(datwin(ind,:));
        v(ii,:) = 1.96*std(datwin(ind,:))/sqrt(length(ind));
    end
    subplot(4,1,1:4), hold on; clr = ['k','r','b','g']; % assumes no more than this many clusters will be used
    if pci
        for ii = 1:cln
            patch([t,fliplr(t)],[m(ii,:)+v(ii,:),fliplr(m(ii,:)-v(ii,:))],'k','FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
        end
    end
    for ii = 1:cln
        plot(t,m(ii,:),'Color',clr(ii),'Tag',num2str(ii));
    end
end
axis tight; set(gca,'Box','off','XLim',[t(1) t(end)],'YLim',[-500 500]);
xlabel('ms'); ylabel('\muV');
