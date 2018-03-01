function ram_stx_trials(varargin)
% function ram_stx_trials(varargin)
%   Stacked plot of every stimulus-triggered trial.
%
%   DR 02/2015

% parameters
ddir = 'D:\'; % directory
subj = 'ditka'; % subject
tank = 'ditka_stim_DT1_033015'; % tank
block = 4; % block
chrec = 6; % rec channel
twin = [-50 300]; % time window around train (ms)
scale = 200; % scale (uV)
bpf = [1 1000]; % bandpass filter (Hz)
dns = 10; % downsample factor
rln = 0; % remove line noise? (0 or 1)
clus = 1; % show EP clusters? (0 or 1)
if nargin
    v2struct(varargin{1});
end

% window
try load([ddir subj '\processed\' tank '_' num2str(block) '.mat'],'-mat');
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
indwin = downsample(indwin',dns)'; % downsample
t = downsample(swin/fs*1000,dns);
N = size(indwin,1);

% data
cd([ddir subj '\' tank '\Block-' num2str(block)]);
fid = fopen([tank '_Block-' num2str(block) '_xRAW_ch' num2str(chrec) '.sev'],'r');
fseek(fid,40,'bof');
dat = fread(fid,inf,'single'); % uV
fclose(fid);
if rln
    ln = mtmlinenoise(dat,3,round(fs),round(fs),[60 120]);
    dat = dat - ln;
end
[b,a] = butter(2,bpf/(fs/2)); % bandpass filter
dat = filter(b,a,dat); % causal so post-stim effects don't change pre-stim baseline
datwin = dat(indwin);

% plot
figure('Name',tank,'NumberTitle','off','Units','normalized','Position',[1/3 1/8 1/3 3/4],'Color','w');
clr = ['k','r','b','g']; % assumes no more than this many clusters will be used
if ~exist('epgroup','var') || ~clus, epgroup = ones(N,1); end
for ii = 1:N
    plot(t,datwin(ii,:)/scale+ii,'Color',clr(epgroup(ii)),'Tag',num2str(ii)); hold on;
%     plot(t,datwin(ii,:),'Color',clr(epgroup(ii)),'Tag',num2str(ii)); hold on;
end
set(gca,'Box','off','XLim',twin,'YLim',[0 N(1)+1]); xlabel('time (ms)'); ylabel('stimulus number');
line(twin(end)*ones(1,2),[0 10],'Color','k');
text(twin(end),5,[num2str(scale*10) '\muV'],'HorizontalAlignment','left','VerticalAlignment','middle');
