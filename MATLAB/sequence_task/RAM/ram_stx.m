function ram_stx
% function ram_stx
%   Stacked plot of every stimulus-triggered trial.
%
%   DR 02/2015

% parameters
ddir = 'D:\oblio\'; % data directory
dtank = 'oblio_stim_DT1_030215'; % data tank
dblock = 20; % data block
chrec = 7; % rec channel
twin = [-25 200]; % time window around pulse/burst (ms)
scale = 50; % uV
bpf = [10 500]; % bandpass filter (Hz)
dns = 10; % downsample factor
rln = 0; % remove line noise? (0 or 1)

% constants
fs = 24414.06; % sampling rate (Hz)

% stim triggers (ch25)
cd([ddir dtank '\Block-' num2str(dblock)]);
fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch25.sev'],'r');
fseek(fid,40,'bof');
dat = fread(fid,inf,'single');
fclose(fid);
itrain = find(dat(1:end-1)<0.5 & dat(2:end)>0.5); % trigger for each pulse/burst
itrain = itrain + round(0.3/1000*fs); % adjust trigger time for ~0.3-ms delay until stim voltage starts
swin = round(twin(1)/1000*fs):round(twin(2)/1000*fs);
indwin = itrain*ones(1,length(swin)) + ones(length(itrain),1)*swin;
indwin(indwin(:,1)<1,:) = [];
indwin(indwin(:,end)>length(dat),:) = [];
t = swin/fs*1000;
N = size(indwin,1);

% rec channel
fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch' num2str(chrec) '.sev'],'r');
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
datwin = downsample(datwin',dns)'; % downsample
t = downsample(t,dns);

% trial plot
figure('Name',dtank,'NumberTitle','off','Units','normalized','Position',[1/3 1/8 1/3 3/4],'Color','w');
for ii = 1:N
    plot(t,datwin(ii,:)/scale+ii,'k','Tag',num2str(ii)); hold on;
end
set(gca,'Box','off','XLim',twin,'YLim',[0 N(1)+1]); xlabel('time (ms)'); ylabel('stimulus number');
line(twin(end)*ones(1,2),[0 10],'Color','k');
text(twin(end),5,[num2str(scale*10) '\muV'],'HorizontalAlignment','left','VerticalAlignment','middle');
