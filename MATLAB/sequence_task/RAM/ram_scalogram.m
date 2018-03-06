function ram_scalogram(varargin)
% function ram_scalogram(varargin)
%   Scalogram average across trials.
%
%   DR 03/2015

% parameters
ddir = 'D:\'; % directory
subj = 'ditka'; % subject
dtank = 'ditka_stim_DT1_031115'; % tank
dblock = 11; % block
chrec = 3; % rec channel(s)
twin = [-750 0]; % time window around train (ms)
dns = 10; % downsample factor
rln = 0; % remove line noise? (0 or 1)
clus = 1; % plot EP clusters? (0 or 1)
freq = 1:250; % scalogram frequencies (Hz)
wavl = 'morl'; % wavelet
if nargin
    v2struct(varargin{1});
end

% window
try load([ddir subj '\processed\' dtank '_' num2str(dblock) '.mat'],'-mat');
catch disp('run ''ram_stimtrig'' first'); return; end
swin = round(twin(1)/1000*fs):round(twin(2)/1000*fs);
indwin = itrain*ones(1,length(swin)) + ones(length(itrain),1)*swin;
indwin(indwin(:,1)<1,:) = [];
t = swin/fs*1000;
N = size(indwin,1);

% data
cd([ddir subj '\' dtank '\Block-' num2str(dblock)]);
fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch' num2str(chrec) '.sev'],'r');
fseek(fid,40,'bof');
dat = fread(fid,inf,'single'); % uV
fclose(fid);
if rln
    ln = mtmlinenoise(dat,3,round(fs),round(fs),[60 120 180]);
    dat = dat - ln;
end
if dns>1
    [b,a] = butter(2,round(fs/(2.5*dns))/(fs/2)); % antialiasing filter
    dat = filtfilt(b,a,dat);
end
datwin = dat(indwin);
datwin = downsample(datwin',dns)';
t = downsample(t,dns);
fs = round(fs/dns);

% wavelet transform
scales = centfrq(wavl)./(freq*(1/fs));
Pxx = zeros(N,length(freq),length(t));
for ii = 1:N
    x = datwin(ii,:);
    x = x-mean(x);
    coefs = cwt(x,scales,wavl);
    Pxx(ii,:,:) = reshape(coefs,1,length(freq),length(t));
end

% plot
if ~exist('epgroup') || ~clus, epgroup = ones(N,1); end
cln = length(unique(epgroup));
M = length(t);
P = zeros(length(freq),cln*M);
for ii = 1:cln
    ind = find(epgroup==ii);
    P(:,(ii-1)*M+1:ii*M) = squeeze(trimmean(abs(Pxx(ind,:,:)),5,1));
end
P = zscore(P')';
for ii = 1:cln
    figure('Name',dtank,'NumberTitle','off','Units','normalized','Position',[1/3 1/4 1/3 1/2],'Color','w');
    Pc = P(:,(ii-1)*M+1:ii*M);
    Pc = smooth2(Pc,5,2);
    imagesc(t,freq,Pc,[-3 3]);
    axis_square(gca); set(gca,'Box','off','YDir','normal','XLim',[t(1)+50 t(end)-50]); xlabel('time (s)'); ylabel('frequency (Hz)'); % remove edge effects
    text(max(get(gca,'XLim')),max(get(gca,'YLim')),['N = ' num2str(length(find(epgroup==ii)))],'HorizontalAlignment','right','VerticalAlignment','bottom');
    pos = get(gca,'Position'); hcb = colorbar('East'); cbpos = get(hcb,'Position');
    set(hcb,'Position',[pos(1)+pos(3)+cbpos(3)/2, cbpos(2:4)],'YAxisLocation','right');
    ylabel(hcb,'power (z-score)');
end
