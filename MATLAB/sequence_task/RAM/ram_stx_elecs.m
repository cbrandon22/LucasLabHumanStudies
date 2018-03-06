function ram_stx_elecs(varargin)
% function ram_stx_elecs(varargin)
%   Stacked plots of each stimulus-triggered trial across electrodes.
%
%   DR 03/2015

% parameters
ddir = 'D:\'; % directory
subj = 'ditka'; % subject
dtank = 'ditka_stim_DT1_032715'; % tank
dblock = 2; % block
ch = 1:8; % channels
twin = [-500 750]; % time window around train (ms)
scale = 500; % scale (uV)
bpf = [1 1000]; % bandpass filter (Hz)
dns = 10; % downsample factor
if nargin
    v2struct(varargin{1});
end

% window
try load([ddir subj '\processed\' dtank '_' num2str(dblock) '.mat'],'-mat');
catch
    cd([ddir subj '\' dtank '\Block-' num2str(dblock)]);
    fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch25.sev'],'r');
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
cd([ddir subj '\' dtank '\Block-' num2str(dblock)]);
[b,a] = butter(2,bpf/(fs/2)); % bandpass filter
datwin = zeros(N,size(indwin,2),length(ch));
for ich = 1:length(ch)
    fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch' num2str(ch(ich)) '.sev'],'r');
    fseek(fid,40,'bof');
    dat = fread(fid,inf,'single');
    fclose(fid);
    dat = filter(b,a,dat); % causal so post-stim effects don't change pre-stim baseline
    datwin(:,:,ich) = dat(indwin);
end
clear dat;

% plot
figure('Name',dtank,'NumberTitle','off','Units','normalized','Position',[1/3 1/4 1/3 1/2],'Color','w');
for ii = 1:N
    cla;
    for ich = 1:length(ch)
        plot(t,squeeze(datwin(ii,:,ich))/scale+ich,'k'); hold on;
    end
    set(gca,'Box','off','XLim',twin,'YLim',[0 length(ch)+1],'YTick',1:length(ch)); xlabel('time (ms)'); ylabel('channel'); 
    title(['N = ' num2str(ii)]);
%     title(['N = ' num2str(ii) ', pLG = ' num2str(pLG(ii)) ', pHG = ' num2str(pHG(ii))]);
    pause
end
