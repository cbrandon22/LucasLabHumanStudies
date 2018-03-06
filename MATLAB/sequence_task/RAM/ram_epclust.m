function ram_epclust(varargin)
% function ram_epclust(varargin)
%   Cluster evoked potentials into N groups.
%
%   DR 03/2015

% parameters
ddir = 'D:\'; % directory
subj = 'oblio'; % subject
tank = 'oblio_stim_DT1_040115'; % tank
block = 13; % block
chrec = 7; % rec channel
bpf = [1 1000]; % evoked-potential bandpass (Hz)
twin = [4 20]; % evoked-potential window (ms)
rln = 0; % remove line noise? (0 or 1)
cln = 2; % number of clusters
if nargin
    v2struct(varargin{1});
end

% window
try load([ddir subj '\processed\' tank '_' num2str(block) '.mat'],'-mat');
catch disp('run ''ram_stimtrig'' first'); return; end
swin = round(twin(1)/1000*fs):round(twin(2)/1000*fs);
indwin = itrain*ones(1,length(swin)) + ones(length(itrain),1)*swin;
indwin(indwin(:,1)<1,:) = [];

% EPs
cd([ddir subj '\' tank '\Block-' num2str(block)]);
fid = fopen([tank '_Block-' num2str(block) '_xRAW_ch' num2str(chrec) '.sev'],'r');
fseek(fid,40,'bof');
dat = fread(fid,inf,'single'); % uV
fclose(fid);
if rln
    ln = mtmlinenoise(dat,3,round(fs),round(fs),[60 120]);
    dat = dat - ln;
end
[b,a] = butter(2,bpf/(fs/2));
dat = filtfilt(b,a,dat);
datwin = dat(indwin);

% cluster
epgroup = kmeans(datwin,cln);
mpp = zeros(cln,1); inds = cell(cln,1);
for ii = 1:cln % re-name 'epgroup' in order of ascending mean EP size
    inds{ii} = find(epgroup==ii);
    m = mean(datwin(inds{ii},:));
    mpp(ii) = max(m)-min(m);
end
[~,isort] = sort(mpp);
for ii = 1:cln
    epgroup(inds{ii}) = isort(ii);
end

% save
save([ddir subj '\processed\' tank '_' num2str(block) '.mat'],'epgroup','-append','-mat');