function ram_pow(varargin)
% function ram_pow(varargin)
%   Compute power in frequency band and time window relative to stimulus.
%
%   DR 03/2015

% parameters
ddir = 'E:\'; % directory
subj = 'fidel'; % subject
tank = 'fidel_CCDT_DT1_061715'; % tank
block = 8; % block
chrec = 3; % rec channel
twin = [10 300]; % time window around train (ms)
band = [8 15]; % frequency band (Hz)
bnm = 'alphapst'; % band name
dns = 10; % downsample factor
rln = 0; % remove line noise? (0 or 1)
if nargin
    v2struct(varargin{1});
end

% window
try 
    load([ddir subj '\processed\' tank '_' num2str(block) '.mat'],'-mat');
    if spuls>1 && twin(1)>0, itrain = itrain + round((spuls-1)/sfreq*fs); end % last pulse of each train (sample)
catch disp('run ''ram_stimtrig'' first'); return;
end
swin = round(twin(1)/1000*fs):round(twin(2)/1000*fs);
indwin = itrain*ones(1,length(swin)) + ones(length(itrain),1)*swin;
indwin(indwin(:,1)<1,:) = [];
indwin = downsample(indwin',dns)'; % downsample
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
datwin = dat(indwin);
fs = round(fs/dns);

% band power
p = zeros(N,1);
for ii = 1:N
    p(ii) = bandpower(datwin(ii,:),fs,band); 
end
eval([bnm ' = p;']);

% save
save([ddir subj '\processed\' tank '_' num2str(block) '.mat'],bnm,'-append','-mat');