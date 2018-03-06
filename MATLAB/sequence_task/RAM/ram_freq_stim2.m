function ram_freq_stim2
% function ram_freq_stim2
%   Compute power spectrum from stimulation files: compare spectrum pre and
%   post stimulus
% 
%   DR 02/2015

% parameters
ddir = 'D:\oblio\'; % directory
dtank = 'oblio_stim_DT1_021715'; % data tank
dblock = 4; % data block
ch = 7; % channel
twin1 = [-1024 -25]; % time window relative to stimulus (ms) - pre-stim
twin2 = [1251 2250]; % time window relative to stimulus (ms) - post-stim
frq = 1:300; % frequencies (Hz)
rln = 1; % remove line noise? (0 or 1)

% constants
fs = 24414.06; % sampling rate (Hz)

% stim triggers (ch25)
cd([ddir dtank '\Block-' num2str(dblock)]);
fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch25.sev'],'r');
fseek(fid,40,'bof');
dat = fread(fid,inf,'single');
fclose(fid);
itrain = find(dat(1:end-1)<0.5 & dat(2:end)>0.5); % trigger for each pulse/burst
itrain = itrain + round(0.7/1000*fs); % adjust trigger time for ~0.7-ms delay until stim voltage starts
swin1 = round(twin1(1)/1000*fs):round(twin1(2)/1000*fs);
indwin1 = itrain*ones(1,length(swin1)) + ones(length(itrain),1)*swin1;
indwin1(indwin1(:,1)<1,:) = [];
indwin1(indwin1(:,end)>length(dat),:) = [];
N = size(indwin1,1);
swin2 = round(twin2(1)/1000*fs):round(twin2(2)/1000*fs);
indwin2 = itrain*ones(1,length(swin2)) + ones(length(itrain),1)*swin2;
indwin2(indwin2(:,1)<1,:) = [];
indwin2(indwin2(:,end)>length(dat),:) = [];

% data
cd([ddir dtank '\Block-' num2str(dblock)]);
[blp,alp] = butter(4,500/(fs/2)); 
fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch' num2str(ch) '.sev'],'r');
fseek(fid,40,'bof');
dat = fread(fid,inf,'single');
fclose(fid);
dat = filtfilt(blp,alp,dat); % anti-aliasing filter (for downsampling)
datwin1 = dat(indwin1);
datwin1 = downsample(datwin1',10)'; % downsample
datwin2 = dat(indwin2);
datwin2 = downsample(datwin2',10)';
fs = round(fs/10);

% spectral analysis
figure('Name','oblio','NumberTitle','off','Units','normalized','Position',[1/3 1/4 1/3 1/2],'Color','w');
clrm = [0 0 1; 1 0 0]; clrv = [0.8 0.8 1; 1 0.8 0.8];
m = zeros(2,length(frq)); v = zeros(2,length(frq));
for jj = 1:2
    P = zeros(N,length(frq));
    for ii = 1:N
        eval(['x = datwin' num2str(jj) '(ii,:);']);
        x = x - median(x);
        if rln
            ln = mtmlinenoise(x,3,length(x),fs,60:60:frq(end));
            x = x - ln';
        end
        P(ii,:) = pmtm(x,3,frq,fs);
    end
    m(jj,:) = mean(P); v(jj,:) = 1.96*std(P)/sqrt(N);
    patch([frq, fliplr(frq)],[m(jj,:)-v(jj,:), fliplr(m(jj,:)+v(jj,:))],'k','FaceColor',clrv(jj,:),'EdgeColor','none'); hold on;
end
for jj = 1:2
    h(jj) = plot(frq,m(jj,:),'Color',clrm(jj,:),'LineWidth',2); hold on;
end
set(gca,'Box','off','XScale','log','YScale','log','XLim',[frq(1) frq(end)]); axis tight;
xlabel('Hz'); ylabel('\muV^2/Hz'); axis square;
legend(h,'pre-stim','post-stim'); legend('boxoff');