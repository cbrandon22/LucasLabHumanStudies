function ram_cta
% function ram_cta
%   Compute field potential cycle-triggered averages and mean phase
%   differences for a specific frequency band from recording files.
% 
%   DR 03/2015

% parameters
ddir = 'G:\ditka\'; % directory
dtank = 'ditka_rec_DT1_032715'; % data tank
dblock = 1; % data block
ch = 1:16; % channel(s)
chref = 4; % reference channel
band = [2 6]; % frequency band (Hz)
dns = 10; % downsample factor
twin = [-1000/band(1) 1000/band(1)]; % time window (ms)
hamp = 1; % include all cycles (0) or only those with high ref amplitude (1)?
scale = 30; % uV

% cycle peaks of reference channel
cd([ddir dtank '\Block-' num2str(dblock)]);
fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch' num2str(chref) '.sev'],'r');
fseek(fid,40,'bof');
dat = fread(fid,inf,'single');
fclose(fid);
fs = 24414.06; % sampling rate (Hz)
[bA,aA] = butter(2,500/(fs/2)); % antialiasing
dat = filtfilt(bA,aA,dat);
dat = downsample(dat',dns)';
fs = round(fs/dns);
[bB,aB] = butter(2,band/(fs/2)); % bandpass
an = hilbert(filtfilt(bB,aB,dat));
A = abs(an);
Pref = angle(an);
Pu = unwrap(Pref);
ipeaks = round(interp1(Pu,1:length(Pref),0:2*pi:floor(Pu(end)/(2*pi))*2*pi))';
ipeaks(isnan(ipeaks)) = [];
if hamp
    Ap = A(ipeaks);
    ipeaks(Ap<median(Ap)) = [];
end
swin = round(twin(1)/1000*fs):round(twin(2)/1000*fs);
indwin = ipeaks*ones(1,length(swin)) + ones(length(ipeaks),1)*swin;
indwin(indwin(:,1)<1,:) = [];
indwin(indwin(:,end)>length(Pref),:) = [];
t = swin/fs*1000;

% cycle triggered averages and mean phase differences
CTA = zeros(length(ch),size(indwin,2));
MPD = zeros(length(ch),1);
for jj = 1:length(ch)
    fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch' num2str(ch(jj)) '.sev'],'r');
    fseek(fid,40,'bof');
    dat = fread(fid,inf,'single');
    fclose(fid);
    dat = filtfilt(bA,aA,dat);
    dat = downsample(dat',dns)';
    dat = filtfilt(bB,aB,dat);
    datwin = dat(indwin);
    CTA(jj,:) = mean(datwin);
    P = angle(hilbert(dat));
    dtheta = angle(exp(i*(P-Pref)));
    c = nansum(cos(dtheta));
    s = nansum(sin(dtheta));
    MPD(jj) = atan2(s,c)*180/pi;
end

% plot
figure('Name',dtank,'NumberTitle','off','Units','normalized','Position',[1/3 1/8 1/3 3/4],'Color','w');
for jj = 1:length(ch)
    hl = plot(t,CTA(jj,:)/scale+jj,'k','LineWidth',1,'Tag',num2str(jj)); hold on;
    if ch(jj) == chref, set(hl,'LineWidth',2); end
end
set(gca,'Box','off','XLim',twin,'YLim',[0 length(ch)+1],'YTick',1:length(ch));
xlabel('ms'); ylabel('channel'); title([num2str(band(1)) ' - ' num2str(band(2)) ' Hz']);
line([0 0],get(gca,'YLim'),'Color','k')
text(max(get(gca,'XLim')),max(get(gca,'YLim')),['N = ' num2str(length(ipeaks))],'HorizontalAlignment','right','VerticalAlignment','bottom');
line(max(get(gca,'XLim'))*ones(1,2),[0 1],'Color','k');
text(max(get(gca,'XLim')),0.5,[num2str(scale) '\muV'],'HorizontalAlignment','left','VerticalAlignment','middle');
figure('Name',dtank,'NumberTitle','off','Units','normalized','Position',[1/4 1/4 1/2 1/2],'Color','w');
plot(ch,MPD,'ko','MarkerFaceColor','k');
set(gca,'Box','off','YLim',[-180 180]);
xlabel('channel'); ylabel('mean phase difference (deg)');