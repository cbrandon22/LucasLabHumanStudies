function ram_freq_stim(varargin)
% function ram_freq_stim(varargin)
%   Compute power spectrum from stimulation files.
% 
%   DR 02/2015

% parameters
ddir = 'G:\'; % directory
subj = 'ditka'; % subject
dtank = 'ditka_stim_DT1_032515'; % tank
dblock = 2; % block
ch = 4; % channels
twin = [-1000 -5]; % time window relative to stimulus (ms)
frq = 1:300; % frequencies (Hz)
dns = 10; % downsample factor
car = 0; % common average reference? (0 or 1)
rln = 0; % remove line noise? (0 or 1)
clus = 1; % plot EP clusters? (0 or 1)
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
N = size(indwin,1);

% data
cd([ddir subj '\' dtank '\Block-' num2str(dblock)]);
[blp,alp] = butter(4,500/(fs/2)); 
for ich = 1:length(ch)
    fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch' num2str(ch(ich)) '.sev'],'r');
    fseek(fid,40,'bof');
    dat = fread(fid,inf,'single');
    fclose(fid);
    dat = filtfilt(blp,alp,dat); % anti-aliasing filter
    if ich == 1, DAT = zeros(length(dat),length(ch)); end
    DAT(:,ich) = dat;
end
clear dat;
if car && size(DAT,2)>1
    DAT = DAT - mean(DAT,2)*ones(1,size(DAT,2));
end

% spectral analysis
figure('Name',dtank,'NumberTitle','off','Units','normalized','Position',[1/3 1/4 1/3 1/2],'Color','w');
if length(ch)>1 || ~clus || ~exist('epgroup','var')
    clr = colormap(lines(length(ch)));
    fs = round(fs/10);
    for ich = 1:length(ch)
        dat = DAT(:,ich);
        datwin = dat(indwin);
        P = zeros(N,length(frq));
        for ii = 1:N
            x = datwin(ii,:);
            x = x - median(x);
            if rln
                ln = mtmlinenoise(x,3,length(x),fs,60:60:frq(end));
                x = x - ln';
            end
            P(ii,:) = pmtm(x,3,frq,fs);
        end
        plot(frq,median(P),'Color',clr(ich,:),'Tag',num2str(ch(ich))); hold on;
    end
else
    fs = round(fs/10);
    dat = DAT(:,ich);
    datwin = dat(indwin);
    P = zeros(N,length(frq));
    for ii = 1:N
        x = datwin(ii,:);
        x = x - median(x);
        if rln
            ln = mtmlinenoise(x,3,length(x),fs,60:60:frq(end));
            x = x - ln';
        end
        P(ii,:) = pmtm(x,3,frq,fs);
    end
    cln = length(unique(epgroup));
    m = zeros(cln,length(frq));
    v = zeros(cln,length(frq));
    for ii = 1:cln
        ind = find(epgroup==ii);
        m(ii,:) = mean(P(ind,:));
        v(ii,:) = 1.96*std(P(ind,:))/sqrt(length(ind));
    end
    clr = colormap(lines(cln));
    for ii = 1:cln
        patch([frq,fliplr(frq)],[m(ii,:)+v(ii,:),fliplr(m(ii,:)-v(ii,:))],'k','FaceColor',[0.8 0.8 0.8],'EdgeColor','none'); hold on;
    end
    for ii = 1:cln
        plot(frq,m(ii,:),'Color',clr(ii,:),'Tag',num2str(ii));
    end
end
set(gca,'Box','off','XScale','log','YScale','log','XLim',[frq(1) frq(end)]); axis tight;
xlabel('Hz'); ylabel('\muV^2/Hz'); title([num2str(twin(1)) ' to ' num2str(twin(2)) ' ms']); axis square;

% % spectral amplitude (computed directly from FFT)
% M = size(indwin,2);
% A = fft(datwin',M)';
% A = abs(A(:,1:floor(M/2)+1))/M;
% f = fs/M*(0:floor(M/2));
% ind = find(f>freq(1) & f<freq(2));
% f = f(ind);
% A = A(:,ind);