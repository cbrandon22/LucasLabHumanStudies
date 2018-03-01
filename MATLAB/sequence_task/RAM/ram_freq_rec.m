function ram_freq_rec
% function oblio_freq_rec
%   Compute power spectrum from recording files (no stim).
% 
%   DR 01/2015

% parameters
ddir = 'G:\ditka\'; % directory
dtank = 'ditka_rec_DT1_033015'; % data tank
dblock = 1; % data block
ch = [3:5 10:14]; % channels
frq = [0.3 300]; % frequency range (Hz)
rln = 0; % remove line noise? (0 or 1)
car = 0; % common average reference? (0 or 1)

% constants
fs = 24414.06; % sampling rate (Hz)

% data
cd([ddir dtank '\Block-' num2str(dblock)]);
[blp,alp] = butter(4,500/(fs/2)); 
for ich = 1:length(ch)
    fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch' num2str(ch(ich)) '.sev'],'r');
    fseek(fid,40,'bof');
    dat = fread(fid,inf,'single');
    fclose(fid);
    dat = filtfilt(blp,alp,dat); % anti-aliasing filter
    dat = downsample(dat,10); % downsample
    if ich == 1, DAT = zeros(length(dat),length(ch)); end
    DAT(:,ich) = dat;
end
fs = round(fs/10);
clear dat;

% condition
if rln
    for ich = 1:length(ch)
        ln = mtmlinenoise(DAT(:,ich),3,fs,fs,60:60:300);
        DAT(:,ich) = DAT(:,ich) - ln;
    end
end
if car
    DAT = DAT - mean(DAT,2)*ones(1,size(DAT,2));
end

% plot
figure('Name',dtank,'NumberTitle','off','Units','normalized','Position',[1/3 1/4 1/3 1/2],'Color','w');
if length(ch)>1, clr = colormap(jet(length(ch))); else clr = 'k'; end
for ich = 1:length(ch)
    dat = DAT(:,ich);
    dat = dat-median(dat);
%     [P,f] = pwelch(dat,fix(fs),0,fix(fs),fs); % 1-sec windows, no overlap
    [P,f] = pwelch(dat,10*fix(fs),5*fix(fs),30*fix(fs),fs); % 10-sec windows, 50% overlap
    ind = find(f>=frq(1) & f<=frq(2));
    f = f(ind); P = P(ind);
    plot(f,P,'Color',clr(ich,:),'Tag',num2str(ch(ich))); hold on;
end
set(gca,'Box','off','XScale','log','YScale','log','XLim',frq); axis tight;
xlabel('Hz'); ylabel('\muV^2/Hz'); axis square;
