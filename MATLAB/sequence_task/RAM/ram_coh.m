function ram_coh
% function ram_coh
%   Compute coherence between pairs of recording sites (no stim).
% 
%   DR 03/2015

% parameters
ddir = 'D:\ditka\'; % directory
dtank = 'ditka_rec_DT1_033015'; % data tank
dblock = 1; % data block
ch = [4 12]; % channel(s)
freq = 1:150; % frequencies (Hz)
dns = 10; % downsample factor

% data
fs = 24414.06; % sampling rate (Hz)
[bA,aA] = butter(2,500/(fs/2)); % antialiasing
Nch = length(ch);
pairs = nchoosek(1:Nch,2);
COH = zeros(length(freq),size(pairs,1));
for jj = 1:size(pairs,1)
    cd([ddir dtank '\Block-' num2str(dblock)]);
    fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch' num2str(ch(pairs(jj,1))) '.sev'],'r');
    fseek(fid,40,'bof');
    dat1 = fread(fid,inf,'single');
    fclose(fid);
    dat1 = filtfilt(bA,aA,dat1);
    dat1 = downsample(dat1',dns)';
    fid = fopen([dtank '_Block-' num2str(dblock) '_xRAW_ch' num2str(ch(pairs(jj,2))) '.sev'],'r');
    fseek(fid,40,'bof');
    dat2 = fread(fid,inf,'single');
    fclose(fid);
    dat2 = filtfilt(bA,aA,dat2);
    dat2 = downsample(dat2',dns)';
    dfs = round(fs/dns);
    
    % coherence
    [Cxy,f] = mscohere(dat1,dat2,dfs,0,dfs,dfs); % 1-sec windows, no overlap
    ind = find(f>=freq(1) & f<=freq(end));
    COH(:,jj) = Cxy(ind);
end
figure('Name',dtank,'NumberTitle','off','Units','normalized','Position',[1/4 1/8 1/2 3/4],'Color','w');
for jj = 1:size(pairs,1)
    plot(freq,COH(:,jj),'k','Tag',[num2str(pairs(jj,1)) '-' num2str(pairs(jj,2))]); hold on;
end
axis square; set(gca,'Box','off','XLim',[freq(1) freq(end)],'XScale','log','YLim',[0 1]);
xlabel('frequency (Hz)'); ylabel('mean-square coherence');